// mendrive_core.cpp  (автономная версия: MenDriveCpp/standalone/)
//
// Копия ../mendrive_core.cpp (C++ порт резонатора МенДрайв: JA / LLG / Hybrid /
// Preisach / Preisach2D, СГС) с МИНИМАЛЬНЫМИ исправлениями. Численная схема и
// физика не изменены. Все отличия от оригинала помечены:
//   // FIX:  -- исправление ошибки (память, индексация, переполнение, UB);
//   // PERF: -- оптимизация БЕЗ изменения арифметики (те же значения);
//   // DIAG: -- новая диагностика (на динамику не влияет);
//   // NOTE: -- замечание к оригинальному коду, оставленному без изменений.
//
// Этот файл можно:
//   (а) включить в main.cpp (#include "mendrive_core.cpp") -- автономная программа;
//   (б) собрать как разделяемую библиотеку для старой ctypes-обёртки
//       (extern "C" API ниже совместим с mendrive_ctypes.py без изменений).
//
// Соответствие Python-реализации:
//   langevin_L        -> langevin_L()
//   JASmoothRK4Robust -> struct JAEngine   (численная dM/dH -- без сингулярности при H->0)
//   LLGVector         -> struct LLGEngine  (Rodrigues-rotation + демпфирование Гильберта)
//   LLGJAHybrid       -> связка внутри MenDriveSim::step_once (ferrite_model==2)
//   MenDriveFDTD      -> struct MenDriveSim

#include <vector>
#include <cmath>
#include <cstdint>
#include <cstdio>
#include <cstring>
#include <algorithm>
#include <memory>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <mutex>
#include <condition_variable>
#include <atomic>
#include <functional>
#ifdef USE_OPENMP
#include <omp.h>
#endif

enum {
    JA         = 0,
    LLG        = 1,
    Hybrid     = 2,
    Preisach   = 3,
    Preisach2D = 4,
} FERRITE_CODE;

// FIX: алгоритм хранения состояния Preisach2D.
//   PR2D_BITS -- бит-упакованное состояние (1 бит на пару вместо int8, память /8),
//                суммирование M в ТОМ ЖЕ порядке, что в оригинале -> те же числа;
//   PR2D_ROWS -- точное «строчное» представление: для каждой строки hc_i сетки
//                состояние реле вдоль hi_j всегда имеет вид (+1..+1,-1..-1)
//                (alpha_ij и beta_ij монотонны по j), т.е. задаётся одним целым
//                K_i. Память O(n_grid) на ячейку вместо O(n_grid^2), время шага
//                O(n_grid) вместо O(n_grid^2). Состояния реле совпадают точно;
//                M отличается от PR2D_BITS только округлением суммы (~1e-15).
enum { PR2D_BITS = 0, PR2D_ROWS = 1 };

static inline long double langevin_L(long double x) {
    if (std::fabs(x) < 1e-6) { return x/3.0 - x*x*x/45.0; }
    return 1.0/std::tanh(x) - 1.0/x;
}

// ----------------------------------------------------------------------------
// PERF: минимальный пул потоков на std::thread (без внешних зависимостей; Apple
// clang не умеет -fopenmp без brew libomp). Используется для параллельного
// обхода ячеек феррита. Каждая ячейка пишет только свои данные, суммы
// собираются ПОСЛЕДОВАТЕЛЬНО в фиксированном порядке -> результат не зависит
// от числа потоков (детерминирован, бит-в-бит как при 1 потоке).
// ----------------------------------------------------------------------------
class ParallelPool {
public:
    explicit ParallelPool(int n_threads) : n_(std::max(1, n_threads)) {
        for (int t = 1; t < n_; ++t) workers_.emplace_back([this] { worker_loop(); });
    }
    ~ParallelPool() {
        { std::lock_guard<std::mutex> lk(m_); stop_ = true; ++gen_; }
        cv_.notify_all();
        for (auto& w : workers_) w.join();
    }
    ParallelPool(const ParallelPool&) = delete;
    ParallelPool& operator=(const ParallelPool&) = delete;
    int size() const { return n_; }

    void parallel_for(int n_items, const std::function<void(int)>& fn) {
        if (n_ <= 1 || n_items <= 1) { for (int i = 0; i < n_items; ++i) fn(i); return; }
        {
            std::lock_guard<std::mutex> lk(m_);
            job_ = &fn; n_items_ = n_items; next_.store(0);
            active_ = (int)workers_.size(); ++gen_;
        }
        cv_.notify_all();
        run_items(fn);
        std::unique_lock<std::mutex> lk(m_);
        cv_done_.wait(lk, [this] { return active_ == 0; });
        job_ = nullptr;
    }
private:
    void run_items(const std::function<void(int)>& fn) {
        for (;;) { int i = next_.fetch_add(1); if (i >= n_items_) break; fn(i); }
    }
    void worker_loop() {
        long long seen = 0;
        for (;;) {
            const std::function<void(int)>* job = nullptr;
            {
                std::unique_lock<std::mutex> lk(m_);
                cv_.wait(lk, [&] { return gen_ != seen; });
                seen = gen_;
                if (stop_) return;
                job = job_;
            }
            run_items(*job);
            {
                std::lock_guard<std::mutex> lk(m_);
                if (--active_ == 0) cv_done_.notify_one();
            }
        }
    }
    int n_;
    std::vector<std::thread> workers_;
    std::mutex m_;
    std::condition_variable cv_, cv_done_;
    const std::function<void(int)>* job_ = nullptr;
    int n_items_ = 0;
    std::atomic<int> next_{0};
    int active_ = 0;
    long long gen_ = 0;
    bool stop_ = false;
};

struct JAEngine {
    long double Ms, a, alpha, k, c, sigma_m_leak, H0_bias;
    int n_sub, n_fp;
    JAEngine(long double Ms_, long double a_, long double alpha_, long double k_, long double c_,
             long double sml, long double H0b, int ns, int nf)
        : Ms(Ms_), a(a_), alpha(alpha_), k(k_), c(c_),
          sigma_m_leak(sml), H0_bias(H0b), n_sub(ns), n_fp(nf) {}

    inline void M_of_Mirr(long double M_irr, long double H, long double M_guess,
                           long double& M_out, long double& Man_out) const {
        long double M = M_guess;
        for (int i = 0; i < n_fp; ++i) {
            long double He = (H + H0_bias) + alpha * M;
            long double Man = Ms * langevin_L(He / a);
            M = (1.0 - c) * M_irr + c * Man;
        }
        long double He = (H + H0_bias) + alpha * M;
        long double Man = Ms * langevin_L(He / a);
        M_out = M; Man_out = Man;
    }

    inline void rhs(long double M_irr, long double H, long double delta, long double M_guess,
                     long double& drhs, long double& M_out) const {
        long double M, Man;
        M_of_Mirr(M_irr, H, M_guess, M, Man);
        long double denom = k * delta - alpha * (Man - M_irr);
        long double denom_safe = (std::fabs(denom) < 1e-9)
                             ? 1e-9 * (denom >= 0 ? 1.0 : -1.0) : denom;
        drhs = (Man - M_irr) / denom_safe;
        M_out = M;
    }

    void integrate_rk4(long double M_irr_old, long double H_old, long double H_new,
                        long double M_guess_old, long double delta,
                        long double& M_irr_out, long double& M_out) const {
        long double h = (H_new - H_old) / n_sub;
        long double M_irr = M_irr_old, M_guess = M_guess_old, H_cur = H_old;
        for (int s = 0; s < n_sub; ++s) {
            long double k1, M1; rhs(M_irr, H_cur, delta, M_guess, k1, M1);
            long double k2, M2; rhs(M_irr + 0.5*h*k1, H_cur + 0.5*h, delta, M1, k2, M2);
            long double k3, M3; rhs(M_irr + 0.5*h*k2, H_cur + 0.5*h, delta, M2, k3, M3);
            long double k4, M4; rhs(M_irr + h*k3, H_cur + h, delta, M3, k4, M4);
            M_irr = M_irr + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4);
            M_irr = std::max(-Ms, std::min(Ms, M_irr));
            H_cur += h; M_guess = M4;
        }
        long double M_final, Man_final;
        M_of_Mirr(M_irr, H_new, M_guess, M_final, Man_final);
        M_irr_out = M_irr;
        M_out = std::max(-Ms, std::min(Ms, M_final));
    }

    // Производная dM/dH считается ЧИСЛЕННО (конечная разность), а не по
    // аналитической формуле -- избегаем сингулярности 1/Heff^2 при Heff->0,
    // которая присутствует в наивной прямой реализации.
    void dM_dH(long double M_irr_old, long double H_old, long double H_new,
               long double M_guess_old, long double delta,
               long double& dMdH_out, long double& M0_out) const {
        long double eps = 1e-6 * std::max((long double)1.0, std::fabs(H_new));
        long double tmp1, M0; integrate_rk4(M_irr_old, H_old, H_new, M_guess_old, delta, tmp1, M0);
        long double tmp2, M1; integrate_rk4(M_irr_old, H_old, H_new + eps, M_guess_old, delta, tmp2, M1);
        dMdH_out = (M1 - M0) / eps; M0_out = M0;
    }

    void newton_branch(long double M_irr_old, long double H_old, long double M_old, long double rhs_extra,
                        long double dt, long double delta, int n_newton, long double max_step,
                        long double& H_new_out, long double& F_out) const {
        long double H_new = H_old - rhs_extra;
        long double leak_coef = 4*M_PI*sigma_m_leak*dt;
        long double dMdH, Mc;
        for (int it = 0; it < n_newton; ++it) {
            dM_dH(M_irr_old, H_old, H_new, M_old, delta, dMdH, Mc);
            long double F = H_new - H_old + 4*M_PI*(Mc - M_old) + leak_coef*H_new + rhs_extra;
            long double Fp = 1.0 + 4*M_PI*dMdH + leak_coef;
            long double Fp_safe = (std::fabs(Fp) < 0.5) ? 0.5*(Fp >= 0 ? 1.0 : -1.0) : Fp;
            long double step = std::max(-max_step, std::min(max_step, F/Fp_safe));
            H_new -= step;
        }
        dM_dH(M_irr_old, H_old, H_new, M_old, delta, dMdH, Mc);
        F_out = H_new - H_old + 4*M_PI*(Mc - M_old) + leak_coef*H_new + rhs_extra;
        H_new_out = H_new;
    }

    // Двусторонняя проверка знака delta=sign(dH/dt), фиксируется единожды по
    // знаку внешнего толчка -- избегает паразитного цикла период-2.
    void implicit_step(long double M_irr_old, long double H_old, long double M_old, long double rhs_extra,
                        long double dt, int n_newton, long double max_step_factor, long double tol,
                        long double& H_new_out, long double& M_irr_out, long double& M_out) const {
        long double max_step = max_step_factor * a;
        long double delta_pred = (-rhs_extra >= 0) ? 1.0 : -1.0;
        long double H_a, F_a;
        newton_branch(M_irr_old, H_old, M_old, rhs_extra, dt, delta_pred, n_newton, max_step, H_a, F_a);
        long double H_new = H_a, delta_final = delta_pred;
        if (std::fabs(F_a) > tol) {
            long double delta_alt = -delta_pred;
            long double H_b, F_b;
            newton_branch(M_irr_old, H_old, M_old, rhs_extra, dt, delta_alt, n_newton, max_step, H_b, F_b);
            if (std::fabs(F_b) < std::fabs(F_a)) { H_new = H_b; delta_final = delta_alt; }
        }
        long double M_irr_f, M_f;
        integrate_rk4(M_irr_old, H_old, H_new, M_old, delta_final, M_irr_f, M_f);
        H_new_out = H_new; M_irr_out = M_irr_f; M_out = M_f;
    }
};

struct Vec3 { long double x, y, z; };

struct LLGEngine {
    long double Ms, gamma, alpha, H0_bias;
    int bias_axis; // 0='none', 1='x', 2='y', 3='z'
    LLGEngine(long double Ms_, long double gamma_, long double alpha_, long double H0_bias_, int bias_axis_)
        : Ms(Ms_), gamma(gamma_), alpha(alpha_), H0_bias(H0_bias_), bias_axis(bias_axis_) {}

    inline Vec3 build_Heff(long double Hy_rf, long double Hz_rf) const {
        Vec3 H{0,0,0};
        if (bias_axis == 1) { H.x = H0_bias; H.y = Hy_rf; H.z = Hz_rf; }
        else if (bias_axis == 2) { H.y = H0_bias + Hy_rf; H.z = Hz_rf; }
        else if (bias_axis == 3) { H.y = Hy_rf; H.z = H0_bias + Hz_rf; }
        else { H.y = Hy_rf; H.z = Hz_rf; }
        return H;
    }
    inline Vec3 initial_M() const {
        Vec3 M{0,0,0};
        if (bias_axis == 1) M.x = Ms; else if (bias_axis == 2) M.y = Ms; else M.z = Ms;
        return M;
    }
    static inline Vec3 cross(const Vec3&a, const Vec3&b) {
        return Vec3{a.y*b.z-a.z*b.y, a.z*b.x-a.x*b.z, a.x*b.y-a.y*b.x};
    }
    static inline long double dot(const Vec3&a, const Vec3&b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
    static inline long double norm(const Vec3&a) { return std::sqrt(dot(a,a)); }

    Vec3 rodrigues_rotate(const Vec3& M, const Vec3& Heff, long double dt) const {
        long double Hn = norm(Heff);
        if (Hn <= 1e-300) return M;
        Vec3 axis{Heff.x/Hn, Heff.y/Hn, Heff.z/Hn};
        long double theta = gamma * Hn * dt;
        long double ct = std::cos(theta), st = std::sin(theta);
        long double d = dot(axis, M);
        Vec3 crs = cross(axis, M);
        return Vec3{ M.x*ct + crs.x*st + axis.x*d*(1-ct),
                     M.y*ct + crs.y*st + axis.y*d*(1-ct),
                     M.z*ct + crs.z*st + axis.z*d*(1-ct) };
    }
    Vec3 llg_step_core(const Vec3& M, const Vec3& Heff, long double dt) const {
        Vec3 M1 = rodrigues_rotate(M, Heff, dt);
        Vec3 cross1 = cross(M1, Heff);
        Vec3 damp = cross(M1, cross1);
        Vec3 M2{ M1.x - dt*(gamma*alpha/Ms)*damp.x,
                 M1.y - dt*(gamma*alpha/Ms)*damp.y,
                 M1.z - dt*(gamma*alpha/Ms)*damp.z };
        long double n2 = norm(M2);
        long double n2s = (n2 > 1e-300) ? n2 : 1.0;
        return Vec3{ M2.x*(Ms/n2s), M2.y*(Ms/n2s), M2.z*(Ms/n2s) };
    }
    void step(const Vec3& M, long double Hy_rf, long double Hz_rf, long double dt,
              Vec3& M_new_out, long double& dMy_dt_out, long double& dMz_dt_out) const {
        Vec3 Heff = build_Heff(Hy_rf, Hz_rf);
        Vec3 M_new = llg_step_core(M, Heff, dt);
        dMy_dt_out = (M_new.y - M.y) / dt;
        dMz_dt_out = (M_new.z - M.z) / dt;
        M_new_out = M_new;
    }
};

// Классическая скалярная модель Прейзаха: ансамбль элементарных реле
// (гистеронов) gamma_{alpha,beta}[H], пороги включения/выключения alpha_k>=beta_k
// заданы через физически интерпретируемые Hc=(alpha-beta)/2 (коэрцитивность)
// и Hb=(alpha+beta)/2 (поле взаимодействия/смещения гистерона).
// NOTE: std::normal_distribution зависит от реализации стандартной библиотеки
// (libstdc++ в Linux g++ и libc++ в macOS clang дают РАЗНЫЕ выборки при том же
// seed) -> результаты Preisach (1D) на разных платформах статистически
// эквивалентны, но не совпадают бит-в-бит. Preisach2D от этого свободен.
struct PreisachEngine {
    long double Ms, Hc_mean, Hc_sigma, Hb_sigma;
    int n_hyst;
    std::vector<long double> alpha_k, beta_k, weight_k;

    PreisachEngine(long double Ms_, long double Hc_mean_, long double Hc_sigma_,
                   long double Hb_sigma_, int n_hyst_)
        : Ms(Ms_), Hc_mean(Hc_mean_), Hc_sigma(Hc_sigma_), Hb_sigma(Hb_sigma_), n_hyst(n_hyst_)
    {
        // FIX: защита от n_hyst<=0 и sigma<=0 (UB в std::normal_distribution).
        if (n_hyst < 1) throw std::invalid_argument("Preisach: n_hyst_pr должно быть >= 1");
        if (!(Hc_sigma > 0) || !(Hb_sigma > 0))
            throw std::invalid_argument("Preisach: Hc_sigma_pr и Hb_sigma_pr должны быть > 0");
        alpha_k.resize(n_hyst); beta_k.resize(n_hyst); weight_k.resize(n_hyst);
        std::mt19937 rng(12345); // детерминированный посев -- воспроизводимость сканов
        std::normal_distribution<double> distHc(static_cast<double>(Hc_mean), static_cast<double>(Hc_sigma));
        std::normal_distribution<double> distHb(0.0, static_cast<double>(Hb_sigma));
        long double wsum = 0.0;
        for (int k = 0; k < n_hyst; ++k) {
            long double Hc = std::fabs((long double)distHc(rng));
            long double Hb = (long double)distHb(rng);
            alpha_k[k] = Hb + Hc;
            beta_k[k] = Hb - Hc;
            weight_k[k] = 1.0;
            wsum += 1.0;
        }
        for (auto& w : weight_k) w *= (Ms / wsum);
    }

    // Квазистатическое (rate-independent) переключение реле по направлению
    // движения поля H_old->H_new; state[k] -- текущий выход гистерона (+-1).
    void step(int8_t* state, long double H_old, long double H_new,
              long double& M_out) const {
        bool rising = H_new >= H_old;
        for (int k = 0; k < n_hyst; ++k) {
            if (rising) { if (H_new >= alpha_k[k]) state[k] = 1; }
            else { if (H_new <= beta_k[k]) state[k] = -1; }
        }
        long double M = 0.0;
        for (int k = 0; k < n_hyst; ++k) M += weight_k[k] * state[k];
        M_out = M;
    }

    // PERF: пробное вычисление M без копирования вектора состояния (в оригинале
    // на каждой итерации делалась копия std::vector<int8_t> trial = state --
    // лишняя аллокация n_hyst байт). Арифметика и порядок суммирования те же.
    long double trial_M(const int8_t* state, long double H_old, long double H_new) const {
        bool rising = H_new >= H_old;
        long double M = 0.0;
        for (int k = 0; k < n_hyst; ++k) {
            int8_t s = state[k];
            if (rising) { if (H_new >= alpha_k[k]) s = 1; }
            else { if (H_new <= beta_k[k]) s = -1; }
            M += weight_k[k] * s;
        }
        return M;
    }

    // Неявная состыковка с уравнением ротора H по аналогии с
    // JAEngine::newton_branch, но без аналитической dM/dH -- ансамбль реле
    // rate-independent, поэтому используем демпфированную квазиньютоновскую
    // итерацию с консервативным знаменателем (1 + leak_coef).
    void implicit_step(int8_t* state, long double H_old, long double M_old,
                        long double rhs_extra, long double leak_coef, int n_iter,
                        long double& H_new_out, long double& M_new_out) const {
        long double H_new = H_old - rhs_extra;
        long double M_new = M_old;
        const long double tol = 1e-10;
        for (int it = 0; it < n_iter; ++it) {
            M_new = trial_M(state, H_old, H_new);   // PERF: вместо копии trial
            long double F = H_new - H_old + 4*M_PI*(M_new - M_old) + leak_coef*H_new + rhs_extra;
            long double dM_dH = (M_new - M_old) / (H_new - H_old + 1e-15);
            long double denom = 1.0 + 4*M_PI*dM_dH + leak_coef;
            long double damping = 0.5;  // коэффициент демпфирования
            long double delta = F / denom;
            H_new -= damping * delta;

            if (std::abs(delta) < tol) break;  // выход при сходимости
        }
        step(state, H_old, H_new, M_new); // фиксируем состояния реле только после сходимости
        H_new_out = H_new; M_new_out = M_new;
    }
};

// ============================================================================
// PreisachEngine2D -- детерминированная сеточная реализация классического
// 2D оператора Прейзаха: пары (hc_i, hi_j), alpha=hi+hc, beta=hi-hc,
// вес w_ij = exp(-z_hc^2/2)*exp(-z_hi^2/2), нормировка sum w = Ms.
//
// FIX (память): в оригинале на КАЖДУЮ ячейку феррита хранился std::vector<int8_t>
//   из n_grid^2 элементов, в каждой итерации Ньютона делалась его копия, а сам
//   движок хранил 3 массива long double по n_grid^2. При n_grid = n_hyst_pr =
//   2001 это 4 004 001 пар на ячейку. Теперь:
//     * alpha/beta/weight НЕ хранятся (вычисляются из сеток hc, hi тем же
//       выражением -> те же значения);
//     * состояние -- либо бит-упаковка (PR2D_BITS), либо одно целое на строку
//       (PR2D_ROWS, см. описание enum выше);
//     * копия состояния в итерациях Ньютона удалена.
// FIX (индексация): n_pairs и индексы пар -- int64_t (в оригинале int; при
//   n_grid > 46340 n_grid*n_grid переполнял int -> UB/выход за границы).
// ============================================================================
struct PreisachEngine2D {
    long double Ms;
    int n_grid;        // разрешение сетки по каждой из осей (hc, hi)
    int64_t n_pairs;   // n_grid*n_grid (FIX: int64_t)
    std::vector<long double> hc_grid, hi_grid, w_hc, w_hi;
    long double scale;             // Ms / wsum; вес пары = (w_hc[i]*w_hi[j])*scale (как в оригинале)
    std::vector<long double> S_hi; // префиксные суммы w_hi (для PR2D_ROWS)
    long double S_tot;
    long double hi_lo, dhi;

    PreisachEngine2D(long double Ms_, long double Hc_mean, long double Hc_sigma,
                      long double Hb_sigma, int n_grid_)
        : Ms(Ms_), n_grid(n_grid_)
    {
        // FIX: проверки входа (n_grid<1 -> пустая сетка и деление на wsum=0;
        // sigma<=0 -> деление на ноль в gauss -> NaN-веса).
        if (n_grid < 1) throw std::invalid_argument("Preisach2D: n_grid должно быть >= 1");
        if (n_grid > 1000000) throw std::invalid_argument("Preisach2D: n_grid слишком велико (>1e6)");
        if (!(Hc_sigma > 0) || !(Hb_sigma > 0))
            throw std::invalid_argument("Preisach2D: Hc_sigma_pr и Hb_sigma_pr должны быть > 0");
        n_pairs = (int64_t)n_grid * (int64_t)n_grid;

        long double hc_lo = std::max((long double)0.0, Hc_mean - 4*Hc_sigma);
        long double hc_hi = Hc_mean + 4*Hc_sigma;
        hi_lo = -4*Hb_sigma; long double hi_hi = 4*Hb_sigma;
        hc_grid.resize(n_grid); hi_grid.resize(n_grid);
        for (int i = 0; i < n_grid; ++i) {
            long double f = (n_grid > 1) ? (long double)i/(n_grid-1) : 0.0;
            hc_grid[i] = hc_lo + (hc_hi-hc_lo)*f;
            hi_grid[i] = hi_lo + (hi_hi-hi_lo)*f;
        }
        dhi = (n_grid > 1) ? (hi_hi - hi_lo)/(n_grid - 1) : 0.0;
        auto gauss = [](long double x, long double mu, long double sigma) {
            long double z = (x-mu)/sigma;
            return std::exp((long double)(-0.5)*z*z);
        };
        w_hc.resize(n_grid); w_hi.resize(n_grid);
        for (int i = 0; i < n_grid; ++i) w_hc[i] = gauss(hc_grid[i], Hc_mean, Hc_sigma);
        for (int j = 0; j < n_grid; ++j) w_hi[j] = gauss(hi_grid[j], 0.0, Hb_sigma);
        // wsum накапливается в том же порядке (i внешний, j внутренний) и тем же
        // разбиением на операторы, что в оригинале (w = w_hc*w_hi; wsum += w).
        long double wsum = 0.0;
        for (int i = 0; i < n_grid; ++i)
            for (int j = 0; j < n_grid; ++j) {
                long double w = w_hc[i] * w_hi[j];
                wsum += w;
            }
        scale = Ms / wsum;
        S_hi.assign(n_grid + 1, 0.0);
        for (int j = 0; j < n_grid; ++j) S_hi[j+1] = S_hi[j] + w_hi[j];
        S_tot = S_hi[n_grid];
    }

    // ---------------- PR2D_BITS: бит-упакованное состояние ----------------
    static size_t bit_words(int64_t n_pairs_) { return (size_t)((n_pairs_ + 63) / 64); }
    static inline bool get_bit(const uint64_t* st, int64_t k) { return (st[k >> 6] >> (k & 63)) & 1ULL; }
    static inline void set_bit(uint64_t* st, int64_t k, bool v) {
        if (v) st[k >> 6] |= (1ULL << (k & 63)); else st[k >> 6] &= ~(1ULL << (k & 63));
    }
    // commit=false -- пробное вычисление (state не меняется), commit=true -- фиксация.
    // Порядок суммирования M идентичен оригинальному step() (k = i*n_grid + j).
    long double eval_bits(uint64_t* st, long double H_old, long double H_new, bool commit) const {
        bool rising = H_new >= H_old;
        long double M = 0.0;
        int64_t k = 0;
        for (int i = 0; i < n_grid; ++i) {
            const long double hc = hc_grid[i];
            const long double wh = w_hc[i];
            for (int j = 0; j < n_grid; ++j, ++k) {
                bool s = get_bit(st, k);
                const long double hi = hi_grid[j];
                const long double alpha = hi + hc, beta = hi - hc;
                if (rising) { if (H_new >= alpha) s = true; }
                else        { if (H_new <= beta)  s = false; }
                if (commit) set_bit(st, k, s);
                long double w = wh * w_hi[j];
                w = w * scale;
                M += w * (s ? 1 : -1);
            }
        }
        return M;
    }

    // ---------------- PR2D_ROWS: одно целое K_i на строку ----------------
    // Число j с H >= hi_j + hc_i (реле, включённые ростом поля до H).
    inline int count_up(int i, long double H) const {
        const long double hc = hc_grid[i];
        int g = 0;
        if (dhi > 0) {
            long double x = (H - hc - hi_lo) / dhi;
            if (!(x > -1)) g = 0;                  // также ловит NaN
            else if (x >= n_grid) g = n_grid;
            else g = (int)std::floor(x) + 1;
            g = std::max(0, std::min(n_grid, g));
        }
        // точная коррекция теми же сравнениями, что в оригинале
        while (g < n_grid && H >= hi_grid[g] + hc) ++g;
        while (g > 0 && !(H >= hi_grid[g-1] + hc)) --g;
        return g;
    }
    // Наименьшее j с H <= hi_j - hc_i (реле, выключаемые спадом поля до H).
    inline int first_dn(int i, long double H) const {
        const long double hc = hc_grid[i];
        int g = n_grid;
        if (dhi > 0) {
            long double x = (H + hc - hi_lo) / dhi;
            if (!(x > 0)) g = 0;                   // также ловит NaN
            else if (x >= n_grid) g = n_grid;
            else g = (int)std::ceil(x);
            g = std::max(0, std::min(n_grid, g));
        }
        while (g > 0 && H <= hi_grid[g-1] - hc) --g;
        while (g < n_grid && !(H <= hi_grid[g] - hc)) ++g;
        return g;
    }
    long double eval_rows(int32_t* K, long double H_old, long double H_new, bool commit) const {
        bool rising = H_new >= H_old;
        long double acc = 0.0;
        for (int i = 0; i < n_grid; ++i) {
            int k = K[i];
            int kn = rising ? std::max(k, count_up(i, H_new)) : std::min(k, first_dn(i, H_new));
            if (commit) K[i] = kn;
            acc += w_hc[i] * (2*S_hi[kn] - S_tot);
        }
        return acc * scale;
    }

    long double eval(uint64_t* st, long double H_old, long double H_new, bool commit) const {
        return eval_bits(st, H_old, H_new, commit);
    }
    long double eval(int32_t* K, long double H_old, long double H_new, bool commit) const {
        return eval_rows(K, H_old, H_new, commit);
    }

    // Интерфейс идентичен PreisachEngine::implicit_step -- прямая замена.
    // StateT = uint64_t (PR2D_BITS) или int32_t (PR2D_ROWS).
    template <class StateT>
    void implicit_step(StateT* state, long double H_old, long double M_old,
                        long double rhs_extra, long double leak_coef, int n_iter,
                        long double& H_new_out, long double& M_new_out) const {
        long double H_new = H_old - rhs_extra;
        long double M_new = M_old;
        const long double tol = 1e-10;
        for (int it = 0; it < n_iter; ++it) {
            M_new = eval(state, H_old, H_new, false);   // FIX/PERF: без копии состояния
            long double F = H_new - H_old + 4*M_PI*(M_new - M_old) + leak_coef*H_new + rhs_extra;
            long double dM_dH = (M_new - M_old) / (H_new - H_old + 1e-15);
            long double denom = 1.0 + 4*M_PI*dM_dH + leak_coef;
            long double damping = 0.5;  // коэффициент демпфирования
            long double delta = F / denom;
            H_new -= damping * delta;

            if (std::abs(delta) < tol) break;  // выход при сходимости
        }
        M_new = eval(state, H_old, H_new, true);
        H_new_out = H_new; M_new_out = M_new;
    }
};

// ---------------- помощники бинарной сериализации (чекпойнты) ----------------
template <class T>
static bool ck_write_vec(FILE* f, const std::vector<T>& v) {
    uint64_t n = v.size();
    if (fwrite(&n, sizeof(n), 1, f) != 1) return false;
    if (n && fwrite(v.data(), sizeof(T), n, f) != n) return false;
    return true;
}
template <class T>
static bool ck_read_vec(FILE* f, std::vector<T>& v) {
    uint64_t n = 0;
    if (fread(&n, sizeof(n), 1, f) != 1) return false;
    if (n != v.size()) return false;   // размер обязан совпасть с текущей конфигурацией
    if (n && fread(v.data(), sizeof(T), n, f) != n) return false;
    return true;
}
template <class T>
static bool ck_write_pod(FILE* f, const T& x) { return fwrite(&x, sizeof(T), 1, f) == 1; }
template <class T>
static bool ck_read_pod(FILE* f, T& x) { return fread(&x, sizeof(T), 1, f) == 1; }

struct MenDriveSim {
    int N;
    long double a, h_l, h_r, dx, dt;
    int ferrite_model;   // 0='JA', 1='LLG', 2='Hybrid', 3='Preisach', 4='Preisach2D'
    int excitation_mode; // 0 ='magnetic_right', - MenDrive_qnm-mode3_jm_z_TE.ipynb
                         // 1 ='electric_left'   - MenDrive_qnm-mode2_je_z_TM.ipynb
    int n_newton;

    std::vector<long double> xA, xB;
    std::vector<uint8_t> leftA, rightA, leftB, rightB;
    int iL_A, iR_A, n_fer, n_left;
    std::vector<int> fer_cells;   // PERF: глобальные индексы ячеек феррита (тот же порядок, что fi)

    std::vector<long double> Ca_e, Cb_e, src_mag_B, src_el_A, Ey, Ez, Hy, Hz;
    long double t;
    std::vector<long double> sigma_e_profile;  // диагностика: sigma_e(x) на A-сетке

    std::vector<long double> M_irr, M_field;          // JA
    std::vector<Vec3> M_llg;                          // LLG
    std::vector<long double> hyb_M_irr, hyb_M_field;  // Hybrid (внутренний JA)
    std::vector<Vec3> hyb_M;                          // Hybrid (LLG-часть)
    std::vector<long double> last_HzJA;               // Hybrid (неискажённое поле)

    std::unique_ptr<JAEngine> ja;
    std::unique_ptr<LLGEngine> llg;
    std::unique_ptr<PreisachEngine> pr;                 // Preisach
    std::unique_ptr<PreisachEngine2D> pr2d;             // Preisach 2D
    std::vector<std::vector<int8_t>> pr_state;          // Preisach (состояния гистеронов; теперь только 1D)
    std::vector<std::vector<uint64_t>> pr2d_bits;       // FIX: Preisach2D, PR2D_BITS (1 бит на пару)
    std::vector<std::vector<int32_t>> pr2d_rows;        // FIX: Preisach2D, PR2D_ROWS (K_i на строку)
    int pr2d_algo;
    std::vector<long double> pr_M;                      // Preisach (намагниченность)
    long double sigma_m_leak_pr;                        // Preisach (утечка энергии в стенку)
    long double p_mag_dynamics;  // диагностика: H*dM/dt за последний шаг (энергия к намагниченности)
    long double p_num_residual;  // DIAG: мощность «численного источника» из невязки неявного шага

    // PERF: рабочие буферы шага (в оригинале выделялись заново на каждом шаге)
    std::vector<long double> j_e_y_A, j_m_z_B, Ey_new, Ez_new, rotE_z, Hz_new, Hy_new;
    std::vector<long double> p_cell, pnum_cell;
    int n_threads;
    std::unique_ptr<ParallelPool> pool;

    // bias_orientation: 0='none', 1='x', 2='y', 3='z', 4='parallel'
    MenDriveSim(int N_, long double a_, long double h_l_, long double h_r_, long double dt_frac,
                int excitation_mode_, int ferrite_model_, int bias_orientation, long double H0_bias,
                long double sigma_m_leak, long double sigma_e_left, long double sigma_e_ferrite,
                long double Ms, long double a_JA, long double alpha_JA, long double k_JA, long double c_JA,
                int n_sub, int n_fp, int n_newton_,
                long double Ms_llg, long double gamma_llg, long double alpha_llg, long double skin_depth,
                int n_hyst_pr, long double Hc_mean_pr, long double Hc_sigma_pr, long double Hb_sigma_pr,
                long double sigma_m_leak_pr_, int pr2d_algo_ = PR2D_BITS, int n_threads_ = 1)
        : N(N_), a(a_), h_l(h_l_), h_r(h_r_),
          ferrite_model(ferrite_model_), excitation_mode(excitation_mode_), n_newton(n_newton_),
          pr2d_algo(pr2d_algo_), sigma_m_leak_pr(sigma_m_leak_pr_), n_threads(std::max(1, n_threads_))
    {
        // FIX: проверка входа (N<2 -> выход за границы в step_once).
        if (N < 2) throw std::invalid_argument("N должно быть >= 2");
        if (ferrite_model < JA || ferrite_model > Preisach2D) throw std::invalid_argument("неизвестный ferrite_model");
        if (pr2d_algo != PR2D_BITS && pr2d_algo != PR2D_ROWS) throw std::invalid_argument("неизвестный pr2d_algo");

        long double x_min = -(a + h_l), x_max = (a + h_r);
        dx = (x_max - x_min) / N; dt = dt_frac * dx;

        xA.resize(N+1); xB.resize(N);
        leftA.resize(N+1); rightA.resize(N+1); leftB.resize(N); rightB.resize(N);
        for (int i = 0; i <= N; ++i) { xA[i] = x_min + dx*i; leftA[i] = xA[i] < -a; rightA[i] = xA[i] > a; }
        for (int i = 0; i < N; ++i) { xB[i] = x_min + dx*(i + 0.5); leftB[i] = xB[i] < -a; rightB[i] = xB[i] > a; }

        iL_A = 0; iR_A = 0;
        long double bestL = 1e300, bestR = 1e300;
        for (int i = 0; i <= N; ++i) {
            long double dl = std::fabs(xA[i] + a); if (dl < bestL) { bestL = dl; iL_A = i; }
            long double dr = std::fabs(xA[i] - a); if (dr < bestR) { bestR = dr; iR_A = i; }
        }
        n_fer = 0; for (int i = 0; i < N; ++i) if (rightB[i]) n_fer++;
        n_left = 0; for (int i = 0; i < N; ++i) if (leftB[i]) n_left++;
        fer_cells.clear(); for (int i = 0; i < N; ++i) if (rightB[i]) fer_cells.push_back(i);

        if (ferrite_model == JA) {
            long double H0_eff = (bias_orientation == 4 || bias_orientation == 3) ? H0_bias : 0.0;
            ja.reset(new JAEngine(Ms, a_JA, alpha_JA, k_JA, c_JA, sigma_m_leak, H0_eff, n_sub, n_fp));
        } else if (ferrite_model == LLG) {
            int axis = (bias_orientation >= 1 && bias_orientation <= 3) ? bias_orientation : 0;
            llg.reset(new LLGEngine(Ms_llg, gamma_llg, alpha_llg, H0_bias, axis));
        } else if (ferrite_model == Hybrid) {
            ja.reset(new JAEngine(Ms, a_JA, alpha_JA, k_JA, c_JA, sigma_m_leak, 0.0, n_sub, n_fp));
            int axis = (bias_orientation >= 1 && bias_orientation <= 3) ? bias_orientation : 0;
            llg.reset(new LLGEngine(Ms_llg, gamma_llg, alpha_llg, H0_bias, axis));
        } else if (ferrite_model == Preisach) {
            pr.reset(new PreisachEngine(Ms, Hc_mean_pr, Hc_sigma_pr, Hb_sigma_pr, n_hyst_pr));
        } else if (ferrite_model == Preisach2D) {
            // NOTE: как и в оригинале, для Preisach2D параметр n_hyst_pr трактуется
            // как n_grid (число узлов сетки по КАЖДОЙ оси) -> n_grid^2 пар!
            // n_hyst_pr=2001 (типичное значение для 1D Preisach) даёт 4 004 001 пар.
            pr2d.reset(new PreisachEngine2D(Ms, Hc_mean_pr, Hc_sigma_pr, Hb_sigma_pr, n_hyst_pr));
        }

        Ca_e.resize(N+1); Cb_e.resize(N+1); sigma_e_profile.resize(N+1);
        for (int i = 0; i <= N; ++i) {
            long double sig = leftA[i] ? sigma_e_left : (rightA[i] ? sigma_e_ferrite : 0.0);
            Ca_e[i] = (1 - 2*M_PI*sig*dt) / (1 + 2*M_PI*sig*dt);
            Cb_e[i] = (dt/dx) / (1 + 2*M_PI*sig*dt);
            sigma_e_profile[i] = sig;
        }
        src_mag_B.resize(N); src_el_A.resize(N+1);
        for (int i = 0; i < N; ++i) src_mag_B[i] = rightB[i] ? std::exp(-(xB[i]-a)/skin_depth) : 0.0;
        for (int i = 0; i <= N; ++i) src_el_A[i] = leftA[i] ? std::exp(-(-a-xA[i])/skin_depth) : 0.0;

        j_e_y_A.assign(N+1, 0.0); j_m_z_B.assign(N, 0.0);
        Ey_new.assign(N+1, 0.0); Ez_new.assign(N+1, 0.0); rotE_z.assign(N, 0.0);
        Hz_new.assign(N, 0.0); Hy_new.assign(N, 0.0);
        p_cell.assign(n_fer, 0.0); pnum_cell.assign(n_fer, 0.0);
        if (n_threads > 1) pool.reset(new ParallelPool(n_threads));
        reset_state();
    }

    void reset_state() {
        // TODO: switch TE or TM depending on excitation_mode
        Ey.assign(N+1, 0.0); Ez.assign(N+1, 0.0); Hy.assign(N, 0.0); Hz.assign(N, 0.0); t = 0.0;
        p_mag_dynamics = 0.0; p_num_residual = 0.0;
        if (ferrite_model == JA) { M_irr.assign(n_fer, 0.0); M_field.assign(n_fer, 0.0); }
        else if (ferrite_model == LLG) { Vec3 M0 = llg->initial_M(); M_llg.assign(n_fer, M0); }
        else if (ferrite_model == Hybrid) {
            hyb_M_irr.assign(n_fer, 0.0); hyb_M_field.assign(n_fer, 0.0);
            Vec3 M0 = llg->initial_M(); hyb_M.assign(n_fer, M0);
            last_HzJA.assign(n_fer, 0.0);
        }
        else if (ferrite_model == Preisach) {
            pr_state.assign(n_fer, std::vector<int8_t>(pr->n_hyst, -1));
            pr_M.assign(n_fer, -pr->Ms);
        }
        else if (ferrite_model == Preisach2D) {
            // FIX: было pr_state.assign(n_fer, std::vector<int8_t>(n_pairs, -1)) --
            // n_fer*n_grid^2 байт. Бит 0 == состояние -1 (начальное), бит 1 == +1.
            if (pr2d_algo == PR2D_BITS)
                pr2d_bits.assign(n_fer, std::vector<uint64_t>(PreisachEngine2D::bit_words(pr2d->n_pairs), 0ULL));
            else
                pr2d_rows.assign(n_fer, std::vector<int32_t>(pr2d->n_grid, 0));
            pr_M.assign(n_fer, -pr2d->Ms);
        }
    }

    // Память под состояние ферритовой модели (байт) -- для печати оценки.
    double state_bytes() const {
        if (ferrite_model == Preisach) return (double)n_fer * pr->n_hyst;
        if (ferrite_model == Preisach2D) {
            if (pr2d_algo == PR2D_BITS) return (double)n_fer * 8.0 * PreisachEngine2D::bit_words(pr2d->n_pairs);
            return (double)n_fer * 4.0 * pr2d->n_grid;
        }
        return (double)n_fer * 9.0 * sizeof(long double);
    }

    long double Txx_at(int iA) const {
        long double Hy_A = (iA > 0 && iA < N) ? 0.5*(Hy[iA-1]+Hy[iA]) : 0.0;
        long double Hz_A = (iA > 0 && iA < N) ? 0.5*(Hz[iA-1]+Hz[iA]) : 0.0;
        return (Ey[iA]*Ey[iA] + Ez[iA]*Ez[iA] + Hy_A*Hy_A + Hz_A*Hz_A) / (8*M_PI);
    }

    // PERF: параллельный (или последовательный) обход ячеек феррита.
    template <class F>
    void for_each_fer(F&& f) {
#ifdef USE_OPENMP
        if (n_threads > 1) {
            #pragma omp parallel for schedule(static) num_threads(n_threads)
            for (int fi = 0; fi < n_fer; ++fi) f(fi);
            return;
        }
#endif
        if (pool && n_threads > 1) { std::function<void(int)> fn = f; pool->parallel_for(n_fer, fn); }
        else { for (int fi = 0; fi < n_fer; ++fi) f(fi); }
    }

    long double step_once(long double omega0, long double amp, long double ramp_periods, long double T) {
        long double ramp = std::min((long double)1.0, t / (ramp_periods * T));
        long double s = amp * ramp * std::sin(omega0 * t);
        // PERF: буферы переиспользуются; обнуление эквивалентно созданию vector(N, 0.0)
        std::fill(j_e_y_A.begin(), j_e_y_A.end(), (long double)0.0);
        std::fill(j_m_z_B.begin(), j_m_z_B.end(), (long double)0.0);
        if (excitation_mode == 0) { for (int i = 0; i < N; ++i) j_m_z_B[i] = s * src_mag_B[i]; }
        else { for (int i = 0; i <= N; ++i) j_e_y_A[i] = s * src_el_A[i]; }

        std::fill(Ey_new.begin(), Ey_new.end(), (long double)0.0);
        std::fill(Ez_new.begin(), Ez_new.end(), (long double)0.0);
        for (int i = 1; i < N; ++i) {
            Ey_new[i] = Ca_e[i]*Ey[i] - Cb_e[i]*(Hz[i]-Hz[i-1]) - dt*4*M_PI*j_e_y_A[i];
            Ez_new[i] = Ca_e[i]*Ez[i] + Cb_e[i]*(Hy[i]-Hy[i-1]);
        }
        Ey_new[0]=Ey_new[N]=0.0; Ez_new[0]=Ez_new[N]=0.0;
        long double P_e = 0.0;
        for (int i = 0; i <= N; ++i) P_e -= j_e_y_A[i] * 0.5*(Ey[i]+Ey_new[i]) * dx;
        Ey.swap(Ey_new); Ez.swap(Ez_new);   // PERF: вместо копирования Ey = Ey_new

        for (int i = 0; i < N; ++i) rotE_z[i] = (Ey[i+1]-Ey[i]) / dx;

        std::fill(Hz_new.begin(), Hz_new.end(), (long double)0.0);
        std::fill(Hy_new.begin(), Hy_new.end(), (long double)0.0);
        for (int i = 0; i < N; ++i) Hy_new[i] = Hy[i] + (dt/dx)*(Ez[i+1]-Ez[i]);
        for (int i = 0; i < N; ++i) { if (!rightB[i]) Hz_new[i] = Hz[i] - dt*rotE_z[i] - dt*4*M_PI*j_m_z_B[i]; }

        if (ferrite_model == JA) {
            const long double leak_coef = 4*M_PI*ja->sigma_m_leak*dt;  // DIAG (тот же leak, что в newton_branch)
            for_each_fer([&](int fi) {
                int i = fer_cells[fi];
                long double rhs_extra = dt*rotE_z[i] + dt*4*M_PI*j_m_z_B[i];
                long double Hz_f, M_irr_new, M_new;
                ja->implicit_step(M_irr[fi], Hz[i], M_field[fi], rhs_extra, dt, n_newton, 0.25, 1e-5, Hz_f, M_irr_new, M_new);
                p_cell[fi] = 0.5*(Hz[i]+Hz_f) * (M_new-M_field[fi]) / dt * dx;
                // DIAG: невязка дискретного закона Фарадея после неявного шага
                long double F_res = Hz_f - Hz[i] + 4*M_PI*(M_new - M_field[fi]) + leak_coef*Hz_f + rhs_extra;
                pnum_cell[fi] = 0.5*(Hz[i]+Hz_f) * F_res / (4*M_PI*dt) * dx;
                Hz_new[i] = Hz_f; M_irr[fi] = M_irr_new; M_field[fi] = M_new;
            });
            long double p_mag_acc = 0.0, p_num_acc = 0.0;
            for (int fi = 0; fi < n_fer; ++fi) { p_mag_acc += p_cell[fi]; p_num_acc += pnum_cell[fi]; }
            p_mag_dynamics = p_mag_acc; p_num_residual = p_num_acc;
        } else if (ferrite_model == LLG) {
            int fi = 0;
            long double p_mag_acc = 0.0;
            for (int i = 0; i < N; ++i) {
                if (!rightB[i]) continue;
                Vec3 M_new; long double dMy_dt, dMz_dt;
                llg->step(M_llg[fi], Hy[i], Hz[i], dt, M_new, dMy_dt, dMz_dt);
                p_mag_acc += (Hy[i]*dMy_dt + Hz[i]*dMz_dt) * dx;
                M_llg[fi] = M_new;
                Hy_new[i] = Hy[i] + (dt/dx)*(Ez[i+1]-Ez[i]) - dt*4*M_PI*dMy_dt;
                Hz_new[i] = Hz[i] - dt*rotE_z[i] - dt*4*M_PI*j_m_z_B[i] - dt*4*M_PI*dMz_dt;
                fi++;
            }
            p_mag_dynamics = p_mag_acc; p_num_residual = 0.0;
        } else if (ferrite_model == Hybrid) {
            int fi = 0;
            long double p_mag_acc = 0.0;
            for (int i = 0; i < N; ++i) {
                if (!rightB[i]) continue;
                long double rhs_extra_z = dt*rotE_z[i] + dt*4*M_PI*j_m_z_B[i];
                long double Hz_after_JA, M_irr_new, M_field_new;
                ja->implicit_step(hyb_M_irr[fi], Hz[i], hyb_M_field[fi], rhs_extra_z, dt, n_newton, 0.25, 1e-5, Hz_after_JA, M_irr_new, M_field_new);
                p_mag_acc += 0.5*(Hz[i]+Hz_after_JA) * (M_field_new-hyb_M_field[fi]) / dt * dx;
                hyb_M_irr[fi] = M_irr_new; hyb_M_field[fi] = M_field_new;
                last_HzJA[fi] = Hz_after_JA;
                Vec3 M_new; long double dMy_dt, dMz_dt;
                llg->step(hyb_M[fi], Hy[i], Hz_after_JA, dt, M_new, dMy_dt, dMz_dt);
                p_mag_acc += (Hy[i]*dMy_dt + Hz_after_JA*dMz_dt) * dx;
                hyb_M[fi] = M_new;
                Hy_new[i] = Hy[i] + (dt/dx)*(Ez[i+1]-Ez[i]) - dt*4*M_PI*dMy_dt;
                Hz_new[i] = Hz_after_JA - dt*4*M_PI*dMz_dt;
                fi++;
            }
            p_mag_dynamics = p_mag_acc; p_num_residual = 0.0;
        } else if (ferrite_model == Preisach || ferrite_model == Preisach2D) {
            const long double leak_coef = 4*M_PI*sigma_m_leak_pr*dt;
            for_each_fer([&](int fi) {
                int i = fer_cells[fi];
                long double rhs_extra = dt*rotE_z[i] + dt*4*M_PI*j_m_z_B[i];
                long double Hz_f, M_new;
                if (ferrite_model == Preisach)
                    pr->implicit_step(pr_state[fi].data(), Hz[i], pr_M[fi], rhs_extra, leak_coef, n_newton, Hz_f, M_new);
                else if (pr2d_algo == PR2D_BITS)
                    pr2d->implicit_step(pr2d_bits[fi].data(), Hz[i], pr_M[fi], rhs_extra, leak_coef, n_newton, Hz_f, M_new);
                else
                    pr2d->implicit_step(pr2d_rows[fi].data(), Hz[i], pr_M[fi], rhs_extra, leak_coef, n_newton, Hz_f, M_new);
                p_cell[fi] = 0.5*(Hz[i]+Hz_f) * (M_new-pr_M[fi]) / dt * dx;
                // DIAG: невязка F после фиксации состояния реле. Реальное приращение
                // B = H+4piM в ячейке равно (-rhs_extra - leak*H_new + F_res), т.е.
                // F_res действует как фиктивный магнитный ток -> мощность H*F_res/(4pi dt).
                long double F_res = Hz_f - Hz[i] + 4*M_PI*(M_new - pr_M[fi]) + leak_coef*Hz_f + rhs_extra;
                pnum_cell[fi] = 0.5*(Hz[i]+Hz_f) * F_res / (4*M_PI*dt) * dx;
                Hz_new[i] = Hz_f; pr_M[fi] = M_new;
            });
            long double p_mag_acc = 0.0, p_num_acc = 0.0;
            for (int fi = 0; fi < n_fer; ++fi) { p_mag_acc += p_cell[fi]; p_num_acc += pnum_cell[fi]; }
            p_mag_dynamics = p_mag_acc; p_num_residual = p_num_acc;
        }

        long double P_m = 0.0;
        for (int i = 0; i < N; ++i) P_m -= j_m_z_B[i] * 0.5*(Hz[i]+Hz_new[i]) * dx;
        Hy.swap(Hy_new); Hz.swap(Hz_new); t += dt;   // PERF: вместо копирования
        // NOTE: возвращается МОЩНОСТЬ источника -j.E (или -j_m.H), проинтегрированная
        // по x (энергия / время / площадь). Множителя dt здесь нет -- см. mendrive_run.
        return P_e + P_m;
    }

    bool fields_finite() const { for (long double v : Hz) if (!std::isfinite(v)) return false; return true; }

    // Диагностика энергобаланса: прямой временной интеграл закона Джоуля
    // sigma_e*E^2 по обеим проводящим стенкам (профиль sigma_e_profile на
    // A-сетке) + sigma_m_leak*H^2 по ферритовой стенке (линейный вязкий
    // канал JA/Hybrid/Preisach) + H*dM/dt (энергия, передаваемая
    // намагниченности -- площадь петли гистерезиса для JA/Preisach/Hybrid,
    // работа демпфирования Гильберта для LLG).
    long double dissipated_power_density() const {
        long double p_diss = 0.0;
        for (int i = 0; i <= N; ++i) {
            p_diss += sigma_e_profile[i] * (Ey[i]*Ey[i] + Ez[i]*Ez[i]) * dx;
        }
        long double sigma_m_wall = 0.0;
        if (ferrite_model == JA || ferrite_model == Hybrid) sigma_m_wall = ja->sigma_m_leak;
        else if (ferrite_model == Preisach || ferrite_model == Preisach2D) sigma_m_wall = sigma_m_leak_pr;
        if (sigma_m_wall != 0.0) {
            for (int i = 0; i < N; ++i) {
                if (!rightB[i]) continue;
                p_diss += sigma_m_wall * (Hy[i]*Hy[i] + Hz[i]*Hz[i]) * dx;
            }
        }
        p_diss += p_mag_dynamics;
        return p_diss;
    }

    // DIAG: те же три слагаемых по отдельности (для разбора баланса).
    void dissipated_power_parts(long double& p_e, long double& p_m_wall, long double& p_mag) const {
        p_e = 0.0;
        for (int i = 0; i <= N; ++i) p_e += sigma_e_profile[i] * (Ey[i]*Ey[i] + Ez[i]*Ez[i]) * dx;
        long double sigma_m_wall = 0.0;
        if (ferrite_model == JA || ferrite_model == Hybrid) sigma_m_wall = ja->sigma_m_leak;
        else if (ferrite_model == Preisach || ferrite_model == Preisach2D) sigma_m_wall = sigma_m_leak_pr;
        p_m_wall = 0.0;
        if (sigma_m_wall != 0.0)
            for (int i = 0; i < N; ++i) if (rightB[i]) p_m_wall += sigma_m_wall * (Hy[i]*Hy[i] + Hz[i]*Hz[i]) * dx;
        p_mag = p_mag_dynamics;
    }

    // DIAG: полная энергия ЭМ поля (E^2+H^2)/(8 pi) на единицу площади.
    long double field_energy() const {
        long double W = 0.0;
        for (int i = 0; i <= N; ++i) W += (Ey[i]*Ey[i] + Ez[i]*Ez[i]) * dx;
        for (int i = 0; i < N; ++i) W += (Hy[i]*Hy[i] + Hz[i]*Hz[i]) * dx;
        return W / (8*M_PI);
    }

    // Глобальный индекс B-ячейки пробы (fi-я ячейка феррита) или -1.
    int probe_global_index(int probe_idx) const {
        return (probe_idx >= 0 && probe_idx < n_fer) ? fer_cells[probe_idx] : -1;
    }

    // ---------------- чекпойнт: полное динамическое состояние ----------------
    bool save_state(FILE* f) const {
        const uint32_t magic = 0x4D445354; // 'MDST'
        uint32_t ld = sizeof(long double);
        bool ok = ck_write_pod(f, magic) && ck_write_pod(f, ld) && ck_write_pod(f, ferrite_model)
               && ck_write_pod(f, pr2d_algo) && ck_write_pod(f, N) && ck_write_pod(f, n_fer)
               && ck_write_pod(f, t) && ck_write_pod(f, p_mag_dynamics) && ck_write_pod(f, p_num_residual)
               && ck_write_vec(f, Ey) && ck_write_vec(f, Ez) && ck_write_vec(f, Hy) && ck_write_vec(f, Hz);
        if (!ok) return false;
        if (ferrite_model == JA) ok = ck_write_vec(f, M_irr) && ck_write_vec(f, M_field);
        else if (ferrite_model == LLG) ok = ck_write_vec(f, M_llg);
        else if (ferrite_model == Hybrid)
            ok = ck_write_vec(f, hyb_M_irr) && ck_write_vec(f, hyb_M_field) && ck_write_vec(f, hyb_M) && ck_write_vec(f, last_HzJA);
        else if (ferrite_model == Preisach) {
            ok = ck_write_vec(f, pr_M);
            for (int fi = 0; ok && fi < n_fer; ++fi) ok = ck_write_vec(f, pr_state[fi]);
        } else if (ferrite_model == Preisach2D) {
            ok = ck_write_vec(f, pr_M);
            for (int fi = 0; ok && fi < n_fer; ++fi)
                ok = (pr2d_algo == PR2D_BITS) ? ck_write_vec(f, pr2d_bits[fi]) : ck_write_vec(f, pr2d_rows[fi]);
        }
        return ok;
    }
    bool load_state(FILE* f) {
        uint32_t magic = 0, ld = 0; int fm = -1, algo = -1, n = -1, nf = -1;
        if (!(ck_read_pod(f, magic) && ck_read_pod(f, ld) && ck_read_pod(f, fm) && ck_read_pod(f, algo)
              && ck_read_pod(f, n) && ck_read_pod(f, nf))) return false;
        if (magic != 0x4D445354 || ld != sizeof(long double) || fm != ferrite_model || algo != pr2d_algo
            || n != N || nf != n_fer) return false;
        bool ok = ck_read_pod(f, t) && ck_read_pod(f, p_mag_dynamics) && ck_read_pod(f, p_num_residual)
               && ck_read_vec(f, Ey) && ck_read_vec(f, Ez) && ck_read_vec(f, Hy) && ck_read_vec(f, Hz);
        if (!ok) return false;
        if (ferrite_model == JA) ok = ck_read_vec(f, M_irr) && ck_read_vec(f, M_field);
        else if (ferrite_model == LLG) ok = ck_read_vec(f, M_llg);
        else if (ferrite_model == Hybrid)
            ok = ck_read_vec(f, hyb_M_irr) && ck_read_vec(f, hyb_M_field) && ck_read_vec(f, hyb_M) && ck_read_vec(f, last_HzJA);
        else if (ferrite_model == Preisach) {
            ok = ck_read_vec(f, pr_M);
            for (int fi = 0; ok && fi < n_fer; ++fi) ok = ck_read_vec(f, pr_state[fi]);
        } else if (ferrite_model == Preisach2D) {
            ok = ck_read_vec(f, pr_M);
            for (int fi = 0; ok && fi < n_fer; ++fi)
                ok = (pr2d_algo == PR2D_BITS) ? ck_read_vec(f, pr2d_bits[fi]) : ck_read_vec(f, pr2d_rows[fi]);
        }
        return ok;
    }
};

// ---------------- C API для ctypes ----------------
#ifndef MENDRIVE_NO_C_API
extern "C" {

// Сигнатура совпадает с оригиналом (mendrive_ctypes.py работает без изменений).
// Preisach2D здесь использует PR2D_BITS (те же числа, что оригинальный алгоритм).
void* mendrive_create(int N, double a, double h_l, double h_r, double dt_frac,
                       int excitation_mode, int ferrite_model, int bias_orientation, double H0_bias,
                       double sigma_m_leak, double sigma_e_left, double sigma_e_ferrite,
                       double Ms, double a_JA, double alpha_JA, double k_JA, double c_JA,
                       int n_sub, int n_fp, int n_newton,
                       double Ms_llg, double gamma_llg, double alpha_llg, double skin_depth,
                       int n_hyst_pr, double Hc_mean_pr, double Hc_sigma_pr, double Hb_sigma_pr,
                       double sigma_m_leak_pr) {
    // FIX: исключение из конструктора (в т.ч. std::bad_alloc при нехватке памяти)
    // не должно пересекать границу extern "C" -- иначе std::terminate и падение
    // процесса Python/ядра Jupyter. Возвращаем nullptr и печатаем причину.
    try {
        return new MenDriveSim(N, a, h_l, h_r, dt_frac, excitation_mode, ferrite_model,
                                bias_orientation, H0_bias, sigma_m_leak, sigma_e_left, sigma_e_ferrite,
                                Ms, a_JA, alpha_JA, k_JA, c_JA, n_sub, n_fp, n_newton,
                                Ms_llg, gamma_llg, alpha_llg, skin_depth,
                                n_hyst_pr, Hc_mean_pr, Hc_sigma_pr, Hb_sigma_pr, sigma_m_leak_pr,
                                PR2D_BITS, 1);
    } catch (const std::exception& e) {
        fprintf(stderr, "mendrive_create: %s\n", e.what());
        return nullptr;
    }
}
// Расширенный конструктор: выбор алгоритма Preisach2D (0=bits, 1=rows) и числа потоков.
void* mendrive_create_ex(int N, double a, double h_l, double h_r, double dt_frac,
                       int excitation_mode, int ferrite_model, int bias_orientation, double H0_bias,
                       double sigma_m_leak, double sigma_e_left, double sigma_e_ferrite,
                       double Ms, double a_JA, double alpha_JA, double k_JA, double c_JA,
                       int n_sub, int n_fp, int n_newton,
                       double Ms_llg, double gamma_llg, double alpha_llg, double skin_depth,
                       int n_hyst_pr, double Hc_mean_pr, double Hc_sigma_pr, double Hb_sigma_pr,
                       double sigma_m_leak_pr, int pr2d_algo, int n_threads) {
    try {
        return new MenDriveSim(N, a, h_l, h_r, dt_frac, excitation_mode, ferrite_model,
                                bias_orientation, H0_bias, sigma_m_leak, sigma_e_left, sigma_e_ferrite,
                                Ms, a_JA, alpha_JA, k_JA, c_JA, n_sub, n_fp, n_newton,
                                Ms_llg, gamma_llg, alpha_llg, skin_depth,
                                n_hyst_pr, Hc_mean_pr, Hc_sigma_pr, Hb_sigma_pr, sigma_m_leak_pr,
                                pr2d_algo, n_threads);
    } catch (const std::exception& e) {
        fprintf(stderr, "mendrive_create_ex: %s\n", e.what());
        return nullptr;
    }
}
void mendrive_destroy(void* handle) { delete static_cast<MenDriveSim*>(handle); }
int mendrive_get_n_fer(void* handle) { return handle ? static_cast<MenDriveSim*>(handle)->n_fer : 0; }
double mendrive_get_dt(void* handle) { return handle ? (double)static_cast<MenDriveSim*>(handle)->dt : 0.0; }
double mendrive_get_dx(void* handle) { return handle ? (double)static_cast<MenDriveSim*>(handle)->dx : 0.0; }

// Прогон nsteps шагов, запись с record_start. probe_idx -- индекс внутри
// феррита для скалярных Hn/M/HzJA/MzJA. Буферы out_* должны быть выделены
// вызывающей стороной размером >= (nsteps - record_start).
// NOTE: интерфейс сохранён ради совместимости; автономная программа main.cpp
// эту функцию НЕ использует (она пишет результаты потоково на диск).
int mendrive_run(void* handle, double omega0, int nsteps, int record_start,
                  double amp, double ramp_periods, int probe_idx,
                  double* out_t, double* out_dTxx, double* out_Hn, double* out_P,
                  double* out_Mx, double* out_My, double* out_Mz,
                  double* out_HzJA, double* out_MzJA, double* out_Pdiss, int* out_blew_up) {
    MenDriveSim* sim = static_cast<MenDriveSim*>(handle);
    if (!sim) { if (out_blew_up) *out_blew_up = 1; return 0; }   // FIX: nullptr-handle
    double T = 2*M_PI/omega0;
    int rec_count = 0;
    *out_blew_up = 0;
    // FIX: в оригинале probe_idx не проверялся -> чтение за границей M_field/pr_M/...
    // при probe_idx >= n_fer (UB, возможное падение). Теперь ограничиваем.
    if (probe_idx < 0 || probe_idx >= sim->n_fer) {
        fprintf(stderr, "mendrive_run: probe_idx=%d вне [0,%d), использую 0\n", probe_idx, sim->n_fer);
        probe_idx = 0;
    }
    const bool have_probe = sim->n_fer > 0;
    const int gi = sim->probe_global_index(probe_idx);  // PERF: вычисляется один раз
    for (int n = 0; n < nsteps; ++n) {
        double P_step = sim->step_once(omega0, amp, ramp_periods, T);
        if (n >= record_start) {
            out_t[rec_count] = sim->t;
            // NOTE (ОШИБКА МАСШТАБА, оставлена ради совместимости с ноутбуком):
            // P_step -- уже мощность; деление на dt завышает P в 1/dt раз
            // (при N=160, dt=5e-4 -- в 2000 раз). main.cpp выводит обе величины:
            // P (как здесь) и P_src = P_step.
            out_P[rec_count] = P_step / sim->dt;
            out_dTxx[rec_count] = sim->Txx_at(sim->iR_A) - sim->Txx_at(sim->iL_A);
            out_Pdiss[rec_count] = (double)sim->dissipated_power_density();
            out_Hn[rec_count] = (gi >= 0) ? (double)sim->Hz[gi] : 0.0;
            if (!have_probe) {
                out_Mx[rec_count] = 0.0; out_My[rec_count] = 0.0; out_Mz[rec_count] = 0.0;
                out_HzJA[rec_count] = 0.0; out_MzJA[rec_count] = 0.0;
            } else if (sim->ferrite_model == JA) {
                out_My[rec_count] = sim->M_field[probe_idx];
                out_Mx[rec_count] = 0.0; out_Mz[rec_count] = 0.0;
                out_HzJA[rec_count] = 0.0; out_MzJA[rec_count] = 0.0;
            } else if (sim->ferrite_model == LLG) {
                Vec3 M = sim->M_llg[probe_idx];
                out_Mx[rec_count] = M.x; out_My[rec_count] = M.y; out_Mz[rec_count] = M.z;
                out_HzJA[rec_count] = 0.0; out_MzJA[rec_count] = 0.0;
            } else if (sim->ferrite_model == Hybrid) {
                Vec3 M = sim->hyb_M[probe_idx];
                out_Mx[rec_count] = M.x; out_My[rec_count] = M.y; out_Mz[rec_count] = M.z;
                out_HzJA[rec_count] = sim->last_HzJA[probe_idx];
                out_MzJA[rec_count] = sim->hyb_M_field[probe_idx];
            } else { // Preisach Preisach2D
                out_My[rec_count] = (double)sim->pr_M[probe_idx];
                out_Mx[rec_count] = 0.0; out_Mz[rec_count] = 0.0;
                out_HzJA[rec_count] = 0.0; out_MzJA[rec_count] = 0.0;
            }
            rec_count++;
        }
        if (!sim->fields_finite()) { *out_blew_up = 1; break; }
    }
    return rec_count;
}

} // extern "C"
#endif // MENDRIVE_NO_C_API
