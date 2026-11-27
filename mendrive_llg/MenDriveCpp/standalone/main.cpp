// main.cpp -- автономная программа расчёта резонатора МенДрайв (без Python/ctypes/Jupyter).
//
// Порт логики main_cpp() / scan_resonance_cpp() / force_per_power() /
// check_energy_balance() / shoelace_area() из MenDriveCpp_run.ipynb и
// mendrive_ctypes.py. Ядро -- mendrive_core.cpp (включается целиком).
//
// Главные отличия от ноутбука (ради работы на ограниченном железе):
//   * результаты НЕ копятся в памяти: средние считаются на лету (суммы Кэхэна),
//     временной ряд пишется на диск потоково и с прореживанием;
//   * в памяти держится только последний период (для петли гистерезиса) и
//     аккумуляторы ДПФ на выбранных частотах (вместо полного FFT по 29 млн точек);
//   * периодические чекпойнты + продолжение прогона ключом --resume;
//   * прогресс в stdout; Ctrl+C / SIGTERM -> чекпойнт и аккуратный выход;
//   * многопоточность по ячейкам феррита (threads=K), детерминированная.
//
// Запуск:  ./mendrive [-c файл.cfg] [ключ=значение ...] [--resume]
//          ./mendrive --help        -- список всех параметров

#define MENDRIVE_NO_C_API
#include "mendrive_core.cpp"

#include <cerrno>
#include <cinttypes>
#include <complex>
#include <csignal>
#include <chrono>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <map>
#include <sstream>
#include <string>
#include <vector>
#include <sys/resource.h>

namespace fs = std::filesystem;

static volatile std::sig_atomic_t g_stop = 0;
static void on_signal(int) { g_stop = 1; }

static const double NaN = std::numeric_limits<double>::quiet_NaN();

// ============================================================================
// Параметры
// ============================================================================
// kind = 0: влияет на результат (при --resume обязан совпадать);
// kind = 1: служебный (можно менять между запусками);
// kind = 2: MODE: новый параметр, влияющий на результат; входит в подпись чекпойнта, только если
//           отличается от значения по умолчанию (старые чекпойнты продолжаются с --resume);
// kind = 3: FIELDS: экспорт полей (на динамику не влияет; своя подпись в checkpoint_fields.bin).
struct ParamDef { const char* key; const char* def; int kind; const char* help; };

static const ParamDef PARAMS[] = {
    // --- геометрия и сетка
    {"N",               "20",       0, "число ячеек сетки (в ноутбуке обычно 160)"},
    {"a",               "0.2",      0, "полуширина вакуумного промежутка"},
    {"h_l",             "0.2",      0, "толщина левой (проводящей) стенки"},
    {"h_r",             "0.2",      0, "толщина правой (ферритовой) стенки"},
    {"dt_frac",         "0.1",      0, "dt = dt_frac*dx"},
    // --- возбуждение и модель феррита
    {"excitation_mode", "magnetic_right", 0, "magnetic_right(=jm_z_right) | electric_left(=je_y_left) | je_z_left | jm_y_right | je_y_right | je_z_right | jm_z_left | jm_y_left (README §8)"},
    {"ferrite_model",   "Preisach2D", 0, "JA | LLG | Hybrid | Preisach | Preisach2D | Linear (B=mu*H, эталон для сверки)"},
    // --- MODE: поляризация и линейные параметры стенок (по умолчанию -- как в исходном ядре)
    {"polarization",    "auto",     2, "auto | TE (Ey,Hz) | TM (Ez,Hy) | both; auto = пара источника; LLG/Hybrid -> всегда both"},
    {"src_k_phase",     "0",        2, "Re K источника: j=amp*exp(-xi/skin_depth)*sin(w t - k*xi) (профиль exp(iK xi) ноутбука)"},
    {"src_sigma_fix",   "0",        2, "1 = je делится на (eps+2*pi*sigma*dt) (согласованная полунеявная схема; 0 = как в исходном ядре)"},
    {"eps_left",        "1",        2, "диэлектрическая проницаемость левой стенки"},
    {"eps_ferrite",     "1",        2, "диэлектрическая проницаемость феррита (в ноутбуках epsilon_r = 16)"},
    {"mu_yy_lin",       "1",        2, "mu для Hy в феррите (JA/Preisach/Preisach2D/Linear; в ноутбуках mu_r_yy)"},
    {"mu_zz_lin",       "1",        2, "mu для Hz в феррите (только Linear; в ноутбуках mu_r_zz)"},
    {"bias_orientation","x",        0, "none | x | y | z | parallel"},
    {"H0_bias",         "4.0",      0, "поле подмагничивания"},
    {"sigma_m_leak",    "3.0",      0, "магнитная утечка JA/Hybrid"},
    {"sigma_e_left",    "3.0",      0, "проводимость левой стенки"},
    {"sigma_e_ferrite", "0.0",      0, "электропроводность феррита"},
    {"Ms",              "1.0",      0, "намагниченность насыщения (JA/Preisach)"},
    {"a_JA",            "0.3",      0, "JA: a"},
    {"alpha_JA",        "0.001",    0, "JA: alpha"},
    {"k_JA",            "0.15",     0, "JA: k"},
    {"c_JA",            "0.15",     0, "JA: c"},
    {"n_sub",           "2",        0, "JA: подшаги RK4"},
    {"n_fp",            "2",        0, "JA: итерации неподвижной точки"},
    {"n_newton",        "4",        0, "итерации Ньютона неявного шага"},
    {"Ms_llg",          "0.3",      0, "LLG: Ms"},
    {"gamma_llg",       "1.0",      0, "LLG: gamma"},
    {"alpha_llg",       "0.1",      0, "LLG: alpha"},
    {"skin_depth",      "0.06",     0, "глубина источника в стенке"},
    {"n_hyst_pr",       "14",       0, "Preisach: число гистеронов; Preisach2D: n_grid (если n_grid=0)"},
    {"n_grid",          "0",        0, "Preisach2D: узлов сетки на ось (>0 переопределяет n_hyst_pr)"},
    {"Hc_mean_pr",      "1.0",      0, "Preisach: средняя коэрцитивность"},
    {"Hc_sigma_pr",     "0.3",      0, "Preisach: разброс Hc"},
    {"Hb_sigma_pr",     "0.5",      0, "Preisach: разброс поля взаимодействия"},
    {"sigma_m_leak_pr", "3.0",      0, "Preisach: магнитная утечка (ноутбук НЕ передаёт sigma_m_leak сюда!)"},
    {"pr2d_algo",       "rows",     0, "Preisach2D: rows (быстро, O(n_grid)) | bits (как оригинал, O(n_grid^2))"},
    // --- скан резонанса (если omega0 <= 0)
    {"omega0",          "0",        0, ">0: частота задана, скан пропускается"},
    {"scan_lo",         "0.005",    0, "грубый скан: np.arange(scan_lo, scan_hi, scan_step)"},
    {"scan_hi",         "2.005",    0, ""},
    {"scan_step",       "0.005",    0, ""},
    {"fine_halfwidth",  "0.9",      0, "уточняющий скан: +-полуширина вокруг пика"},
    {"fine_step",       "0.001",    0, "шаг уточняющего скана"},
    {"scan_periods",    "3",        0, "периодов на точку скана"},
    {"scan_record_from","1",        0, "запись отклика с периода"},
    {"scan_amp",        "0.3",      0, "амплитуда в скане"},
    {"scan_ramp",       "1.0",      0, "разгон в скане, периодов"},
    // --- длинный прогон
    {"n_periods_total", "80",       0, "всего периодов"},
    {"record_from_period","15",     0, "запись с периода"},
    {"amp",             "1.0",      0, "амплитуда источника"},
    {"ramp_periods",    "2.0",      0, "плавный разгон, периодов"},
    {"probe_idx",       "0",        0, "индекс ячейки-пробы внутри феррита"},
    {"last_frac",       "0.5",      0, "доля записанных периодов для средних (как last_frac_force)"},
    {"length_unit_cm",  "1.0",      0, "калибровка force_per_kW (как в mendrive_ctypes.py)"},
    {"area_cm2",        "1.0",      0, ""},
    {"time_unit_s",     "1.0",      0, ""},
    // --- вывод
    {"out_dir",         "out",      1, "папка результатов"},
    {"ts_every",        "10",       0, "временной ряд: писать каждый k-й шаг (0 = не писать)"},
    {"ts_last_periods", "20",       0, "временной ряд: только последние K периодов (0 = вся запись)"},
    {"ts_format",       "csv",      0, "csv | bin (float64, столбцы в ts_columns.txt)"},
    {"fft_harmonics",   "8",        0, "ДПФ dTxx на гармониках k*omega0, k=0..K"},
    {"fft_bins",        "400",      0, "ДПФ на равномерной сетке [0, fft_wmax_factor*omega0]"},
    {"fft_wmax_factor", "4",        0, ""},
    // --- FIELDS: экспорт профилей полей для сверки с линейными моделями (README §9)
    {"fields_snap_phases", "",      3, "снимки профилей в фазах периода, напр. 0,0.25,0.5,0.75 (пусто = нет)"},
    {"fields_snap_every_periods","0", 3, "снимки каждые K периодов начиная с fields_snap_from_period (0 = только последний период)"},
    {"fields_snap_from_period","-1", 3, "с какого периода снимки (-1 = record_from_period)"},
    {"fields_harm",     "",         3, "номера гармоник для комплексных профилей, напр. 1 или 0,1,2,3 (пусто = нет)"},
    {"fields_harm_periods","0",     3, "накопление по последним K периодам (0 = то же окно, что для средних: last_frac)"},
    {"time_convention", "minus",    3, "minus: F=Re[F~ exp(-i w t)] (как в ноутбуках MenDrive_qnm*) | plus: exp(+i w t)"},
    {"fields_harm_normalize","none",3, "none | source: разделить на комплексную амплитуду источника j0 (= i*amp при minus)"},
    {"fields_every_x",  "1",        3, "прореживание узлов по x в CSV полей"},
    // --- служебные
    {"threads",         "1",        1, "потоков (ячейки феррита в прогоне; частоты в скане)"},
    {"checkpoint_every_periods","100", 1, "чекпойнт каждые K периодов"},
    {"checkpoint_every_sec","1800", 1, "...и не реже, чем раз в S секунд"},
    {"progress_every_periods","0",  1, "печать прогресса каждые K периодов (0 = авто, ~1%)"},
    {"force_resume",    "0",        1, "1 = продолжить, даже если конфиг отличается от чекпойнта"},
};

static std::string trim(const std::string& s) {
    size_t b = s.find_first_not_of(" \t\r\n"), e = s.find_last_not_of(" \t\r\n");
    return (b == std::string::npos) ? std::string() : s.substr(b, e - b + 1);
}

struct Config {
    std::map<std::string, std::string> kv;
    Config() { for (const auto& p : PARAMS) kv[p.key] = p.def; }
    static const ParamDef* def_of(const std::string& k) {
        for (const auto& p : PARAMS) if (k == p.key) return &p;
        return nullptr;
    }
    void set(const std::string& k, const std::string& v) {
        if (!def_of(k)) throw std::runtime_error("неизвестный параметр '" + k + "' (см. --help)");
        kv[k] = v;
    }
    const std::string& s(const std::string& k) const {
        auto it = kv.find(k);
        if (it == kv.end()) throw std::runtime_error("нет параметра " + k);
        return it->second;
    }
    double d(const std::string& k) const {
        const std::string& v = s(k);
        char* end = nullptr; errno = 0;
        double x = std::strtod(v.c_str(), &end);
        if (end == v.c_str() || *end != '\0' || errno == ERANGE)
            throw std::runtime_error("параметр " + k + ": не число '" + v + "'");
        return x;
    }
    long long i(const std::string& k) const {
        double x = d(k);
        if (x != std::floor(x) || std::fabs(x) > 9e15)
            throw std::runtime_error("параметр " + k + ": ожидалось целое, получено '" + s(k) + "'");
        return (long long)x;
    }
    void load_file(const std::string& path) {
        std::ifstream in(path);
        if (!in) throw std::runtime_error("не удалось открыть конфиг " + path);
        std::string line; int ln = 0;
        while (std::getline(in, line)) {
            ++ln;
            size_t h = line.find('#'); if (h != std::string::npos) line = line.substr(0, h);
            line = trim(line);
            if (line.empty()) continue;
            size_t eq = line.find('=');
            if (eq == std::string::npos)
                throw std::runtime_error(path + ":" + std::to_string(ln) + ": ожидалось ключ=значение");
            set(trim(line.substr(0, eq)), trim(line.substr(eq + 1)));
        }
    }
    std::string signature() const {
        std::ostringstream o;
        // MODE: kind 2 -- только если отличается от умолчания (совместимость старых чекпойнтов)
        for (const auto& p : PARAMS)
            if (p.kind == 0 || (p.kind == 2 && s(p.key) != p.def)) o << p.key << '=' << s(p.key) << '\n';
        return o.str();
    }
    // FIELDS: подпись параметров экспорта (kind 3) -- хранится в checkpoint_fields.bin
    std::string fields_signature() const {
        std::ostringstream o;
        for (const auto& p : PARAMS) if (p.kind == 3) o << p.key << '=' << s(p.key) << '\n';
        return o.str();
    }
    void dump(FILE* f) const { for (const auto& p : PARAMS) fprintf(f, "%s=%s\n", p.key, s(p.key).c_str()); }
};

static int code_of(const std::string& v, const std::vector<std::pair<std::string, int>>& m, const char* what) {
    for (const auto& p : m) if (p.first == v) return p.second;
    std::string opts;
    for (const auto& p : m) { if (!opts.empty()) opts += ", "; opts += p.first; }
    throw std::runtime_error(std::string(what) + ": неизвестное значение '" + v + "' (допустимо: " + opts + ")");
}
static int ferrite_code(const Config& c) {
    return code_of(c.s("ferrite_model"), {{"JA", JA}, {"LLG", LLG}, {"Hybrid", Hybrid},
                   {"Preisach", Preisach}, {"Preisach2D", Preisach2D}, {"Linear", Linear}}, "ferrite_model");
}
static int excitation_code(const Config& c) {
    // MODE: прежние имена + новые (таблица соответствия driving_mode -- README §8)
    return code_of(c.s("excitation_mode"), {{"magnetic_right", EXC_JM_Z_RIGHT}, {"electric_left", EXC_JE_Y_LEFT},
                   {"jm_z_right", EXC_JM_Z_RIGHT}, {"je_y_left", EXC_JE_Y_LEFT}, {"je_z_left", EXC_JE_Z_LEFT},
                   {"jm_y_right", EXC_JM_Y_RIGHT}, {"je_y_right", EXC_JE_Y_RIGHT}, {"je_z_right", EXC_JE_Z_RIGHT},
                   {"jm_z_left", EXC_JM_Z_LEFT}, {"jm_y_left", EXC_JM_Y_LEFT}}, "excitation_mode");
}
static int polarization_code(const Config& c) {
    return code_of(c.s("polarization"), {{"auto", POL_AUTO}, {"TE", POL_TE}, {"TM", POL_TM}, {"both", POL_BOTH}}, "polarization");
}
static int bias_code(const Config& c) {
    return code_of(c.s("bias_orientation"), {{"none", 0}, {"x", 1}, {"y", 2}, {"z", 3}, {"parallel", 4}}, "bias_orientation");
}
static int pr2d_algo_code(const Config& c) {
    return code_of(c.s("pr2d_algo"), {{"bits", PR2D_BITS}, {"rows", PR2D_ROWS}}, "pr2d_algo");
}
static int effective_n_hyst(const Config& c) {
    long long n = c.i("n_hyst_pr");
    if (ferrite_code(c) == Preisach2D && c.i("n_grid") > 0) n = c.i("n_grid");
    if (n < 1 || n > 1000000) throw std::runtime_error("n_hyst_pr/n_grid вне диапазона [1, 1e6]");
    return (int)n;
}

static std::unique_ptr<MenDriveSim> make_sim(const Config& c, int threads) {
    ModeOptions mo;   // MODE
    mo.polarization = polarization_code(c);
    mo.src_k_phase = c.d("src_k_phase");
    mo.src_sigma_fix = (int)c.i("src_sigma_fix");
    mo.eps_left = c.d("eps_left"); mo.eps_ferrite = c.d("eps_ferrite");
    mo.mu_yy_lin = c.d("mu_yy_lin"); mo.mu_zz_lin = c.d("mu_zz_lin");
    return std::unique_ptr<MenDriveSim>(new MenDriveSim(
        (int)c.i("N"), c.d("a"), c.d("h_l"), c.d("h_r"), c.d("dt_frac"),
        excitation_code(c), ferrite_code(c), bias_code(c), c.d("H0_bias"),
        c.d("sigma_m_leak"), c.d("sigma_e_left"), c.d("sigma_e_ferrite"),
        c.d("Ms"), c.d("a_JA"), c.d("alpha_JA"), c.d("k_JA"), c.d("c_JA"),
        (int)c.i("n_sub"), (int)c.i("n_fp"), (int)c.i("n_newton"),
        c.d("Ms_llg"), c.d("gamma_llg"), c.d("alpha_llg"), c.d("skin_depth"),
        effective_n_hyst(c), c.d("Hc_mean_pr"), c.d("Hc_sigma_pr"), c.d("Hb_sigma_pr"),
        c.d("sigma_m_leak_pr"), pr2d_algo_code(c), threads, mo));
}

// MODE: тангенциальное H «своей» пары в ячейке-пробе: Hz для TE (и both -- как раньше), Hy для TM.
static inline double probe_H(const MenDriveSim& sim, int gi) {
    if (gi < 0) return 0.0;
    return (sim.upd_TE) ? (double)sim.Hz[gi] : (double)sim.Hy[gi];
}

// ============================================================================
// Утилиты
// ============================================================================
struct Kahan {   // компенсированная сумма (точность ~ np.mean с попарным суммированием)
    double s = 0.0, c = 0.0;
    void add(double x) { double y = x - c; double t = s + y; c = (t - s) - y; s = t; }
};

static double peak_rss_mb() {
    struct rusage ru; getrusage(RUSAGE_SELF, &ru);
#ifdef __APPLE__
    return ru.ru_maxrss / (1024.0 * 1024.0);   // macOS: байты
#else
    return ru.ru_maxrss / 1024.0;              // Linux: килобайты
#endif
}

static std::string fmt_dur(double sec) {
    if (!(sec >= 0) || !std::isfinite(sec)) return "?";
    long long s = (long long)(sec + 0.5);
    char buf[64];
    if (s >= 3600) snprintf(buf, sizeof buf, "%lldч %02lldм %02lldс", s / 3600, (s / 60) % 60, s % 60);
    else if (s >= 60) snprintf(buf, sizeof buf, "%lldм %02lldс", s / 60, s % 60);
    else snprintf(buf, sizeof buf, "%lldс", s);
    return buf;
}

static std::string fmt_bytes(double b) {
    char buf[64];
    if (b >= 1024.0 * 1024 * 1024) snprintf(buf, sizeof buf, "%.2f ГБ", b / (1024.0 * 1024 * 1024));
    else if (b >= 1024.0 * 1024) snprintf(buf, sizeof buf, "%.1f МБ", b / (1024.0 * 1024));
    else if (b >= 1024.0) snprintf(buf, sizeof buf, "%.1f КБ", b / 1024.0);
    else snprintf(buf, sizeof buf, "%.0f Б", b);
    return buf;
}

// np.arange(lo, hi, step): длина ceil((hi-lo)/step), элементы lo + k*step
static std::vector<double> arange(double lo, double hi, double step) {
    std::vector<double> v;
    if (!(step > 0)) return v;
    double n = std::ceil((hi - lo) / step);
    if (!(n > 0)) return v;
    for (long long k = 0; k < (long long)n; ++k) v.push_back(lo + k * step);
    return v;
}

static long long nanargmax(const std::vector<double>& v) {   // первый максимум, NaN пропускаются
    long long best = -1;
    for (size_t i = 0; i < v.size(); ++i)
        if (!std::isnan(v[i]) && (best < 0 || v[i] > v[(size_t)best])) best = (long long)i;
    return best;
}

// int(n_periods*T/dt) как в mendrive_ctypes.py
static int64_t n_steps_for(double n_periods, double T, double dt) { return (int64_t)(n_periods * T / dt); }

struct ProbeVals { double Hn, Mx, My, Mz, HzJA, MzJA; };
// Значения пробы -- в точности как пишет mendrive_run().
static ProbeVals probe_values(const MenDriveSim& sim, int probe_idx, int gi) {
    ProbeVals v{0, 0, 0, 0, 0, 0};
    v.Hn = probe_H(sim, gi);   // MODE: для TM-режимов -- Hy (для TE/both -- Hz, как раньше)
    if (sim.n_fer <= 0) return v;
    if (sim.ferrite_model == Linear) return v;   // MODE: у Linear нет состояния M
    if (sim.ferrite_model == JA) v.My = (double)sim.M_field[probe_idx];
    else if (sim.ferrite_model == LLG) {
        const Vec3& M = sim.M_llg[probe_idx]; v.Mx = (double)M.x; v.My = (double)M.y; v.Mz = (double)M.z;
    } else if (sim.ferrite_model == Hybrid) {
        const Vec3& M = sim.hyb_M[probe_idx]; v.Mx = (double)M.x; v.My = (double)M.y; v.Mz = (double)M.z;
        v.HzJA = (double)sim.last_HzJA[probe_idx]; v.MzJA = (double)sim.hyb_M_field[probe_idx];
    } else v.My = (double)sim.pr_M[probe_idx];
    return v;
}

// ============================================================================
// Скан резонанса (порт scan_resonance_cpp)
// ============================================================================
// Отклик = (max Hn - min Hn)/2 по записанной части; NaN при «взрыве» или пустой записи.
// Как и в ноутбуке, проба в скане -- ячейка 0 феррита (probe_idx не передаётся).
static double scan_response(const Config& c, double w, bool& stopped) {
    auto sim = make_sim(c, 1);
    const double dt = (double)sim->dt, T = 2 * M_PI / w;
    const int64_t nsteps = n_steps_for((double)c.d("scan_periods"), T, dt);
    const int64_t rs = n_steps_for((double)c.d("scan_record_from"), T, dt);
    if (nsteps - rs <= 0) return NaN;
    const int gi = sim->probe_global_index(0);
    const double amp = c.d("scan_amp"), ramp = c.d("scan_ramp");
    double mx = -std::numeric_limits<double>::infinity(), mn = std::numeric_limits<double>::infinity();
    int64_t nrec = 0; bool blew = false;
    for (int64_t n = 0; n < nsteps; ++n) {
        sim->step_once(w, amp, ramp, T);
        if (n >= rs) {
            double h = probe_H(*sim, gi);   // MODE: Hz (TE/both) или Hy (TM)
            mx = std::max(mx, h); mn = std::min(mn, h); ++nrec;
        }
        if (!sim->fields_finite()) { blew = true; break; }
        if ((n & 1023) == 0 && g_stop) { stopped = true; return NaN; }
    }
    if (blew || nrec == 0) return NaN;
    return (mx - mn) / 2;
}

// Скан набора частот; результаты построчно дописываются в scan_<name>.csv
// (сразу после каждой частоты), поэтому прерванный скан продолжается с места.
static bool run_scan(const Config& c, const std::string& name, const std::vector<double>& freqs,
                     int threads, const fs::path& out_dir, bool reuse, std::vector<double>& resp) {
    resp.assign(freqs.size(), NaN);
    std::vector<char> done(freqs.size(), 0);
    const fs::path path = out_dir / ("scan_" + name + ".csv");
    size_t n_reused = 0;
    const bool exists = fs::exists(path);
    if (reuse && exists) {
        std::ifstream in(path.string());
        std::string line;
        while (std::getline(in, line)) {
            if (line.empty() || line[0] == 'i') continue;
            std::istringstream ss(line);
            std::string a, b, r;
            if (!std::getline(ss, a, ',') || !std::getline(ss, b, ',') || !std::getline(ss, r, ',')) continue;
            size_t idx = (size_t)std::strtoull(a.c_str(), nullptr, 10);
            double w = std::strtod(b.c_str(), nullptr), v = std::strtod(r.c_str(), nullptr);
            if (idx < freqs.size() && std::fabs(w - freqs[idx]) <= 1e-9 * std::max(1.0, std::fabs(freqs[idx])) && !done[idx]) {
                resp[idx] = v; done[idx] = 1; ++n_reused;
            }
        }
    }
    FILE* f = fopen(path.string().c_str(), (reuse && exists) ? "ab" : "wb");
    if (!f) throw std::runtime_error("не удалось открыть " + path.string());
    if (!(reuse && exists)) { fprintf(f, "idx,omega,resp\n"); fflush(f); }

    std::vector<size_t> pending;
    for (size_t i = 0; i < freqs.size(); ++i) if (!done[i]) pending.push_back(i);
    // оценка объёма работы: сумма шагов по всем частотам
    double dt = 0.0;
    { auto probe = make_sim(c, 1); dt = (double)probe->dt; }
    double total_steps = 0;
    for (size_t i : pending) total_steps += (double)n_steps_for(c.d("scan_periods"), 2 * M_PI / freqs[i], dt);
    printf("[scan %s] частот: %zu (уже посчитано: %zu), шагов FDTD к выполнению: %.3g, потоков: %d\n",
           name.c_str(), freqs.size(), n_reused, total_steps, threads);

    std::atomic<size_t> next{0};
    std::mutex mtx;
    size_t n_now = 0;
    const auto t0 = std::chrono::steady_clock::now();
    const size_t print_every = std::max<size_t>(1, pending.size() / 40);
    auto worker = [&]() {
        for (;;) {
            if (g_stop) return;
            size_t q = next.fetch_add(1);
            if (q >= pending.size()) return;
            size_t idx = pending[q];
            bool st = false;
            double r = scan_response(c, freqs[idx], st);
            if (st) return;
            std::lock_guard<std::mutex> lk(mtx);
            resp[idx] = r; done[idx] = 1;
            fprintf(f, "%zu,%.17g,%.17g\n", idx, freqs[idx], r); fflush(f);
            ++n_now;
            if (n_now % print_every == 0 || n_now == pending.size()) {
                double el = std::chrono::duration<double>(std::chrono::steady_clock::now() - t0).count();
                printf("[scan %s] %zu/%zu  omega=%.6g  отклик=%.6g  прошло %s\n", name.c_str(), n_now,
                       pending.size(), freqs[idx], r, fmt_dur(el).c_str());
            }
        }
    };
    std::vector<std::thread> th;
    for (int t = 1; t < threads; ++t) th.emplace_back(worker);
    worker();
    for (auto& x : th) x.join();
    fclose(f);
    return !g_stop;
}

// ============================================================================
// Длинный прогон: потоковая статистика
// ============================================================================
struct RunState {
    // план (однозначно задаётся конфигом и omega)
    double omega = 0, T = 0, dt = 0;
    int64_t nsteps = 0, record_start = 0, n_rec = 0, spp = 1, n_use = 0, win_start = 0, ts_from = 0;
    // прогресс
    int64_t n_done = 0; int32_t blew_up = 0; int64_t rec_count = 0; double wall_prev = 0;
    // окно усреднения (последние n_use записанных шагов, как в force_per_power/check_energy_balance)
    Kahan w_dTxx, w_P, w_Psrc, w_Pdiss, w_Pe, w_Pm, w_Pmag, w_Pnum; int64_t w_count = 0;
    double w_Hmin = std::numeric_limits<double>::infinity(), w_Hmax = -std::numeric_limits<double>::infinity();
    // текущий блок из spp записанных шагов (период) -- для таблицы сходимости
    Kahan b_dTxx, b_P, b_Psrc, b_Pdiss, b_Pnum; int64_t b_count = 0, blocks_done = 0;
    double b_Hmin = std::numeric_limits<double>::infinity(), b_Hmax = -std::numeric_limits<double>::infinity();
    Kahan cum_F, cum_P, cum_Psrc;
    // кольцевой буфер последнего периода (петля гистерезиса)
    std::vector<double> ring_H, ring_M; int64_t ring_pos = 0, ring_filled = 0;
    // ДПФ dTxx с окном Ханна по всей записи, на выбранных частотах
    std::vector<double> fft_w, fft_re, fft_im, fft_zr, fft_zi;
    // размеры выходных файлов на момент чекпойнта
    uint64_t per_bytes = 0, ts_bytes = 0;
};

static void plan_run(const Config& c, double omega, double dt, RunState& R) {
    R.omega = omega; R.dt = dt; R.T = 2 * M_PI / omega;
    R.nsteps = n_steps_for(c.d("n_periods_total"), R.T, dt);
    R.record_start = n_steps_for(c.d("record_from_period"), R.T, dt);
    R.n_rec = std::max<int64_t>(0, R.nsteps - R.record_start);
    R.spp = std::max<int64_t>(1, (int64_t)std::nearbyint(R.T / dt));      // = round(T/dt)
    const int64_t n_periods = std::max<int64_t>(1, R.n_rec / R.spp);
    R.n_use = std::max<int64_t>(1, (int64_t)((double)n_periods * c.d("last_frac"))) * R.spp;
    R.n_use = std::min(R.n_use, R.n_rec);
    R.win_start = R.nsteps - R.n_use;
    const long long tl = c.i("ts_last_periods");
    R.ts_from = (tl > 0) ? std::max<int64_t>(R.record_start, R.nsteps - (int64_t)tl * R.spp) : R.record_start;
    R.ring_H.assign((size_t)R.spp, 0.0); R.ring_M.assign((size_t)R.spp, 0.0);
    R.fft_w.clear();
    const long long K = c.i("fft_harmonics");
    for (long long k = 0; k <= K; ++k) R.fft_w.push_back((double)k * omega);
    const long long B = c.i("fft_bins");
    const double wmax = c.d("fft_wmax_factor") * omega;
    if (B >= 2) for (long long m = 0; m < B; ++m) R.fft_w.push_back(wmax * (double)m / (double)(B - 1));
    R.fft_re.assign(R.fft_w.size(), 0.0); R.fft_im.assign(R.fft_w.size(), 0.0);
    R.fft_zr.assign(R.fft_w.size(), 1.0); R.fft_zi.assign(R.fft_w.size(), 0.0);
}

template <class F> static bool io_run(RunState& R, F&& P) {
    bool ok = true;
    ok = ok && P(R.omega) && P(R.T) && P(R.dt) && P(R.nsteps) && P(R.record_start) && P(R.n_rec)
            && P(R.spp) && P(R.n_use) && P(R.win_start) && P(R.ts_from);
    ok = ok && P(R.n_done) && P(R.blew_up) && P(R.rec_count) && P(R.wall_prev);
    ok = ok && P(R.w_dTxx) && P(R.w_P) && P(R.w_Psrc) && P(R.w_Pdiss) && P(R.w_Pe) && P(R.w_Pm)
            && P(R.w_Pmag) && P(R.w_Pnum) && P(R.w_count) && P(R.w_Hmin) && P(R.w_Hmax);
    ok = ok && P(R.b_dTxx) && P(R.b_P) && P(R.b_Psrc) && P(R.b_Pdiss) && P(R.b_Pnum) && P(R.b_count)
            && P(R.blocks_done) && P(R.b_Hmin) && P(R.b_Hmax) && P(R.cum_F) && P(R.cum_P) && P(R.cum_Psrc);
    ok = ok && P(R.ring_pos) && P(R.ring_filled) && P(R.per_bytes) && P(R.ts_bytes);
    return ok;
}
static bool save_run(FILE* f, RunState& R) {
    bool ok = io_run(R, [&](auto& x) { return ck_write_pod(f, x); });
    return ok && ck_write_vec(f, R.ring_H) && ck_write_vec(f, R.ring_M) && ck_write_vec(f, R.fft_w)
              && ck_write_vec(f, R.fft_re) && ck_write_vec(f, R.fft_im) && ck_write_vec(f, R.fft_zr)
              && ck_write_vec(f, R.fft_zi);
}
static bool load_run(FILE* f, RunState& R) {
    bool ok = io_run(R, [&](auto& x) { return ck_read_pod(f, x); });
    return ok && ck_read_vec(f, R.ring_H) && ck_read_vec(f, R.ring_M) && ck_read_vec(f, R.fft_w)
              && ck_read_vec(f, R.fft_re) && ck_read_vec(f, R.fft_im) && ck_read_vec(f, R.fft_zr)
              && ck_read_vec(f, R.fft_zi);
}

static const char CK_MAGIC[8] = {'M', 'D', 'C', 'K', 'P', 'T', '0', '1'};

static void write_checkpoint(const fs::path& path, const std::string& sig, RunState& R,
                             const MenDriveSim& sim, FILE* fper, FILE* fts) {
    if (fper) { fflush(fper); R.per_bytes = (uint64_t)ftell(fper); }
    if (fts)  { fflush(fts);  R.ts_bytes = (uint64_t)ftell(fts); }
    const fs::path tmp = path.string() + ".tmp";
    FILE* f = fopen(tmp.string().c_str(), "wb");
    if (!f) { fprintf(stderr, "ВНИМАНИЕ: не удалось записать чекпойнт %s\n", tmp.string().c_str()); return; }
    uint32_t ld = sizeof(long double);
    uint64_t n = sig.size();
    bool ok = fwrite(CK_MAGIC, 1, 8, f) == 8 && ck_write_pod(f, ld) && ck_write_pod(f, n)
              && fwrite(sig.data(), 1, sig.size(), f) == sig.size() && ck_write_pod(f, R.omega)
              && save_run(f, R) && sim.save_state(f);
    ok = (fclose(f) == 0) && ok;
    if (!ok) { fprintf(stderr, "ВНИМАНИЕ: ошибка записи чекпойнта (диск полон?)\n"); return; }
    std::error_code ec;
    fs::rename(tmp, path, ec);   // атомарная замена: старый чекпойнт цел до успешной записи нового
    if (ec) fprintf(stderr, "ВНИМАНИЕ: не удалось переименовать чекпойнт: %s\n", ec.message().c_str());
}

// Читает omega и проверяет подпись конфига; остальное -- read_checkpoint_rest().
static FILE* open_checkpoint(const fs::path& path, const Config& c, double& omega) {
    FILE* f = fopen(path.string().c_str(), "rb");
    if (!f) return nullptr;
    char mg[8]; uint32_t ld = 0; uint64_t n = 0;
    if (fread(mg, 1, 8, f) != 8 || std::memcmp(mg, CK_MAGIC, 8) != 0 || !ck_read_pod(f, ld) || ld != sizeof(long double)
        || !ck_read_pod(f, n) || n > (1u << 20)) {
        fclose(f); throw std::runtime_error("чекпойнт повреждён или с другой платформы: " + path.string());
    }
    std::string sig(n, '\0');
    if (fread(&sig[0], 1, n, f) != n || !ck_read_pod(f, omega)) { fclose(f); throw std::runtime_error("чекпойнт обрезан"); }
    const std::string cur = c.signature();
    if (sig != cur) {
        fprintf(stderr, "Конфигурация отличается от сохранённой в чекпойнте:\n");
        std::istringstream a(sig), b(cur); std::string la, lb;
        while (std::getline(a, la) && std::getline(b, lb)) if (la != lb) fprintf(stderr, "   было: %s   стало: %s\n", la.c_str(), lb.c_str());
        if (c.i("force_resume") != 1) {
            fclose(f);
            throw std::runtime_error("продолжение невозможно (верните прежние параметры или force_resume=1)");
        }
        fprintf(stderr, "force_resume=1: продолжаю, несмотря на различия\n");
    }
    return f;
}

// ============================================================================
// FIELDS: экспорт профилей полей по x для сверки с линейными моделями MenDrive_qnm*.ipynb
// ============================================================================
// Все величины выводятся на A-сетке xA[i] = -(a+h_l) + i*dx (узлы E), x -- в единицах кода с тем же
// началом, что в ноутбуках (центр вакуумного зазора; стенки x = -a и x = +a), x_cm = x*length_unit_cm.
// Величины H-типа (B-сетка) интерполируются в узлы: (F[i-1]+F[i])/2, на концах -- ближайшая ячейка.
// Истинные моменты отсчёта после шага (схема Йи): E -- t; H, M, B -- t + dt/2; токи j -- t - dt.
// Гармоники (накопление на лету, без хранения ряда):
//     F~_n(x) = (c_n/T_rec) * sum_k F(x, tau_k) exp(+s i n w tau_k) dt,  c_0 = 1, c_{n>0} = 2,
// s = +1 при time_convention=minus (F = Re[F~ exp(-i n w t)] -- конвенция ноутбуков), s = -1 при plus;
// tau_k -- истинное время отсчёта каждой величины (фазовые поправки dt/2 и dt учитываются точно).
// Квадратичные величины (w, Sx, Txx) считаются с H, центрированным по времени: (H^{n+1/2}+H^{n-1/2})/2.
//   w   = (eps E^2 + H.B)/(8 pi)          -- плотность энергии (линейное приближение)
//   Sx  = (Ey Hz - Ez Hy)/(4 pi)          -- поток Пойнтинга (c = 1)
//   Txx = (Ey^2+Ez^2+Hy^2+Hz^2)/(8 pi)    -- как Txx_at() ядра (= -sigma_xx, давление; в вакууме = w)
enum { FQ_EY, FQ_EZ, FQ_HY, FQ_HZ, FQ_MX, FQ_MY, FQ_MZ, FQ_BY, FQ_BZ, FQ_JEY, FQ_JEZ, FQ_JMY, FQ_JMZ,
       FQ_W, FQ_SX, FQ_TXX, FQ_COUNT };
static const char* FQ_NAMES[FQ_COUNT] = {"Ey", "Ez", "Hy", "Hz", "Mx", "My", "Mz", "By", "Bz",
                                          "je_y", "je_z", "jm_y", "jm_z", "w", "Sx", "Txx"};
static const int FQ_TCLASS[FQ_COUNT] = {0, 0, 1, 1, 1, 1, 1, 1, 1, 2, 2, 2, 2, 0, 0, 0};  // 0: t, 1: t+dt/2, 2: t-dt

static std::vector<double> parse_list(const std::string& v, const char* what) {
    std::vector<double> out;
    std::string cur;
    std::istringstream ss(v);
    while (std::getline(ss, cur, ',')) {
        cur = trim(cur);
        if (cur.empty()) continue;
        char* end = nullptr;
        double x = std::strtod(cur.c_str(), &end);
        if (end == cur.c_str() || *end != '\0') throw std::runtime_error(std::string(what) + ": не число '" + cur + "'");
        out.push_back(x);
    }
    return out;
}

struct FieldsExport {
    bool snap_on = false, harm_on = false, harm_valid = true;
    std::vector<double> phases;            // доли периода, [0,1)
    long long snap_every = 0, snap_from = 0, p_last = -1;
    std::vector<int> harm;                 // номера гармоник
    int64_t harm_start = 0, harm_count = 0;
    double sgn = 1.0;                      // s: +1 (minus), -1 (plus)
    bool norm_src = false;
    int every_x = 1;
    double omega = 0, T = 0, dt = 0, amp = 0, x_unit = 1, a_geom = 0;
    std::string fsig, exc_name, pol_name;
    std::vector<std::vector<double>> acc_re, acc_im;   // [h*FQ_COUNT + q][i]
    std::vector<std::vector<double>> q;                // текущие значения [q][i]
    std::vector<long double> Mx, My, Mz;               // буферы
    FILE* fsnap = nullptr;
    uint64_t snap_bytes = 0;
    fs::path snap_path, ck_path;

    void init(const Config& c, const MenDriveSim& sim, const RunState& R, const fs::path& out_dir) {
        phases = parse_list(c.s("fields_snap_phases"), "fields_snap_phases");
        for (double& p : phases) { p = p - std::floor(p); }
        std::sort(phases.begin(), phases.end());
        snap_on = !phases.empty();
        snap_every = c.i("fields_snap_every_periods");
        const double sf = c.d("fields_snap_from_period");
        snap_from = (sf < 0) ? (long long)std::floor(c.d("record_from_period")) : (long long)std::floor(sf);
        for (double h : parse_list(c.s("fields_harm"), "fields_harm")) {
            if (h < 0 || h != std::floor(h)) throw std::runtime_error("fields_harm: нужны целые >= 0");
            harm.push_back((int)h);
        }
        harm_on = !harm.empty();
        const std::string tc = c.s("time_convention");
        if (tc != "minus" && tc != "plus") throw std::runtime_error("time_convention: minus | plus");
        sgn = (tc == "minus") ? 1.0 : -1.0;
        const std::string nm = c.s("fields_harm_normalize");
        if (nm != "none" && nm != "source") throw std::runtime_error("fields_harm_normalize: none | source");
        norm_src = (nm == "source");
        every_x = (int)std::max(1LL, c.i("fields_every_x"));
        omega = R.omega; T = R.T; dt = R.dt; amp = c.d("amp"); x_unit = c.d("length_unit_cm"); a_geom = (double)sim.a;
        exc_name = c.s("excitation_mode"); pol_name = c.s("polarization");
        fsig = c.fields_signature();
        const long long hp = c.i("fields_harm_periods");
        harm_start = (hp > 0) ? std::max<int64_t>(0, R.nsteps - (int64_t)hp * R.spp) : R.win_start;
        p_last = (long long)std::floor((double)R.nsteps * dt / T + 1e-12) - 1;
        if (harm_on) {
            acc_re.assign(harm.size() * FQ_COUNT, std::vector<double>(sim.N + 1, 0.0));
            acc_im.assign(harm.size() * FQ_COUNT, std::vector<double>(sim.N + 1, 0.0));
            if ((double)harm_start * dt < c.d("ramp_periods") * T)
                printf("[fields] ВНИМАНИЕ: накопление гармоник начинается до окончания разгона (ramp_periods)\n");
        }
        q.assign(FQ_COUNT, std::vector<double>(sim.N + 1, 0.0));
        snap_path = out_dir / "fields_snap.csv";
        ck_path = out_dir / "checkpoint_fields.bin";
        if (snap_on || harm_on)
            printf("[fields] снимки: %s; гармоники: %s (окно с шага %lld, %.1f периодов), time_convention=%s, normalize=%s\n",
                   snap_on ? c.s("fields_snap_phases").c_str() : "нет", harm_on ? c.s("fields_harm").c_str() : "нет",
                   (long long)harm_start, (double)(R.nsteps - harm_start) / (double)R.spp, tc.c_str(), nm.c_str());
    }

    static void head_cols(FILE* f) { for (int k = 0; k < FQ_COUNT; ++k) fprintf(f, ",%s", FQ_NAMES[k]); }

    void open_snap(bool append) {
        if (!snap_on) return;
        if (append && fs::exists(snap_path)) {
            std::error_code ec; fs::resize_file(snap_path, snap_bytes, ec);
            fsnap = fopen(snap_path.string().c_str(), "ab");
        } else {
            fsnap = fopen(snap_path.string().c_str(), "wb");
            if (fsnap) {
                fprintf(fsnap, "# MenDrive fields_snap: E, w, Sx, Txx -- в момент t; Hy,Hz,M,B -- t+dt/2; j -- t-dt "
                               "(x: центр зазора, стенки x=-a, x=+a; x_cm = x*length_unit_cm) omega0=%.12g a=%.12g dt=%.12g\n",
                        omega, a_geom, dt);
                fprintf(fsnap, "period,phase_target,phase,t,i,x,x_cm,region"); head_cols(fsnap); fputc('\n', fsnap);
            }
        }
        if (!fsnap) throw std::runtime_error("не удалось открыть " + snap_path.string());
    }

    void sample(const MenDriveSim& sim) {
        const int N = sim.N;
        sim.magnetization_B(Mx, My, Mz);
        auto toA = [N](const std::vector<long double>& b, int i) -> long double {
            if (i <= 0) return b[0];
            if (i >= N) return b[N - 1];
            return 0.5L * (b[i - 1] + b[i]);
        };
        for (int i = 0; i <= N; ++i) {
            const long double Ey = sim.Ey[i], Ez = sim.Ez[i];
            const long double Hy = toA(sim.Hy, i), Hz = toA(sim.Hz, i);
            // H^{n-1/2} лежит в буферах *_new после swap (неактивная пара -- нули)
            const long double Hyc = 0.5L * (Hy + toA(sim.Hy_new, i)), Hzc = 0.5L * (Hz + toA(sim.Hz_new, i));
            const long double mx = toA(Mx, i), my = toA(My, i), mz = toA(Mz, i);
            q[FQ_EY][i] = (double)Ey; q[FQ_EZ][i] = (double)Ez; q[FQ_HY][i] = (double)Hy; q[FQ_HZ][i] = (double)Hz;
            q[FQ_MX][i] = (double)mx; q[FQ_MY][i] = (double)my; q[FQ_MZ][i] = (double)mz;
            q[FQ_BY][i] = (double)(Hy + 4 * M_PI * my); q[FQ_BZ][i] = (double)(Hz + 4 * M_PI * mz);
            q[FQ_JEY][i] = (double)sim.j_e_y_A[i]; q[FQ_JEZ][i] = (double)sim.j_e_z_A[i];
            q[FQ_JMY][i] = (double)toA(sim.j_m_y_B, i); q[FQ_JMZ][i] = (double)toA(sim.j_m_z_B, i);
            q[FQ_W][i] = (double)((sim.eps_A[i] * (Ey * Ey + Ez * Ez) + Hyc * (Hyc + 4 * M_PI * my) + Hzc * (Hzc + 4 * M_PI * mz)) / (8 * M_PI));
            q[FQ_SX][i] = (double)((Ey * Hzc - Ez * Hyc) / (4 * M_PI));
            q[FQ_TXX][i] = (double)((Ey * Ey + Ez * Ez + Hyc * Hyc + Hzc * Hzc) / (8 * M_PI));
        }
    }

    static int region_of(const MenDriveSim& sim, int i) { return sim.leftA[i] ? -1 : (sim.rightA[i] ? 1 : 0); }

    // вызывается после КАЖДОГО шага n (sim.t = (n+1) dt)
    void on_step(int64_t n, const MenDriveSim& sim) {
        if (!snap_on && !harm_on) return;
        const double t = (double)sim.t, t_prev = t - dt;
        bool sampled = false;
        if (snap_on) {
            for (double ph : phases) {
                const double pf = std::floor(t / T - ph);
                if (pf < 0) continue;
                const long long p = (long long)pf;
                if (!((pf + ph) * T > t_prev)) continue;   // момент (p+ph)T попал в (t-dt, t]
                if (p < snap_from) continue;
                if (snap_every > 0) { if ((p - snap_from) % snap_every != 0) continue; }
                else if (p != p_last) continue;
                if (!sampled) { sample(sim); sampled = true; }
                const double phase = t / T - (double)p;
                for (int i = 0; i <= sim.N; i += every_x) write_snap_row(sim, p, ph, phase, t, i);
                if ((sim.N % every_x) != 0) write_snap_row(sim, p, ph, phase, t, sim.N);
            }
        }
        if (harm_on && harm_valid && n >= harm_start) {
            if (!sampled) { sample(sim); sampled = true; }
            for (size_t h = 0; h < harm.size(); ++h) {
                const double th = std::fmod(sgn * harm[h] * omega * t, 2 * M_PI);
                const double zr = std::cos(th), zi = std::sin(th);
                for (int k = 0; k < FQ_COUNT; ++k) {
                    std::vector<double>& ar = acc_re[h * FQ_COUNT + k];
                    std::vector<double>& ai = acc_im[h * FQ_COUNT + k];
                    const std::vector<double>& v = q[k];
                    for (int i = 0; i <= sim.N; ++i) { ar[i] += v[i] * zr; ai[i] += v[i] * zi; }
                }
            }
            harm_count++;
        }
    }

    void write_snap_row(const MenDriveSim& sim, long long p, double ph, double phase, double t, int i) {
        fprintf(fsnap, "%lld,%.6g,%.10g,%.12e,%d,%.12e,%.12e,%d", p, ph, phase, t, i, (double)sim.xA[i],
                (double)sim.xA[i] * x_unit, region_of(sim, i));
        for (int k = 0; k < FQ_COUNT; ++k) fprintf(fsnap, ",%.10e", q[k][i]);
        fputc('\n', fsnap);
    }

    // ---- чекпойнт экспорта (отдельный файл: формат checkpoint.bin не менялся) ----
    void save() {
        if (!snap_on && !harm_on) return;
        if (fsnap) { fflush(fsnap); snap_bytes = (uint64_t)ftell(fsnap); }
        const fs::path tmp = ck_path.string() + ".tmp";
        FILE* f = fopen(tmp.string().c_str(), "wb");
        if (!f) return;
        const char mg[8] = {'M', 'D', 'F', 'X', '0', '0', '0', '1'};
        uint64_t ns = fsig.size();
        int32_t hv = harm_valid ? 1 : 0;
        bool ok = fwrite(mg, 1, 8, f) == 8 && ck_write_pod(f, ns) && fwrite(fsig.data(), 1, ns, f) == ns
                  && ck_write_pod(f, harm_count) && ck_write_pod(f, snap_bytes) && ck_write_pod(f, hv);
        for (size_t k = 0; ok && k < acc_re.size(); ++k) ok = ck_write_vec(f, acc_re[k]) && ck_write_vec(f, acc_im[k]);
        ok = (fclose(f) == 0) && ok;
        if (ok) { std::error_code ec; fs::rename(tmp, ck_path, ec); }
    }
    // resume: true, если состояние экспорта восстановлено
    bool load(int64_t n_done) {
        if (!snap_on && !harm_on) return true;
        FILE* f = fopen(ck_path.string().c_str(), "rb");
        bool ok = false;
        if (f) {
            char mg[8]; uint64_t ns = 0; int32_t hv = 1;
            if (fread(mg, 1, 8, f) == 8 && std::memcmp(mg, "MDFX0001", 8) == 0 && ck_read_pod(f, ns) && ns < (1u << 20)) {
                std::string sg(ns, '\0');
                if (fread(&sg[0], 1, ns, f) == ns && sg == fsig && ck_read_pod(f, harm_count) && ck_read_pod(f, snap_bytes)
                    && ck_read_pod(f, hv)) {
                    ok = true;
                    for (size_t k = 0; ok && k < acc_re.size(); ++k) ok = ck_read_vec(f, acc_re[k]) && ck_read_vec(f, acc_im[k]);
                    harm_valid = (hv != 0);
                }
            }
            fclose(f);
        }
        if (!ok) {
            printf("[fields] ВНИМАНИЕ: нет совместимого checkpoint_fields.bin (другие параметры экспорта или старый чекпойнт)\n");
            for (auto& v : acc_re) std::fill(v.begin(), v.end(), 0.0);
            for (auto& v : acc_im) std::fill(v.begin(), v.end(), 0.0);
            harm_count = 0; snap_bytes = 0;
            if (harm_on && n_done > harm_start) {
                harm_valid = false;
                printf("[fields] окно гармоник уже началось -> гармоники в этом прогоне НЕ будут записаны\n");
            }
            open_snap(false);
            return false;
        }
        open_snap(true);
        return true;
    }

    void finish(const MenDriveSim& sim, const fs::path& out_dir) {
        if (fsnap) { fclose(fsnap); fsnap = nullptr; printf("[fields] снимки: %s\n", snap_path.string().c_str()); }
        if (!harm_on) return;
        if (!harm_valid || harm_count == 0) { printf("[fields] гармоники не записаны (нет отсчётов)\n"); return; }
        const double T_rec = (double)harm_count * dt;
        for (size_t h = 0; h < harm.size(); ++h) {
            const int n = harm[h];
            const fs::path path = out_dir / ("fields_harm_n" + std::to_string(n) + ".csv");
            FILE* f = fopen(path.string().c_str(), "w");
            if (!f) continue;
            const double cn = (n == 0 ? 1.0 : 2.0) * dt / T_rec;
            // нормировка: j0 = i*s*amp -- комплексная амплитуда источника (sin(w t) = Re[i e^{-i w t}])
            std::complex<double> j0(0.0, sgn * amp);
            const bool do_norm = norm_src && n >= 1 && amp != 0;
            fprintf(f, "# MenDrive fields_harm: n=%d omega_n=%.12g omega0=%.12g T_rec=%.12g periods=%.3f samples=%lld "
                       "a=%.12g h_l=%.12g h_r=%.12g amp=%.12g\n",
                    n, n * omega, omega, T_rec, T_rec / T, (long long)harm_count,
                    (double)sim.a, (double)sim.h_l, (double)sim.h_r, amp);
            fprintf(f, "# time_convention=%s: F(x,t) = Re[F~(x) exp(%si n omega0 t)] (n=0: F~ = среднее)\n",
                    sgn > 0 ? "minus" : "plus", sgn > 0 ? "-" : "+");
            fprintf(f, "# normalize=%s%s; excitation_mode=%s; polarization=%s(%s); x: центр зазора, стенки x=-a,+a; "
                       "x_cm=x*%.6g\n", do_norm ? "source (F~/j0, j0=" : "none", do_norm ? (sgn > 0 ? "+i*amp)" : "-i*amp)") : "",
                    exc_name.c_str(), pol_name.c_str(), sim.upd_TE && sim.upd_TM ? "both" : (sim.upd_TE ? "TE" : "TM"), x_unit);
            fprintf(f, "n,omega_n,i,x,x_cm,region");
            for (int k = 0; k < FQ_COUNT; ++k) fprintf(f, ",%s_re,%s_im,%s_abs,%s_arg", FQ_NAMES[k], FQ_NAMES[k], FQ_NAMES[k], FQ_NAMES[k]);
            fputc('\n', f);
            // фазовые поправки на истинное время отсчёта
            std::complex<double> corr[3];
            corr[0] = 1.0;
            corr[1] = std::polar(1.0, sgn * n * omega * 0.5 * dt);
            corr[2] = std::polar(1.0, -sgn * n * omega * dt);
            for (int i = 0; i <= sim.N; ++i) {
                if (i % every_x != 0 && i != sim.N) continue;
                fprintf(f, "%d,%.12g,%d,%.12e,%.12e,%d", n, n * omega, i, (double)sim.xA[i], (double)sim.xA[i] * x_unit,
                        region_of(sim, i));
                for (int k = 0; k < FQ_COUNT; ++k) {
                    std::complex<double> v(acc_re[h * FQ_COUNT + k][i] * cn, acc_im[h * FQ_COUNT + k][i] * cn);
                    v *= corr[FQ_TCLASS[k]];
                    if (do_norm) v /= j0;
                    fprintf(f, ",%.12e,%.12e,%.12e,%.12e", v.real(), v.imag(), std::abs(v), std::arg(v));
                }
                fputc('\n', f);
            }
            fclose(f);
            printf("[fields] гармоника n=%d: %s\n", n, path.string().c_str());
        }
    }
};

static void print_usage() {
    printf("Использование: mendrive [-c файл.cfg] [ключ=значение ...] [--resume]\n\n");
    printf("  -c, --config F   прочитать параметры из файла (строки ключ=значение, # -- комментарий)\n");
    printf("  --resume         продолжить прогон из <out_dir>/checkpoint.bin (или недосчитанный скан)\n");
    printf("  -h, --help       эта справка\n\n");
    printf("Параметры применяются в порядке появления в командной строке.\n\n");
    printf("%-26s %-16s %s\n", "ключ", "по умолчанию", "описание");
    for (const auto& p : PARAMS) printf("%-26s %-16s %s\n", p.key, p.def, p.help);
}

// ============================================================================
// main
// ============================================================================
static int real_main(int argc, char** argv) {
    setvbuf(stdout, nullptr, _IOLBF, 0);   // построчный вывод (удобно с nohup/tee)
    Config cfg;
    bool resume = false;
    for (int k = 1; k < argc; ++k) {
        std::string a = argv[k];
        if (a == "-h" || a == "--help") { print_usage(); return 0; }
        else if (a == "--resume") resume = true;
        else if (a == "-c" || a == "--config") {
            if (k + 1 >= argc) throw std::runtime_error("после " + a + " нужен путь к файлу");
            cfg.load_file(argv[++k]);
        } else if (a.rfind("--config=", 0) == 0) cfg.load_file(a.substr(9));
        else if (a.find('=') != std::string::npos) {
            size_t eq = a.find('=');
            cfg.set(trim(a.substr(0, eq)), trim(a.substr(eq + 1)));
        } else throw std::runtime_error("непонятный аргумент '" + a + "' (см. --help)");
    }
    // ранняя проверка всех строковых кодов
    const int fm = ferrite_code(cfg); excitation_code(cfg); bias_code(cfg); pr2d_algo_code(cfg); polarization_code(cfg);
    if (cfg.s("ts_format") != "csv" && cfg.s("ts_format") != "bin") throw std::runtime_error("ts_format: csv | bin");
    const int threads = (int)std::max(1LL, cfg.i("threads"));

    const fs::path out_dir = cfg.s("out_dir");
    fs::create_directories(out_dir);
    { FILE* f = fopen((out_dir / "config_used.txt").string().c_str(), "w"); if (f) { cfg.dump(f); fclose(f); } }
    std::signal(SIGINT, on_signal);
    std::signal(SIGTERM, on_signal);

    const auto t_start = std::chrono::steady_clock::now();
    auto elapsed = [&]() { return std::chrono::duration<double>(std::chrono::steady_clock::now() - t_start).count(); };

    // ---------- главная модель и оценка памяти ----------
    auto sim = make_sim(cfg, threads);
    const double dt = (double)sim->dt;
    const int probe_idx = (int)cfg.i("probe_idx");
    if (probe_idx < 0 || probe_idx >= sim->n_fer)
        throw std::runtime_error("probe_idx=" + std::to_string(probe_idx) + " вне [0, n_fer=" + std::to_string(sim->n_fer) + ")");
    const int gi = sim->probe_global_index(probe_idx);
    printf("[init] модель=%s, возбуждение=%s, N=%d, dx=%.6g, dt=%.6g, ячеек феррита n_fer=%d, потоков=%d\n",
           cfg.s("ferrite_model").c_str(), cfg.s("excitation_mode").c_str(), sim->N, (double)sim->dx, dt, sim->n_fer, threads);
    printf("[mode] %s\n", sim->mode_info().c_str());   // MODE
    if (fm == Preisach2D) {
        const double ng = sim->pr2d->n_grid, pairs = ng * ng;
        printf("[init] Preisach2D: n_grid=%d, пар=%.0f, алгоритм=%s\n", sim->pr2d->n_grid, pairs, cfg.s("pr2d_algo").c_str());
        printf("[init]   память состояния: сейчас %s; в исходном коде было бы %s (int8) + %s (alpha/beta/weight)\n",
               fmt_bytes(sim->state_bytes()).c_str(), fmt_bytes(sim->n_fer * pairs).c_str(),
               fmt_bytes(3.0 * pairs * sizeof(long double)).c_str());
        if (cfg.i("n_grid") <= 0 && ng > 300)
            printf("[init]   ВНИМАНИЕ: n_grid взят из n_hyst_pr=%d (как в оригинале) -> n_grid^2 пар. Задайте n_grid явно.\n", (int)ng);
    } else {
        printf("[init] память состояния феррита: %s\n", fmt_bytes(sim->state_bytes()).c_str());
    }

    const fs::path ck_path = out_dir / "checkpoint.bin";
    RunState R;
    double omega_res = NaN;
    FILE* ckf = nullptr;
    if (resume && fs::exists(ck_path)) {
        ckf = open_checkpoint(ck_path, cfg, omega_res);
        printf("[resume] найден чекпойнт, omega0=%.10g\n", omega_res);
    } else if (resume) {
        printf("[resume] чекпойнта нет -- начинаю сначала (готовые точки скана будут переиспользованы)\n");
    }

    // ---------- скан резонанса ----------
    if (!ckf) {
        if (cfg.d("omega0") > 0) {
            omega_res = cfg.d("omega0");
            printf("[scan] пропущен: omega0=%.10g задана явно\n", omega_res);
        } else {
            std::vector<double> freqs_wide = arange(cfg.d("scan_lo"), cfg.d("scan_hi"), cfg.d("scan_step"));
            if (freqs_wide.empty()) throw std::runtime_error("пустая сетка грубого скана");
            std::vector<double> resp_wide;
            if (!run_scan(cfg, "wide", freqs_wide, threads, out_dir, resume, resp_wide)) {
                printf("[scan] прервано; продолжить: тот же запуск с --resume\n"); return 130;
            }
            long long iw = nanargmax(resp_wide);
            if (iw < 0) throw std::runtime_error("грубый скан: все отклики NaN (модель «взрывается»)");
            const double omega_coarse = freqs_wide[(size_t)iw];
            printf("[main] Грубый скан: пик omega0~%.4g, %s\n", omega_coarse, fmt_dur(elapsed()).c_str());
            const double lo = std::max(freqs_wide[0], omega_coarse - cfg.d("fine_halfwidth"));
            const double hi = omega_coarse + cfg.d("fine_halfwidth");
            std::vector<double> freqs_fine = arange(lo, hi + 1e-9, cfg.d("fine_step"));
            std::vector<double> resp_fine;
            if (!run_scan(cfg, "fine", freqs_fine, threads, out_dir, resume, resp_fine)) {
                printf("[scan] прервано; продолжить: тот же запуск с --resume\n"); return 130;
            }
            long long ifn = nanargmax(resp_fine);
            omega_res = (ifn >= 0) ? freqs_fine[(size_t)ifn] : omega_coarse;
            printf("[main] Уточнённый резонанс: omega0=%.6g, %s\n", omega_res, fmt_dur(elapsed()).c_str());
            if (ifn == 0 || ifn == (long long)freqs_fine.size() - 1)
                printf("[main] ВНИМАНИЕ: максимум на краю окна уточняющего скана -- резонанс может быть вне окна!\n");
            if (iw == 0 || iw == (long long)freqs_wide.size() - 1)
                printf("[main] ВНИМАНИЕ: максимум на краю грубого скана -- расширьте scan_lo/scan_hi.\n");
            FILE* f = fopen((out_dir / "scan_result.txt").string().c_str(), "w");
            if (f) { fprintf(f, "omega_coarse=%.17g\nomega_res=%.17g\n", omega_coarse, omega_res); fclose(f); }
        }
    }

    // ---------- план длинного прогона ----------
    plan_run(cfg, omega_res, dt, R);
    const std::string sig = cfg.signature();
    const bool binary_ts = cfg.s("ts_format") == "bin";
    const long long ts_every = cfg.i("ts_every");
    const fs::path per_path = out_dir / "periods.csv";
    const fs::path ts_path = out_dir / (binary_ts ? "timeseries.bin" : "timeseries.csv");
    FILE* fper = nullptr; FILE* fts = nullptr;
    if (ckf) {
        const RunState planned = R;
        bool ok = load_run(ckf, R) && sim->load_state(ckf);
        fclose(ckf); ckf = nullptr;
        if (!ok || R.nsteps != planned.nsteps || R.spp != planned.spp || R.n_use != planned.n_use)
            throw std::runtime_error("чекпойнт не соответствует конфигурации/версии программы");
        std::error_code ec;
        if (fs::exists(per_path)) fs::resize_file(per_path, R.per_bytes, ec);
        if (ts_every > 0 && fs::exists(ts_path)) fs::resize_file(ts_path, R.ts_bytes, ec);
        fper = fopen(per_path.string().c_str(), "ab");
        if (ts_every > 0) fts = fopen(ts_path.string().c_str(), "ab");
        printf("[resume] продолжаю с шага %" PRId64 "/%" PRId64 " (период %.1f)\n", R.n_done, R.nsteps,
               (double)R.n_done / (double)R.spp);
    } else {
        fper = fopen(per_path.string().c_str(), "wb");
        if (fper) fprintf(fper, "block,t_end,dTxx,P,P_src,Pdiss,P_num,Hn_min,Hn_max,W_field,cum_dTxx_over_P,cum_dTxx_over_Psrc\n");
        if (ts_every > 0) {
            fts = fopen(ts_path.string().c_str(), "wb");
            const char* cols = "t,dTxx,Hn,P,P_src,Pdiss,Pe,Pm_wall,Pmag,P_num,Mx,My,Mz,HzJA,MzJA,W_field";
            if (fts && !binary_ts) fprintf(fts, "%s\n", cols);
            if (binary_ts) { FILE* f = fopen((out_dir / "ts_columns.txt").string().c_str(), "w"); if (f) { fprintf(f, "%s\n", cols); fclose(f); } }
        }
    }
    if (!fper || (ts_every > 0 && !fts)) throw std::runtime_error("не удалось открыть выходные файлы в " + out_dir.string());
    // FIELDS: экспорт профилей (снимки и гармоники); состояние -- в checkpoint_fields.bin
    FieldsExport fx;
    fx.init(cfg, *sim, R, out_dir);
    if (R.n_done > 0) fx.load(R.n_done); else fx.open_snap(false);

    {
        const double ts_rows = (ts_every > 0) ? (double)(R.nsteps - R.ts_from) / (double)ts_every : 0.0;
        const double mem_now = sim->state_bytes() + 16.0 * (sim->N + 1) * sizeof(long double)
                               + 2.0 * R.spp * 8 + 5.0 * R.fft_w.size() * 8;
        const double mem_nb = (double)R.n_rec * 8.0 * 13 + (double)R.n_rec * 8.0 * 4;
        printf("[run] omega0=%.10g, T=%.6g, шагов=%" PRId64 ", шагов/период=%" PRId64 ", запись с шага %" PRId64
               " (%" PRId64 " точек), окно средних=%" PRId64 " шагов (%.0f периодов)\n",
               R.omega, R.T, R.nsteps, R.spp, R.record_start, R.n_rec, R.n_use, (double)R.n_use / R.spp);
        printf("[run] память: ~%s в этой программе; ноутбук держал бы ~%s буферов записи+FFT\n",
               fmt_bytes(mem_now).c_str(), fmt_bytes(mem_nb).c_str());
        printf("[run] временной ряд: ~%.0f строк -> %s\n", ts_rows, ts_path.string().c_str());
    }

    // ---------- длинный прогон ----------
    const double amp = cfg.d("amp"), ramp = cfg.d("ramp_periods");
    long long pe = cfg.i("progress_every_periods");
    if (pe <= 0) pe = std::max(1LL, (long long)(cfg.d("n_periods_total") / 100.0));
    const int64_t progress_steps = std::max<int64_t>(1, (int64_t)pe * R.spp);
    const int64_t ck_steps = std::max<int64_t>(1, std::max(1LL, cfg.i("checkpoint_every_periods")) * R.spp);
    const double ck_sec = cfg.d("checkpoint_every_sec");
    const bool loop_from_JA = (fm == Hybrid);
    const int64_t n_session0 = R.n_done;
    double last_ck = elapsed();
    const double t_loop0 = last_ck;
    double last_blk[5] = {NaN, NaN, NaN, NaN, NaN};   // dTxx, P_src, Pdiss, P_num, max|Hn|
    const double dt_d = R.dt;
    // поворот фазора ДПФ за один шаг: exp(-i w dt)
    std::vector<double> rot_r(R.fft_w.size()), rot_i(R.fft_w.size());
    for (size_t m = 0; m < R.fft_w.size(); ++m) {
        rot_r[m] = std::cos(R.fft_w[m] * dt_d); rot_i[m] = -std::sin(R.fft_w[m] * dt_d);
    }

    for (int64_t n = R.n_done; n < R.nsteps; ++n) {
        const double P_step = (double)sim->step_once(R.omega, amp, ramp, R.T);
        if (n >= R.record_start) {
            const int64_t r = n - R.record_start;
            const double P = (double)(P_step / sim->dt);   // как out_P в mendrive_run (с ошибкой 1/dt)
            const double Psrc = P_step;                    // физическая мощность источника
            const double dTxx = (double)(sim->Txx_at(sim->iR_A) - sim->Txx_at(sim->iL_A));
            const double Pdiss = (double)sim->dissipated_power_density();
            long double pe_l, pm_l, pmag_l; sim->dissipated_power_parts(pe_l, pm_l, pmag_l);
            const double Pnum = (double)sim->p_num_residual;
            const ProbeVals pv = probe_values(*sim, probe_idx, gi);
            if (n >= R.win_start) {
                R.w_dTxx.add(dTxx); R.w_P.add(P); R.w_Psrc.add(Psrc); R.w_Pdiss.add(Pdiss);
                R.w_Pe.add((double)pe_l); R.w_Pm.add((double)pm_l); R.w_Pmag.add((double)pmag_l); R.w_Pnum.add(Pnum);
                R.w_count++;
                R.w_Hmin = std::min(R.w_Hmin, pv.Hn); R.w_Hmax = std::max(R.w_Hmax, pv.Hn);
            }
            // блок = период записи
            R.b_dTxx.add(dTxx); R.b_P.add(P); R.b_Psrc.add(Psrc); R.b_Pdiss.add(Pdiss); R.b_Pnum.add(Pnum);
            R.b_Hmin = std::min(R.b_Hmin, pv.Hn); R.b_Hmax = std::max(R.b_Hmax, pv.Hn);
            if (++R.b_count == R.spp) {
                const double sp = (double)R.spp;
                const double mF = R.b_dTxx.s / sp, mP = R.b_P.s / sp, mPs = R.b_Psrc.s / sp;
                R.cum_F.add(mF); R.cum_P.add(mP); R.cum_Psrc.add(mPs);
                const double rP = (std::fabs(R.cum_P.s) > 1e-30) ? R.cum_F.s / R.cum_P.s : NaN;
                const double rPs = (std::fabs(R.cum_Psrc.s) > 1e-300) ? R.cum_F.s / R.cum_Psrc.s : NaN;
                fprintf(fper, "%" PRId64 ",%.10e,%.10e,%.10e,%.10e,%.10e,%.10e,%.6e,%.6e,%.10e,%.10e,%.10e\n",
                        R.blocks_done, (double)sim->t, mF, mP, mPs, R.b_Pdiss.s / sp, R.b_Pnum.s / sp,
                        R.b_Hmin, R.b_Hmax, (double)sim->field_energy(), rP, rPs);
                last_blk[0] = mF; last_blk[1] = mPs; last_blk[2] = R.b_Pdiss.s / sp; last_blk[3] = R.b_Pnum.s / sp;
                last_blk[4] = std::max(std::fabs(R.b_Hmin), std::fabs(R.b_Hmax));
                R.blocks_done++;
                R.b_dTxx = Kahan(); R.b_P = Kahan(); R.b_Psrc = Kahan(); R.b_Pdiss = Kahan(); R.b_Pnum = Kahan();
                R.b_count = 0;
                R.b_Hmin = std::numeric_limits<double>::infinity(); R.b_Hmax = -std::numeric_limits<double>::infinity();
            }
            // последний период для петли: (Hn, My) или (HzJA, MzJA) для Hybrid -- как в ноутбуке
            R.ring_H[(size_t)R.ring_pos] = loop_from_JA ? pv.HzJA : pv.Hn;
            R.ring_M[(size_t)R.ring_pos] = loop_from_JA ? pv.MzJA : pv.My;
            R.ring_pos = (R.ring_pos + 1) % R.spp;
            R.ring_filled = std::min<int64_t>(R.spp, R.ring_filled + 1);
            // ДПФ: X(w) = sum_r hann(r) dTxx_r exp(-i w r dt) -- то же, что np.fft.rfft(dTxx*hanning) в узлах сетки
            if (!R.fft_w.empty()) {
                const double win = (R.n_rec > 1) ? 0.5 - 0.5 * std::cos(2 * M_PI * (double)r / (double)(R.n_rec - 1)) : 1.0;
                const double xw = dTxx * win;
                const size_t nb = R.fft_w.size();
                if ((r & 1023) == 0) {   // точный пересчёт фазора (без накопления ошибки вращения)
                    for (size_t m = 0; m < nb; ++m) {
                        const double th = std::fmod(R.fft_w[m] * dt_d * (double)r, 2 * M_PI);
                        R.fft_zr[m] = std::cos(th); R.fft_zi[m] = -std::sin(th);
                    }
                }
                for (size_t m = 0; m < nb; ++m) {
                    R.fft_re[m] += xw * R.fft_zr[m];
                    R.fft_im[m] += xw * R.fft_zi[m];
                    const double zr = R.fft_zr[m] * rot_r[m] - R.fft_zi[m] * rot_i[m];
                    const double zi = R.fft_zr[m] * rot_i[m] + R.fft_zi[m] * rot_r[m];
                    R.fft_zr[m] = zr; R.fft_zi[m] = zi;
                }
            }
            // временной ряд (прореженный, только хвост)
            if (fts && n >= R.ts_from && ((n - R.ts_from) % ts_every) == 0) {
                const double row[16] = {(double)sim->t, dTxx, pv.Hn, P, Psrc, Pdiss, (double)pe_l, (double)pm_l,
                                        (double)pmag_l, Pnum, pv.Mx, pv.My, pv.Mz, pv.HzJA, pv.MzJA,
                                        (double)sim->field_energy()};
                if (binary_ts) fwrite(row, sizeof(double), 16, fts);
                else {
                    for (int q = 0; q < 16; ++q) fprintf(fts, q ? ",%.10e" : "%.10e", row[q]);
                    fputc('\n', fts);
                }
            }
            R.rec_count++;
        }
        fx.on_step(n, *sim);   // FIELDS
        const bool finite = sim->fields_finite();
        R.n_done = n + 1;
        if (!finite) { R.blew_up = 1; printf("[run] ВНИМАНИЕ: поля стали нечисловыми (blew_up) на шаге %" PRId64 "\n", n); break; }

        if (R.n_done % progress_steps == 0) {
            const double el = elapsed();
            const double rate = (double)(R.n_done - n_session0) / std::max(1e-9, el - t_loop0);
            const double eta = (double)(R.nsteps - R.n_done) / std::max(1e-9, rate);
            printf("[run] период %.0f/%.0f (%.1f%%) | прошло %s | осталось ~%s | %.0f шаг/с | "
                   "посл. период: <dTxx>=%.4e <P_src>=%.4e <Pdiss>=%.4e <P_num>=%.4e max|Hn|=%.3g | RSS=%.0f МБ\n",
                   (double)R.n_done / R.spp, (double)R.nsteps / R.spp, 100.0 * R.n_done / R.nsteps,
                   fmt_dur(R.wall_prev + (el - last_ck)).c_str(), fmt_dur(eta).c_str(), rate,
                   last_blk[0], last_blk[1], last_blk[2], last_blk[3], last_blk[4], peak_rss_mb());
        }
        const bool boundary = (R.n_done % R.spp) == 0;
        if (g_stop || (boundary && (R.n_done % ck_steps == 0 || elapsed() - last_ck >= ck_sec))) {
            const double el = elapsed();
            R.wall_prev += el - last_ck; last_ck = el;
            write_checkpoint(ck_path, sig, R, *sim, fper, fts);
            fx.save();   // FIELDS
            if (g_stop) {
                printf("[run] остановлено на шаге %" PRId64 ", чекпойнт записан. Продолжить: тот же запуск + --resume\n", R.n_done);
                fclose(fper); if (fts) fclose(fts);
                if (fx.fsnap) { fclose(fx.fsnap); fx.fsnap = nullptr; }
                return 130;
            }
        }
    }
    { const double el = elapsed(); R.wall_prev += el - last_ck; last_ck = el; }
    write_checkpoint(ck_path, sig, R, *sim, fper, fts);
    fx.save();   // FIELDS
    fclose(fper); if (fts) fclose(fts);
    fx.finish(*sim, out_dir);   // FIELDS: fields_harm_n*.csv

    // ---------- итоги (порт force_per_power / check_energy_balance / shoelace_area) ----------
    const double wc = (double)R.w_count;
    const double dTxx_avg = wc > 0 ? R.w_dTxx.s / wc : NaN;
    const double P_avg = wc > 0 ? R.w_P.s / wc : NaN;
    const double Psrc_avg = wc > 0 ? R.w_Psrc.s / wc : NaN;
    const double Pdiss_avg = wc > 0 ? R.w_Pdiss.s / wc : NaN;
    const double Pe_avg = wc > 0 ? R.w_Pe.s / wc : NaN, Pm_avg = wc > 0 ? R.w_Pm.s / wc : NaN;
    const double Pmag_avg = wc > 0 ? R.w_Pmag.s / wc : NaN, Pnum_avg = wc > 0 ? R.w_Pnum.s / wc : NaN;
    auto fpp = [&](double Pv, double& code, double& per_kW) {
        code = (std::fabs(Pv) > 1e-300) ? dTxx_avg / Pv : NaN;
        const double force_N = dTxx_avg * cfg.d("area_cm2") * 1e-5;
        const double power_kW = Pv * std::pow(cfg.d("length_unit_cm"), 2) / cfg.d("time_unit_s") * 1e-7 * 1e-3;
        per_kW = (std::fabs(power_kW) > 1e-300) ? force_N / power_kW : NaN;
    };
    double fpp_code, fpkw, fpp_code_src, fpkw_src;
    fpp(P_avg, fpp_code, fpkw);
    fpp(Psrc_avg, fpp_code_src, fpkw_src);
    auto rel = [](double x, double y) { return std::fabs(x - y) / std::max({std::fabs(x), std::fabs(y), 1e-300}); };
    const double rel_dif = rel(P_avg, Pdiss_avg);
    const double rel_src = rel(Psrc_avg, Pdiss_avg);
    const double rel_num = rel(Psrc_avg + Pnum_avg, Pdiss_avg);

    // петля гистерезиса за последний период
    double hyst_area = NaN;
    const int64_t nl = R.ring_filled;
    const bool loop_ok = (fm == JA || fm == Preisach || fm == Preisach2D) ? nl > 0 : (fm == Hybrid ? nl >= R.spp : false);
    {
        FILE* f = fopen((out_dir / "loop.csv").string().c_str(), "w");
        if (f) fprintf(f, "H,M,B\n");
        if (nl > 0) {
            const int64_t start = (R.ring_pos - nl + R.spp) % R.spp;
            double sacc = 0.0;
            for (int64_t q = 0; q < nl; ++q) {
                const size_t k0 = (size_t)((start + q) % R.spp), k1 = (size_t)((start + (q + 1) % nl) % R.spp);
                const double x0 = R.ring_H[k0], y0 = x0 + 4 * M_PI * R.ring_M[k0];
                const double x1 = R.ring_H[k1], y1 = x1 + 4 * M_PI * R.ring_M[k1];
                sacc += x0 * y1 - x1 * y0;
                if (f) fprintf(f, "%.10e,%.10e,%.10e\n", x0, R.ring_M[k0], y0);
            }
            if (loop_ok) hyst_area = 0.5 * std::fabs(sacc);
        }
        if (f) fclose(f);
    }
    // спектр
    {
        const long long K = cfg.i("fft_harmonics");
        FILE* fh = fopen((out_dir / "harmonics.csv").string().c_str(), "w");
        FILE* fs_ = fopen((out_dir / "spectrum.csv").string().c_str(), "w");
        if (fh) fprintf(fh, "k,omega,abs,phase\n");
        if (fs_) fprintf(fs_, "omega,abs\n");
        for (size_t m = 0; m < R.fft_w.size(); ++m) {
            const double ab = std::hypot(R.fft_re[m], R.fft_im[m]);
            if ((long long)m <= K) { if (fh) fprintf(fh, "%zu,%.10e,%.10e,%.10e\n", m, R.fft_w[m], ab, std::atan2(R.fft_im[m], R.fft_re[m])); }
            else if (fs_) fprintf(fs_, "%.10e,%.10e\n", R.fft_w[m], ab);
        }
        if (fh) fclose(fh);
        if (fs_) fclose(fs_);
    }

    const double total_wall = R.wall_prev;
    FILE* fsum = fopen((out_dir / "summary.txt").string().c_str(), "w");
    auto both = [&](const char* fmt, auto... args) {
        printf(fmt, args...);
        if (fsum) fprintf(fsum, fmt, args...);
    };
    both("# Итоги прогона MenDrive (main.cpp)\n");
    both("ferrite_model=%s\nexcitation_mode=%s\namp=%.10g\nN=%d\nn_fer=%d\n", cfg.s("ferrite_model").c_str(),
         cfg.s("excitation_mode").c_str(), amp, sim->N, sim->n_fer);
    both("polarization=%s\n", sim->upd_TE && sim->upd_TM ? "both" : (sim->upd_TE ? "TE" : "TM"));   // MODE
    both("omega_res=%.17g\ndt=%.17g\nspp=%" PRId64 "\nn_rec=%" PRId64 "\nn_window=%" PRId64 "\nblew_up=%d\n",
         R.omega, dt_d, R.spp, R.rec_count, R.w_count, (int)R.blew_up);
    both("# средние по окну (последние last_frac записанных периодов)\n");
    both("dTxx_avg=%.10e\nP_avg=%.10e\nP_src_avg=%.10e\nPdiss_avg=%.10e\n", dTxx_avg, P_avg, Psrc_avg, Pdiss_avg);
    both("Pdiss_sigma_e_E2=%.10e\nPdiss_sigma_m_H2=%.10e\nPdiss_H_dMdt=%.10e\nP_num_residual=%.10e\n", Pe_avg, Pm_avg, Pmag_avg, Pnum_avg);
    both("Hn_min_window=%.6e\nHn_max_window=%.6e\n", R.w_Hmin, R.w_Hmax);
    both("# как в ноутбуке (P = P_src/dt -- завышено в 1/dt раз)\n");
    both("force_per_power_code=%.10e\nforce_per_kW=%.10e\nrel_dif=%.10e\n", fpp_code, fpkw, rel_dif);
    both("# исправленная нормировка мощности (P_src = -<j.E>)\n");
    both("force_per_power_code_src=%.10e\nforce_per_kW_src=%.10e\nrel_dif_src=%.10e\nrel_dif_src_plus_num=%.10e\n",
         fpp_code_src, fpkw_src, rel_src, rel_num);
    both("hysteresis_area=%.10e\n", hyst_area);
    both("wall_time_s=%.1f\npeak_rss_mb=%.1f\n", total_wall, peak_rss_mb());
    if (fsum) fclose(fsum);
    printf("[main] FORCE/POWER: %.4e (код.ед.), %.4e Н/кВт   [с P_src: %.4e, %.4e Н/кВт]\n", fpp_code, fpkw, fpp_code_src, fpkw_src);
    printf("[main] check_energy_balance: P_avg = %.4e, Pdiss_avg=%.4e rel_dif=%.6g\n", P_avg, Pdiss_avg, rel_dif);
    printf("[main] баланс с P_src: P_src=%.4e, P_num=%.4e, Pdiss=%.4e -> rel(P_src,Pdiss)=%.3g, rel(P_src+P_num,Pdiss)=%.3g\n",
           Psrc_avg, Pnum_avg, Pdiss_avg, rel_src, rel_num);
    printf("[main] Итого: %s, пик памяти %.0f МБ, результаты в %s\n", fmt_dur(total_wall).c_str(), peak_rss_mb(), out_dir.string().c_str());
    return 0;
}

int main(int argc, char** argv) {
    try {
        return real_main(argc, argv);
    } catch (const std::bad_alloc&) {
        fprintf(stderr, "ОШИБКА: не хватает памяти (std::bad_alloc). Уменьшите n_grid/N или используйте pr2d_algo=rows.\n");
        return 3;
    } catch (const std::exception& e) {
        fprintf(stderr, "ОШИБКА: %s\n", e.what());
        return 2;
    }
}
