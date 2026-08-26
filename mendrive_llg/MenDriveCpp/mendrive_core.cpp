// mendrive_core.cpp
// C++ порт резонатора МенДрайв (JA / LLG / Hybrid), СГС.
// extern "C" API для загрузки через ctypes (без pybind11/cmake).
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
#include <algorithm>
#include <memory>

static inline double langevin_L(double x) {
    if (std::fabs(x) < 1e-6) { return x/3.0 - x*x*x/45.0; }
    return 1.0/std::tanh(x) - 1.0/x;
}

struct JAEngine {
    double Ms, a, alpha, k, c, sigma_m_leak, H0_bias;
    int n_sub, n_fp;
    JAEngine(double Ms_, double a_, double alpha_, double k_, double c_,
             double sml, double H0b, int ns, int nf)
        : Ms(Ms_), a(a_), alpha(alpha_), k(k_), c(c_),
          sigma_m_leak(sml), H0_bias(H0b), n_sub(ns), n_fp(nf) {}

    inline void M_of_Mirr(double M_irr, double H, double M_guess,
                           double& M_out, double& Man_out) const {
        double M = M_guess;
        for (int i = 0; i < n_fp; ++i) {
            double He = (H + H0_bias) + alpha * M;
            double Man = Ms * langevin_L(He / a);
            M = (1.0 - c) * M_irr + c * Man;
        }
        double He = (H + H0_bias) + alpha * M;
        double Man = Ms * langevin_L(He / a);
        M_out = M; Man_out = Man;
    }

    inline void rhs(double M_irr, double H, double delta, double M_guess,
                     double& drhs, double& M_out) const {
        double M, Man;
        M_of_Mirr(M_irr, H, M_guess, M, Man);
        double denom = k * delta - alpha * (Man - M_irr);
        double denom_safe = (std::fabs(denom) < 1e-9)
                             ? 1e-9 * (denom >= 0 ? 1.0 : -1.0) : denom;
        drhs = (Man - M_irr) / denom_safe;
        M_out = M;
    }

    void integrate_rk4(double M_irr_old, double H_old, double H_new,
                        double M_guess_old, double delta,
                        double& M_irr_out, double& M_out) const {
        double h = (H_new - H_old) / n_sub;
        double M_irr = M_irr_old, M_guess = M_guess_old, H_cur = H_old;
        for (int s = 0; s < n_sub; ++s) {
            double k1, M1; rhs(M_irr, H_cur, delta, M_guess, k1, M1);
            double k2, M2; rhs(M_irr + 0.5*h*k1, H_cur + 0.5*h, delta, M1, k2, M2);
            double k3, M3; rhs(M_irr + 0.5*h*k2, H_cur + 0.5*h, delta, M2, k3, M3);
            double k4, M4; rhs(M_irr + h*k3, H_cur + h, delta, M3, k4, M4);
            M_irr = M_irr + (h/6.0)*(k1 + 2*k2 + 2*k3 + k4);
            M_irr = std::max(-Ms, std::min(Ms, M_irr));
            H_cur += h; M_guess = M4;
        }
        double M_final, Man_final;
        M_of_Mirr(M_irr, H_new, M_guess, M_final, Man_final);
        M_irr_out = M_irr;
        M_out = std::max(-Ms, std::min(Ms, M_final));
    }

    // Производная dM/dH считается ЧИСЛЕННО (конечная разность), а не по
    // аналитической формуле -- избегаем сингулярности 1/Heff^2 при Heff->0,
    // которая присутствует в наивной прямой реализации.
    void dM_dH(double M_irr_old, double H_old, double H_new,
               double M_guess_old, double delta,
               double& dMdH_out, double& M0_out) const {
        double eps = 1e-6 * std::max(1.0, std::fabs(H_new));
        double tmp1, M0; integrate_rk4(M_irr_old, H_old, H_new, M_guess_old, delta, tmp1, M0);
        double tmp2, M1; integrate_rk4(M_irr_old, H_old, H_new + eps, M_guess_old, delta, tmp2, M1);
        dMdH_out = (M1 - M0) / eps; M0_out = M0;
    }

    void newton_branch(double M_irr_old, double H_old, double M_old, double rhs_extra,
                        double dt, double delta, int n_newton, double max_step,
                        double& H_new_out, double& F_out) const {
        double H_new = H_old - rhs_extra;
        double leak_coef = 4*M_PI*sigma_m_leak*dt;
        double dMdH, Mc;
        for (int it = 0; it < n_newton; ++it) {
            dM_dH(M_irr_old, H_old, H_new, M_old, delta, dMdH, Mc);
            double F = H_new - H_old + 4*M_PI*(Mc - M_old) + leak_coef*H_new + rhs_extra;
            double Fp = 1.0 + 4*M_PI*dMdH + leak_coef;
            double Fp_safe = (std::fabs(Fp) < 0.5) ? 0.5*(Fp >= 0 ? 1.0 : -1.0) : Fp;
            double step = std::max(-max_step, std::min(max_step, F/Fp_safe));
            H_new -= step;
        }
        dM_dH(M_irr_old, H_old, H_new, M_old, delta, dMdH, Mc);
        F_out = H_new - H_old + 4*M_PI*(Mc - M_old) + leak_coef*H_new + rhs_extra;
        H_new_out = H_new;
    }

    // Двусторонняя проверка знака delta=sign(dH/dt), фиксируется единожды по
    // знаку внешнего толчка -- избегает паразитного цикла период-2.
    void implicit_step(double M_irr_old, double H_old, double M_old, double rhs_extra,
                        double dt, int n_newton, double max_step_factor, double tol,
                        double& H_new_out, double& M_irr_out, double& M_out) const {
        double max_step = max_step_factor * a;
        double delta_pred = (-rhs_extra >= 0) ? 1.0 : -1.0;
        double H_a, F_a;
        newton_branch(M_irr_old, H_old, M_old, rhs_extra, dt, delta_pred, n_newton, max_step, H_a, F_a);
        double H_new = H_a, delta_final = delta_pred;
        if (std::fabs(F_a) > tol) {
            double delta_alt = -delta_pred;
            double H_b, F_b;
            newton_branch(M_irr_old, H_old, M_old, rhs_extra, dt, delta_alt, n_newton, max_step, H_b, F_b);
            if (std::fabs(F_b) < std::fabs(F_a)) { H_new = H_b; delta_final = delta_alt; }
        }
        double M_irr_f, M_f;
        integrate_rk4(M_irr_old, H_old, H_new, M_old, delta_final, M_irr_f, M_f);
        H_new_out = H_new; M_irr_out = M_irr_f; M_out = M_f;
    }
};

struct Vec3 { double x, y, z; };

struct LLGEngine {
    double Ms, gamma, alpha, H0_bias;
    int bias_axis; // 0='none', 1='x', 2='y', 3='z'
    LLGEngine(double Ms_, double gamma_, double alpha_, double H0_bias_, int bias_axis_)
        : Ms(Ms_), gamma(gamma_), alpha(alpha_), H0_bias(H0_bias_), bias_axis(bias_axis_) {}

    inline Vec3 build_Heff(double Hy_rf, double Hz_rf) const {
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
    static inline double dot(const Vec3&a, const Vec3&b) { return a.x*b.x+a.y*b.y+a.z*b.z; }
    static inline double norm(const Vec3&a) { return std::sqrt(dot(a,a)); }

    Vec3 rodrigues_rotate(const Vec3& M, const Vec3& Heff, double dt) const {
        double Hn = norm(Heff);
        if (Hn <= 1e-300) return M;
        Vec3 axis{Heff.x/Hn, Heff.y/Hn, Heff.z/Hn};
        double theta = gamma * Hn * dt;
        double ct = std::cos(theta), st = std::sin(theta);
        double d = dot(axis, M);
        Vec3 crs = cross(axis, M);
        return Vec3{ M.x*ct + crs.x*st + axis.x*d*(1-ct),
                     M.y*ct + crs.y*st + axis.y*d*(1-ct),
                     M.z*ct + crs.z*st + axis.z*d*(1-ct) };
    }
    Vec3 llg_step_core(const Vec3& M, const Vec3& Heff, double dt) const {
        Vec3 M1 = rodrigues_rotate(M, Heff, dt);
        Vec3 cross1 = cross(M1, Heff);
        Vec3 damp = cross(M1, cross1);
        Vec3 M2{ M1.x - dt*(gamma*alpha/Ms)*damp.x,
                 M1.y - dt*(gamma*alpha/Ms)*damp.y,
                 M1.z - dt*(gamma*alpha/Ms)*damp.z };
        double n2 = norm(M2);
        double n2s = (n2 > 1e-300) ? n2 : 1.0;
        return Vec3{ M2.x*(Ms/n2s), M2.y*(Ms/n2s), M2.z*(Ms/n2s) };
    }
    void step(const Vec3& M, double Hy_rf, double Hz_rf, double dt,
              Vec3& M_new_out, double& dMy_dt_out, double& dMz_dt_out) const {
        Vec3 Heff = build_Heff(Hy_rf, Hz_rf);
        Vec3 M_new = llg_step_core(M, Heff, dt);
        dMy_dt_out = (M_new.y - M.y) / dt;
        dMz_dt_out = (M_new.z - M.z) / dt;
        M_new_out = M_new;
    }
};

struct MenDriveSim {
    int N;
    double a, h_l, h_r, dx, dt;
    int ferrite_model;   // 0='JA', 1='LLG', 2='Hybrid'
    int excitation_mode; // 0='magnetic_right', 1='electric_left'
    int n_newton;

    std::vector<double> xA, xB;
    std::vector<uint8_t> leftA, rightA, leftB, rightB;
    int iL_A, iR_A, n_fer, n_left;

    std::vector<double> Ca_e, Cb_e, src_mag_B, src_el_A, Ey, Ez, Hy, Hz;
    double t;

    std::vector<double> M_irr, M_field;         // JA
    std::vector<Vec3> M_llg;                     // LLG
    std::vector<double> hyb_M_irr, hyb_M_field;  // Hybrid (внутренний JA)
    std::vector<Vec3> hyb_M;                     // Hybrid (LLG-часть)
    std::vector<double> last_HzJA;               // Hybrid (неискажённое поле)

    std::unique_ptr<JAEngine> ja;
    std::unique_ptr<LLGEngine> llg;

    // bias_orientation: 0='none', 1='x', 2='y', 3='z', 4='parallel'
    MenDriveSim(int N_, double a_, double h_l_, double h_r_, double dt_frac,
                int excitation_mode_, int ferrite_model_, int bias_orientation, double H0_bias,
                double sigma_m_leak, double sigma_e_left, double sigma_e_ferrite,
                double Ms, double a_JA, double alpha_JA, double k_JA, double c_JA,
                int n_sub, int n_fp, int n_newton_,
                double Ms_llg, double gamma_llg, double alpha_llg, double skin_depth)
        : N(N_), a(a_), h_l(h_l_), h_r(h_r_),
          ferrite_model(ferrite_model_), excitation_mode(excitation_mode_), n_newton(n_newton_)
    {
        double x_min = -(a + h_l), x_max = (a + h_r);
        dx = (x_max - x_min) / N; dt = dt_frac * dx;

        xA.resize(N+1); xB.resize(N);
        leftA.resize(N+1); rightA.resize(N+1); leftB.resize(N); rightB.resize(N);
        for (int i = 0; i <= N; ++i) { xA[i] = x_min + dx*i; leftA[i] = xA[i] < -a; rightA[i] = xA[i] > a; }
        for (int i = 0; i < N; ++i) { xB[i] = x_min + dx*(i + 0.5); leftB[i] = xB[i] < -a; rightB[i] = xB[i] > a; }

        iL_A = 0; iR_A = 0;
        double bestL = 1e300, bestR = 1e300;
        for (int i = 0; i <= N; ++i) {
            double dl = std::fabs(xA[i] + a); if (dl < bestL) { bestL = dl; iL_A = i; }
            double dr = std::fabs(xA[i] - a); if (dr < bestR) { bestR = dr; iR_A = i; }
        }
        n_fer = 0; for (int i = 0; i < N; ++i) if (rightB[i]) n_fer++;
        n_left = 0; for (int i = 0; i < N; ++i) if (leftB[i]) n_left++;

        if (ferrite_model == 0) {
            double H0_eff = (bias_orientation == 4 || bias_orientation == 3) ? H0_bias : 0.0;
            ja.reset(new JAEngine(Ms, a_JA, alpha_JA, k_JA, c_JA, sigma_m_leak, H0_eff, n_sub, n_fp));
        } else if (ferrite_model == 1) {
            int axis = (bias_orientation >= 1 && bias_orientation <= 3) ? bias_orientation : 0;
            llg.reset(new LLGEngine(Ms_llg, gamma_llg, alpha_llg, H0_bias, axis));
        } else {
            ja.reset(new JAEngine(Ms, a_JA, alpha_JA, k_JA, c_JA, sigma_m_leak, 0.0, n_sub, n_fp));
            int axis = (bias_orientation >= 1 && bias_orientation <= 3) ? bias_orientation : 0;
            llg.reset(new LLGEngine(Ms_llg, gamma_llg, alpha_llg, H0_bias, axis));
        }

        Ca_e.resize(N+1); Cb_e.resize(N+1);
        for (int i = 0; i <= N; ++i) {
            double sig = leftA[i] ? sigma_e_left : (rightA[i] ? sigma_e_ferrite : 0.0);
            Ca_e[i] = (1 - 2*M_PI*sig*dt) / (1 + 2*M_PI*sig*dt);
            Cb_e[i] = (dt/dx) / (1 + 2*M_PI*sig*dt);
        }
        src_mag_B.resize(N); src_el_A.resize(N+1);
        for (int i = 0; i < N; ++i) src_mag_B[i] = rightB[i] ? std::exp(-(xB[i]-a)/skin_depth) : 0.0;
        for (int i = 0; i <= N; ++i) src_el_A[i] = leftA[i] ? std::exp(-(-a-xA[i])/skin_depth) : 0.0;
        reset_state();
    }

    void reset_state() {
        Ey.assign(N+1, 0.0); Ez.assign(N+1, 0.0); Hy.assign(N, 0.0); Hz.assign(N, 0.0); t = 0.0;
        if (ferrite_model == 0) { M_irr.assign(n_fer, 0.0); M_field.assign(n_fer, 0.0); }
        else if (ferrite_model == 1) { Vec3 M0 = llg->initial_M(); M_llg.assign(n_fer, M0); }
        else {
            hyb_M_irr.assign(n_fer, 0.0); hyb_M_field.assign(n_fer, 0.0);
            Vec3 M0 = llg->initial_M(); hyb_M.assign(n_fer, M0);
            last_HzJA.assign(n_fer, 0.0);
        }
    }

    double Txx_at(int iA) const {
        double Hy_A = (iA > 0 && iA < N) ? 0.5*(Hy[iA-1]+Hy[iA]) : 0.0;
        double Hz_A = (iA > 0 && iA < N) ? 0.5*(Hz[iA-1]+Hz[iA]) : 0.0;
        return (Ey[iA]*Ey[iA] + Ez[iA]*Ez[iA] + Hy_A*Hy_A + Hz_A*Hz_A) / (8*M_PI);
    }

    double step_once(double omega0, double amp, double ramp_periods, double T) {
        double ramp = std::min(1.0, t / (ramp_periods * T));
        double s = amp * ramp * std::sin(omega0 * t);
        std::vector<double> j_e_y_A(N+1, 0.0), j_m_z_B(N, 0.0);
        if (excitation_mode == 0) { for (int i = 0; i < N; ++i) j_m_z_B[i] = s * src_mag_B[i]; }
        else { for (int i = 0; i <= N; ++i) j_e_y_A[i] = s * src_el_A[i]; }

        std::vector<double> Ey_new(N+1, 0.0), Ez_new(N+1, 0.0);
        for (int i = 1; i < N; ++i) {
            Ey_new[i] = Ca_e[i]*Ey[i] - Cb_e[i]*(Hz[i]-Hz[i-1]) - dt*4*M_PI*j_e_y_A[i];
            Ez_new[i] = Ca_e[i]*Ez[i] + Cb_e[i]*(Hy[i]-Hy[i-1]);
        }
        Ey_new[0]=Ey_new[N]=0.0; Ez_new[0]=Ez_new[N]=0.0;
        double P_e = 0.0;
        for (int i = 0; i <= N; ++i) P_e -= j_e_y_A[i] * 0.5*(Ey[i]+Ey_new[i]) * dx;
        Ey = Ey_new; Ez = Ez_new;

        std::vector<double> rotE_z(N);
        for (int i = 0; i < N; ++i) rotE_z[i] = (Ey[i+1]-Ey[i]) / dx;

        std::vector<double> Hz_new(N, 0.0), Hy_new(N, 0.0);
        for (int i = 0; i < N; ++i) Hy_new[i] = Hy[i] + (dt/dx)*(Ez[i+1]-Ez[i]);
        for (int i = 0; i < N; ++i) { if (!rightB[i]) Hz_new[i] = Hz[i] - dt*rotE_z[i] - dt*4*M_PI*j_m_z_B[i]; }

        if (ferrite_model == 0) {
            int fi = 0;
            for (int i = 0; i < N; ++i) {
                if (!rightB[i]) continue;
                double rhs_extra = dt*rotE_z[i] + dt*4*M_PI*j_m_z_B[i];
                double Hz_f, M_irr_new, M_new;
                ja->implicit_step(M_irr[fi], Hz[i], M_field[fi], rhs_extra, dt, n_newton, 0.25, 1e-5, Hz_f, M_irr_new, M_new);
                Hz_new[i] = Hz_f; M_irr[fi] = M_irr_new; M_field[fi] = M_new; fi++;
            }
        } else if (ferrite_model == 1) {
            int fi = 0;
            for (int i = 0; i < N; ++i) {
                if (!rightB[i]) continue;
                Vec3 M_new; double dMy_dt, dMz_dt;
                llg->step(M_llg[fi], Hy[i], Hz[i], dt, M_new, dMy_dt, dMz_dt);
                M_llg[fi] = M_new;
                Hy_new[i] = Hy[i] + (dt/dx)*(Ez[i+1]-Ez[i]) - dt*4*M_PI*dMy_dt;
                Hz_new[i] = Hz[i] - dt*rotE_z[i] - dt*4*M_PI*j_m_z_B[i] - dt*4*M_PI*dMz_dt;
                fi++;
            }
        } else {
            int fi = 0;
            for (int i = 0; i < N; ++i) {
                if (!rightB[i]) continue;
                double rhs_extra_z = dt*rotE_z[i] + dt*4*M_PI*j_m_z_B[i];
                double Hz_after_JA, M_irr_new, M_field_new;
                ja->implicit_step(hyb_M_irr[fi], Hz[i], hyb_M_field[fi], rhs_extra_z, dt, n_newton, 0.25, 1e-5, Hz_after_JA, M_irr_new, M_field_new);
                hyb_M_irr[fi] = M_irr_new; hyb_M_field[fi] = M_field_new;
                last_HzJA[fi] = Hz_after_JA;
                Vec3 M_new; double dMy_dt, dMz_dt;
                llg->step(hyb_M[fi], Hy[i], Hz_after_JA, dt, M_new, dMy_dt, dMz_dt);
                hyb_M[fi] = M_new;
                Hy_new[i] = Hy[i] + (dt/dx)*(Ez[i+1]-Ez[i]) - dt*4*M_PI*dMy_dt;
                Hz_new[i] = Hz_after_JA - dt*4*M_PI*dMz_dt;
                fi++;
            }
        }
        double P_m = 0.0;
        for (int i = 0; i < N; ++i) P_m -= j_m_z_B[i] * 0.5*(Hz[i]+Hz_new[i]) * dx;
        Hy = Hy_new; Hz = Hz_new; t += dt;
        return P_e + P_m;
    }

    bool fields_finite() const { for (double v : Hz) if (!std::isfinite(v)) return false; return true; }
};

// ---------------- C API для ctypes ----------------
extern "C" {

void* mendrive_create(int N, double a, double h_l, double h_r, double dt_frac,
                       int excitation_mode, int ferrite_model, int bias_orientation, double H0_bias,
                       double sigma_m_leak, double sigma_e_left, double sigma_e_ferrite,
                       double Ms, double a_JA, double alpha_JA, double k_JA, double c_JA,
                       int n_sub, int n_fp, int n_newton,
                       double Ms_llg, double gamma_llg, double alpha_llg, double skin_depth) {
    return new MenDriveSim(N, a, h_l, h_r, dt_frac, excitation_mode, ferrite_model,
                            bias_orientation, H0_bias, sigma_m_leak, sigma_e_left, sigma_e_ferrite,
                            Ms, a_JA, alpha_JA, k_JA, c_JA, n_sub, n_fp, n_newton,
                            Ms_llg, gamma_llg, alpha_llg, skin_depth);
}
void mendrive_destroy(void* handle) { delete static_cast<MenDriveSim*>(handle); }
int mendrive_get_n_fer(void* handle) { return static_cast<MenDriveSim*>(handle)->n_fer; }
double mendrive_get_dt(void* handle) { return static_cast<MenDriveSim*>(handle)->dt; }
double mendrive_get_dx(void* handle) { return static_cast<MenDriveSim*>(handle)->dx; }

// Прогон nsteps шагов, запись с record_start. probe_idx -- индекс внутри
// феррита для скалярных Hn/M/HzJA/MzJA. Буферы out_* должны быть выделены
// вызывающей стороной размером >= (nsteps - record_start).
int mendrive_run(void* handle, double omega0, int nsteps, int record_start,
                  double amp, double ramp_periods, int probe_idx,
                  double* out_t, double* out_dTxx, double* out_Hn, double* out_P,
                  double* out_Mx, double* out_My, double* out_Mz,
                  double* out_HzJA, double* out_MzJA, int* out_blew_up) {
    MenDriveSim* sim = static_cast<MenDriveSim*>(handle);
    double T = 2*M_PI/omega0;
    int rec_count = 0;
    *out_blew_up = 0;
    for (int n = 0; n < nsteps; ++n) {
        double P_step = sim->step_once(omega0, amp, ramp_periods, T);
        if (n >= record_start) {
            out_t[rec_count] = sim->t;
            out_P[rec_count] = P_step / sim->dt;
            out_dTxx[rec_count] = sim->Txx_at(sim->iR_A) - sim->Txx_at(sim->iL_A);
            int fi = 0, gi = -1;
            for (int i = 0; i < sim->N; ++i) { if (sim->rightB[i]) { if (fi == probe_idx) { gi = i; break; } fi++; } }
            out_Hn[rec_count] = (gi >= 0) ? sim->Hz[gi] : 0.0;
            if (sim->ferrite_model == 0) {
                out_My[rec_count] = sim->M_field[probe_idx];
                out_Mx[rec_count] = 0.0; out_Mz[rec_count] = 0.0;
                out_HzJA[rec_count] = 0.0; out_MzJA[rec_count] = 0.0;
            } else if (sim->ferrite_model == 1) {
                Vec3 M = sim->M_llg[probe_idx];
                out_Mx[rec_count] = M.x; out_My[rec_count] = M.y; out_Mz[rec_count] = M.z;
                out_HzJA[rec_count] = 0.0; out_MzJA[rec_count] = 0.0;
            } else {
                Vec3 M = sim->hyb_M[probe_idx];
                out_Mx[rec_count] = M.x; out_My[rec_count] = M.y; out_Mz[rec_count] = M.z;
                out_HzJA[rec_count] = sim->last_HzJA[probe_idx];
                out_MzJA[rec_count] = sim->hyb_M_field[probe_idx];
            }
            rec_count++;
        }
        if (!sim->fields_finite()) { *out_blew_up = 1; break; }
    }
    return rec_count;
}

} // extern "C"