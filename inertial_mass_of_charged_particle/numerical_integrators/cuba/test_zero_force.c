/* test_zero_force.c: попарная антисимметрия кулоновского поля при omega -> 0 */
#include "calc_transverse_sphere_mass_film.c"
/* Внимание: после Integrand() в файле остаются макросы psi_a, theta_a, ra,
   psi_q, theta_q, rq, observed_ratio_*. Эти имена здесь использовать нельзя. */

static double u01(void) { return rand() / ((double)RAND_MAX + 1.0); }

static void rand_point(double R0, double p[3])   /* p = {r, theta, psi} */
{
    p[0] = R0 * cbrt(u01());
    p[1] = acos(1.0 - 2.0 * u01());
    p[2] = 2.0 * M_PI * u01();
}

/* E_total, как его видит Integrand: компоненты (x, y, z) после rot_to_lab */
static void field(const ProblemParams *P, const double A[3], const double Q[3], double E[3])
{
    double e1r, e1p, e1z, e2r, e2p, e2z;
    if (P->use_delay)
        compute_electric_field_with_delay(A[0], A[1], A[2], Q[0], Q[1], Q[2], P,
            &e1r, &e1p, &e1z, &e2r, &e2p, &e2z, &E[0], &E[1], &E[2]);
    else
        compute_electric_field(A[0], A[1], A[2], Q[0], Q[1], Q[2], P,
            &e1r, &e1p, &e1z, &e2r, &e2p, &e2z, &E[0], &E[1], &E[2]);
}

static double phi_of(const ProblemParams *P, const double p[3])  /* только для use_delay = 0 */
{
    double x, y, z;
    get_source_position(0, p[0], p[1], p[2], P->rho0, P->omega, P->c, &x, &y, &z, 0, 0);
    return atan2(y, x);
}

static double norm3(const double v[3]) { return sqrt(v[0]*v[0] + v[1]*v[1] + v[2]*v[2]); }

static double pair_residual(const ProblemParams *P, const double A[3], const double Q[3])
{
    double Eaq[3], Eqa[3], sum[3];
    field(P, A, Q, Eaq); field(P, Q, A, Eqa);
    for (int k = 0; k < 3; k++) sum[k] = Eaq[k] + Eqa[k];
    double s = norm3(Eaq) + norm3(Eqa);
    return (s < 1e-300) ? -1.0 : norm3(sum) / s;
}

int main(void)
{
    const double rho0s[] = {0.92, 2.0, 100.0};
    const int NPAIRS = 20000;
    ProblemParams P = {0};
    P.R0 = 1.0; P.c = 1.0; P.omega = 1e-8;   /* v/c <= 1e-6 даже при rho0 = 100 */
    srand(12345);

    printf("=== T1: E(a<-q) + E(q<-a) = 0 при omega -> 0 ===\n");
    for (int ir = 0; ir < 3; ir++) for (int d = 0; d < 2; d++) {
        P.rho0 = rho0s[ir]; P.use_delay = d;
        for (int rot = 1; rot >= 0; rot--) {
            P.use_rot_to_lab = rot;
            srand(12345);                       /* одни и те же пары в обоих режимах */
            double mean = 0, mx = 0; int n = 0;
            for (int i = 0; i < NPAIRS; i++) {
                double A[3], Q[3];
                rand_point(P.R0, A); rand_point(P.R0, Q);
                double r = pair_residual(&P, A, Q);
                if (r < 0) continue;
                mean += r; if (r > mx) mx = r; n++;
            }
            printf(" rho0=%6.2f delay=%d rot=%d : mean=%.1e max=%.1e  %s\n",
                   P.rho0, d, rot, mean/n, mx,
                   rot ? ((mean/n < 1e-4) ? "OK" : "FAIL") : "(контроль: невязка должна быть большой)");
        }
    }

    printf("\n=== T2: релятивистские факторы Ферми -> нерелятивистские при v << c ===\n");
    {
        double Rr = 0.4, Rp = 0.3, Rz = 0.2, a = -1e-2, vph = 1e-3, c = 1.0;
        double g  = fermi_correction_general  (0, 0, Rr, Rp, Rz, a, 0, 0, 0, vph, 0, c);
        double g0 = fermi_correction          (Rr, Rp, Rz, a, 0, 0, c);
        double gO = fermi_correction_general_O(0, 0, Rr, Rp, Rz, a, 0, 0, 0, vph, 0, c);
        double O0 = fermi_correction_O        (Rr, Rp, Rz, a, 0, 0, c);
        printf(" general   vs P' : %.6e vs %.6e  (отн. разность %.1e)\n", g,  g0, fabs(g-g0)/fabs(g0));
        printf(" general_O vs O  : %.6e vs %.6e  (отн. разность %.1e)\n", gO, O0, fabs(gO-O0)/fabs(O0));
    }
    return 0;
}
