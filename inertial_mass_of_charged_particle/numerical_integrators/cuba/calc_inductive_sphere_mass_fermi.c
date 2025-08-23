/*
	calc_inductive_sphere_mass_fermi.c
*/

#include "calc_RO.h"

// the charge distribution
static inline cubareal rho_q (cubareal r0, cubareal q)
{
    return 3 * (q) / (4*M_PI*(r0)*(r0)*(r0));
}

// интегрирование по координатам заряда источника потенциала
// E2 - индуктивная компонента массы
// В приближении малых скоростей ${}^{v} \big / {}_{c}\ll 1$
// и малых ускорений $a{{r}_{0}}\ll {{c}^{2}}$ и при игнорировании запаздывания
// Электрическое поле самоиндукции ($z$ компонента)
/*
\[{E}_{2}=\\
\int\limits_{{{r}_{q}}}\int\limits_{{{\varphi}_{q}}}\int\limits_{{{\theta}_{q}}}\\
{\left\{ -\frac{{a_z}}{{{c}^{2}}{{{R}_{0}}}} \right\}\\
{\rho \left( {{r}_{q}} \right){{r}_{q}}^{2}\sin \left( {{\theta }_{q}} \right)}\ }d{{\theta }_{q}}d{{\varphi }_{q}}d{{r}_{q}}\]
*/
static inline cubareal Iq (cubareal r0, cubareal q, cubareal theta_a, cubareal ra, cubareal phi_q, cubareal theta_q, cubareal rq, cubareal c, cubareal a)
{
    cubareal za_minus_zq = ra*cos(theta_a) - rq*cos(theta_q);
    // В приближении малых скоростей ${}^{v}/{}_{c}\ll 1$
    // но при учете запаздывания
    cubareal R = R_a (ra, theta_a, rq, theta_q, phi_q, c, a);
    return - rho_q(r0, q) * Sq(rq) * sin(theta_q) / R * (1.0 + za_minus_zq * a / Sq(c));
}

// интегрирование по координатам точек наблюдения
/*static inline*/ cubareal Ia (cubareal r0, cubareal q, cubareal theta_a, cubareal ra, cubareal phi_q, cubareal theta_q, cubareal rq, cubareal c, cubareal a)
{
    return 2 * M_PI * rho_q(r0, q) * Iq(r0, q, theta_a, ra, phi_q, theta_q, rq, c, a) * sin(theta_a) * Sq(ra) ;
}

int Integrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

    // theta_a, ra, phi_q, theta_q, rq

    #define theta_a xx[0]
    #define ra      xx[1]
    #define phi_q   xx[2]
    #define theta_q xx[3]
    #define rq      xx[4]

    #define f ff[0]

    if ( fabs(theta_a - theta_q) < 1e-12 && fabs(ra - rq) < 1e-12 && fabs(phi_q - 0.0) < 1e-12)
    {
        printf ("theta_a = %e ", theta_a);
        printf ("theta_q = %e ", theta_q);
        printf ("ra = %e ", ra);
        printf ("rq = %e ", rq);
        printf ("phi_q = %e\n", phi_q);

        return -999;
    }

    cubareal r0 = 1.0;
    cubareal q = 1.0;
    cubareal c = 1.0;
    cubareal a = 5.0;

    f = r0 * r0 * M_PI * (2*M_PI) * M_PI * Ia (r0, q, theta_a * M_PI, ra * r0, phi_q * (2*M_PI), theta_q * M_PI, rq * r0 , c, a);

    return 0;
}

