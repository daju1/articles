
#include "calc_RO.h"

// the charge distribution of the proton
// g(r) = exp(-r^2/r__0^2)/(r__0^3*sqrt(pi)^3)
static inline cubareal rho_q (cubareal r0, cubareal r)
{
    return exp(-Sq(r/r0)) / Cb(r0*sqrt(M_PI));
}

// интегрирование по координатам заряда источника потенциала
static inline cubareal Iq (cubareal r0, cubareal theta_a, cubareal ra, cubareal phi_q, cubareal theta_q, cubareal rq )
{
    cubareal R = R0 (ra, theta_a, rq, theta_q, phi_q);
    return rho_q(r0, rq) * Sq(rq) * sin(theta_q) / R;
}

// интегрирование по координатам точек наблюдения
static inline cubareal Ia (cubareal r0, cubareal theta_a, cubareal ra, cubareal phi_q, cubareal theta_q, cubareal rq)
{
    return 2 * M_PI * rho_q(r0, ra) * Iq(r0, theta_a, ra, phi_q, theta_q, rq) * sin(theta_a) * Sq(ra) ;
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

    // http://www.actaphys.uj.edu.pl/fulltext?series=Reg&vol=30&page=119
    // the rms charge radius of the proton being
    // rp_rms = 0.8
    // r0 = 2/3 * rp_rms
    // the charge distribution of the proton
    // g(r) = exp(-r^2/r__0^2)/(r__0^3*sqrt(pi)^3)

    cubareal rp_rms = 0.8, r0 = 2.0/3.0 * rp_rms;
    cubareal q = 1.0;
    cubareal rb = rp_rms * 3;

    f = rb * rb * M_PI * (2*M_PI) * M_PI * Ia (r0, theta_a * M_PI, ra * rb, phi_q * (2*M_PI), theta_q * M_PI, rq * rb );

    return 0;
}
