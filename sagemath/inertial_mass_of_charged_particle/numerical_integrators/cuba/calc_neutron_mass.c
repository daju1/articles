#include "calc_RO.h"

// http://www.actaphys.uj.edu.pl/fulltext?series=Reg&vol=30&page=119

// the rms value of g__n
// r2 := -0.113

// The third parameter r1 is a scaling parameter, which is necessary to define
// a dimensionless quantity (r/r1) in the Gaussian exponent. The results of
// QCD-calculations of the charge density distribution inside the neutron [2]
// are best reproduced by choosing:
// r1 = 0.71*sqrt(2/5) # fm

// the charge density distribution within the neutron
// gn(r) = (-2/3)*(r2 / (r1^2 * (r1*sqrt(pi))^3)) * (r/r1)^2 * (1 - (2/5)*r^2/r1^2) * exp(-r^2/r1^2)

//     rho_q = lambda r1, r2, r : (-2/3)*(r2 / (r1^2 * (r1*sqrt(pi))^3)) * (r/r1)^2 * (1 - (2/5)*r^2/r1^2) * exp(-r^2/r1^2)
static inline cubareal rho_q (cubareal r1, cubareal r2, cubareal r)
{
    return (-2.0/3.0)*(r2 / (Sq(r1) * Cb(r1*sqrt(M_PI)))) * Sq(r/r1) * (1 - (2.0/5.0)*Sq(r/r1)) * exp(-Sq(r/r1));
}


// интегрирование по координатам заряда источника потенциала
static inline cubareal Iq (cubareal r1, cubareal r2, cubareal theta_a, cubareal ra, cubareal phi_q, cubareal theta_q, cubareal rq )
{
    cubareal R = R0 (ra, theta_a, rq, theta_q, phi_q);
    return rho_q(r1, r2, rq) * Sq(rq) * sin(theta_q) / R;
}

// интегрирование по координатам точек наблюдения
static inline cubareal Ia (cubareal r1, cubareal r2, cubareal theta_a, cubareal ra, cubareal phi_q, cubareal theta_q, cubareal rq)
{
    return 2 * M_PI * rho_q(r1, r2, ra) * Iq(r1, r2, theta_a, ra, phi_q, theta_q, rq) * sin(theta_a) * Sq(ra) ;
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

    // the rms value of g__n
    // r2 := -0.113

    // The third parameter r1 is a scaling parameter, which is necessary to define
    // a dimensionless quantity (r/r1) in the Gaussian exponent. The results of
    // QCD-calculations of the charge density distribution inside the neutron [2]
    // are best reproduced by choosing:
    // r1 = 0.71*sqrt(2/5) # fm

    cubareal r2 = -0.113, r1 = 0.71*sqrt(2.0/5.0);
    cubareal q = 1.0;
    cubareal rb = -r2 * 3;

    f = rb * rb * M_PI * (2*M_PI) * M_PI * Ia (r1, r2, theta_a * M_PI, ra * rb, phi_q * (2*M_PI), theta_q * M_PI, rq * rb );

    return 0;
}

