/*
	demo-c.c
		test program for the Cuba library
		last modified 13 Mar 15 th
*/

#include "calc_RO.h"

// the charge distribution
static inline cubareal rho_q (cubareal r0, cubareal q)
{
    return 3 * (q) / (4*M_PI*(r0)*(r0)*(r0));
}

// интегрирование по координатам заряда источника потенциала
static inline cubareal Iq (cubareal r0, cubareal q, cubareal theta_a, cubareal ra, cubareal phi_q, cubareal theta_q, cubareal rq )
{
    cubareal R = R0 (ra, theta_a, rq, theta_q, phi_q);
    return rho_q(r0, q) * Sq(rq) * sin(theta_q) / R;
}

// интегрирование по координатам точек наблюдения
static inline cubareal Ia (cubareal r0, cubareal q, cubareal theta_a, cubareal ra, cubareal phi_q, cubareal theta_q, cubareal rq)
{
    return 2 * M_PI * rho_q(r0, q) * Iq(r0, q, theta_a, ra, phi_q, theta_q, rq) * sin(theta_a) * Sq(ra) ;
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

    f = r0 * r0 * M_PI * (2*M_PI) * M_PI * Ia (r0, q, theta_a * M_PI, ra * r0, phi_q * (2*M_PI), theta_q * M_PI, rq * r0 );

    return 0;
}

