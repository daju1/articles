/*
	demo-c.c
		test program for the Cuba library
		last modified 13 Mar 15 th
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if REALSIZE == 16
#include "cubaq.h"
#elif REALSIZE == 10
#include "cubal.h"
#else
#include "cuba.h"
#endif


static inline cubareal Sq(cubareal x) {
  return x*x;
}

//z - координаты заряда и точки наблюдения
static inline cubareal zq (cubareal rq, cubareal theta_q)
{
    return (rq)*cos(theta_q);
}
static inline cubareal za (cubareal ra, cubareal theta_a)
{
    return (ra)*cos(theta_a);
}

// введём вспомогательные переменные - цилиндрический радиус
// как координата r точки при переходе в цилиндрическую систему коорденат с тем же направлением оси z
static inline cubareal rcq (cubareal rq, cubareal theta_q )
{
    return (rq)*sin(theta_q);
}
static inline cubareal rca (cubareal ra, cubareal theta_a )
{
    return (ra)*sin(theta_a);
}

// выражение для расстояния между точкой заряда и точкой наблюдения примет вид
static inline cubareal R0 (cubareal ra, cubareal theta_a, cubareal rq, cubareal theta_q, cubareal phi_q)
{
    return sqrt(Sq(rca(ra, theta_a)) + Sq(rcq(rq, theta_q)) + Sq(za(ra, theta_a)-zq(rq, theta_q)) - 2*rca(ra, theta_a)*rcq(rq, theta_q)*cos(phi_q));
}

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

static int Integrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {

    // theta_a, ra, phi_q, theta_q, rq

    #define theta_a xx[0]
    #define ra      xx[1]
    #define phi_q   xx[2]
    #define theta_q xx[3]
    #define rq      xx[4]

    #define f ff[0]

    cubareal r0 = 0.1;
    cubareal q = 1.0;

    f = r0 * r0 * M_PI * (2*M_PI) * M_PI * Ia (r0, q, theta_a * M_PI, ra * r0, phi_q * (2*M_PI), theta_q * M_PI, rq * r0 );

    return 0;
}

/*********************************************************************/

#define NDIM 5
#define NCOMP 1
#define USERDATA NULL
#define NVEC 1
#define EPSREL 1e-3
#define EPSABS 1e-12
#define VERBOSE 2
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 50000

#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

#define NNEW 1000
#define NMIN 2
#define FLATNESS 25.

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

#define KEY 0

int main() {
  int comp, nregions, neval, fail;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

#if 0
  printf("-------------------- Vegas test --------------------\n");

  Vegas(NDIM, NCOMP, Integrand, USERDATA, NVEC,
    EPSREL, EPSABS, VERBOSE, SEED,
    MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
    GRIDNO, STATEFILE, SPIN,
    &neval, &fail, integral, error, prob);

  printf("VEGAS RESULT:\tneval %d\tfail %d\n",
    neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);
#endif

#if 1
  printf("\n-------------------- Suave test --------------------\n");

  Suave(NDIM, NCOMP, Integrand, USERDATA, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST, SEED,
    MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);

  printf("SUAVE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("SUAVE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);
#endif

#if 0
  printf("\n------------------- Divonne test -------------------\n");

  Divonne(NDIM, NCOMP, Integrand, USERDATA, NVEC,
    EPSREL, EPSABS, VERBOSE, SEED,
    MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
    BORDER, MAXCHISQ, MINDEVIATION,
    NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);

  printf("DIVONNE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("DIVONNE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);
#endif

#if 0
  printf("\n-------------------- Cuhre test --------------------\n");

  Cuhre(NDIM, NCOMP, Integrand, USERDATA, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);

  printf("CUHRE RESULT:\tnregions %d\tneval %d\tfail %d\n",
    nregions, neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("CUHRE RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);
#endif

  return 0;
}

