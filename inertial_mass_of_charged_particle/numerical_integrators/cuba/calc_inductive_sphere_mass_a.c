/*
	calc_inductive_sphere_mass_a.c
*/

#include "calc_RO.h"

/* Структура для передачи параметров задачи */
typedef struct {
    double r0;    /* Радиус сферы */
    double q;     /* заряд */
    double a;     /* Продольное ускорение */
    double c;     /* Скорость света */
} ProblemParams;


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
    // В приближении малых скоростей ${}^{v}/{}_{c}\ll 1$
    // но при учете запаздывания
    cubareal R = R_a (ra, theta_a, rq, theta_q, phi_q, c, a);
    return - rho_q(r0, q) * Sq(rq) * sin(theta_q) / R;
}

// интегрирование по координатам точек наблюдения
static inline cubareal Ia (cubareal r0, cubareal q, cubareal theta_a, cubareal ra, cubareal phi_q, cubareal theta_q, cubareal rq, cubareal c, cubareal a)
{
    return 2 * M_PI * rho_q(r0, q) * Iq(r0, q, theta_a, ra, phi_q, theta_q, rq, c, a) * sin(theta_a) * Sq(ra) ;
}

int Integrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {
    ProblemParams *params = (ProblemParams *)userdata;

    // theta_a, ra, phi_q, theta_q, rq

    #define theta_a xx[0]
    #define ra      xx[1]
    #define phi_q   xx[2]
    #define theta_q xx[3]
    #define rq      xx[4]

    #define f ff[0]

    if ( fabs(theta_a - theta_q) < 1e-18 && fabs(ra - rq) < 1e-18 && fabs(phi_q - 0.0) < 1e-18)
    {
        printf ("theta_a = %e ", theta_a);
        printf ("theta_q = %e ", theta_q);
        printf ("ra = %e ", ra);
        printf ("rq = %e ", rq);
        printf ("phi_q = %e\n", phi_q);

        return -999;
    }

#if 1
    cubareal r0 = params->r0;
    cubareal q = params->q;
    cubareal c = params->c;
    cubareal a = params->a;
#else
    cubareal r0 = 1.0;
    cubareal q = 1.0;
    cubareal c = 1.0;
    cubareal a = 5.0;
#endif

    f = r0 * r0 * M_PI * (2*M_PI) * M_PI * Ia (r0, q, theta_a * M_PI, ra * r0, phi_q * (2*M_PI), theta_q * M_PI, rq * r0 , c, a);

    return 0;
}

#ifdef INCLUDE_CUBA_H

/*********************************************************************/

#define NDIM 5
#define NCOMP 1
#define USERDATA NULL
#define NVEC 1
#if 1
#define EPSREL 1e-16
#define EPSABS 1e-32
#else
#define EPSREL 1e-3
#define EPSABS 1e-12
#endif
#define VERBOSE 2
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 500000

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

int integrate(
    double r0,    /* Радиус сферы */
    double q,     /* заряд */
    double a,     /* Продольное ускорение */
    double c,     /* Скорость света */
    cubareal* integral, cubareal* error, cubareal* prob
) {
    int comp, nregions, neval, fail;
    // cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

    /* Параметры задачи */
    ProblemParams params;
    params.r0 = r0; /* радиус сферы */
    params.q = q;   /* заряд */
    params.c = c;   /* скорость света в м/с */
    params.a = a;   /* продольное ускорение */

#if 1
  printf("-------------------- Vegas test --------------------\n");

  Vegas(NDIM, NCOMP, Integrand, &params, NVEC,
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

#if 0
  printf("\n-------------------- Suave test --------------------\n");

  Suave(NDIM, NCOMP, Integrand, &params, NVEC,
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

  Divonne(NDIM, NCOMP, Integrand, &params, NVEC,
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

  Cuhre(NDIM, NCOMP, Integrand, &params, NVEC,
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

#endif // #ifdef INCLUDE_CUBA_H

