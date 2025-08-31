#define INCLUDE_CUBA_H
#include "calc_RO.h"

int Integrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata);

/*********************************************************************/

#define NCOMP 9

int integrate(
    double R0,       /* Радиус сферы */
    double rho0,     /* Радиус орбиты */
    double omega,    /* Угловая скорость */
    double c,        /* Скорость света */
    int use_delay,
    int use_lorentz_factor,
    int use_lorentz_general_factor,
    int use_fermi_factor,
    int use_fermi_general_factor,
    int use_fast_integrand,
    cubareal* integral, cubareal* error, cubareal* prob);

int main() {
  int comp, nregions, neval, fail;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

  /* Параметры задачи */

  double R0 = 1;      /* радиус сферы (половина R1) */
  double rho0 = 100;  /* расстояние от оси орбитального движения */
  double c = 1;       /* скорость света в м/с */
  // m = 1
  double omega = sqrt(5 * R0 * Sq(c) / (3*2*4* Cb(rho0)));  /* угловая скорость орбитального движения*/
  // m = 4/3
  //double omega = sqrt(3*5 * R0 * Sq(c) / (4*3*2*4* Cb(rho0)));  /* угловая скорость орбитального движения*/
  double v_c = omega * rho0 / c;  /* отношение скорости заряда к скорости света*/

  int use_delay                  = 0;
  int use_lorentz_factor         = 0;
  int use_lorentz_general_factor = 0;
  int use_fermi_factor           = 0;
  int use_fermi_general_factor   = 0;
  int use_fast_integrand         = 0;

  integrate(
    R0,    /* Радиус сферы */
    rho0,  /* Начальная продольная скорость */
    omega, /* Продольное ускорение */
    c,     /* Скорость света */
    use_delay,
    use_lorentz_factor,
    use_lorentz_general_factor,
    use_fermi_factor,
    use_fermi_general_factor,
    use_fast_integrand,
    integral, error, prob);

  printf("VEGAS RESULT:\tneval %d\tfail %d\n",
    neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);

  return 0;
}
