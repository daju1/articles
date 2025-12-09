#define INCLUDE_CUBA_H
#include "calc_RO.h"

int Integrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata);

/*********************************************************************/

#define NCOMP 9

int integrate(
    double R0,    /* Радиус сферы */
    double v0,    /* Начальная продольная скорость */
    double a,     /* Продольное ускорение */
    double c,     /* Скорость света */
    double t,     /* Текущее время */
    double t0,    /* Начальное время */
    int use_delay,
    int use_lorentz_factor,
    int use_lorentz_general_factor,
    int use_fermi_factor_O,
    int use_fermi_factor,
    int use_fermi_general_factor,
    int use_fast_integrand,
    cubareal* integral, cubareal* error, cubareal* prob);

int main() {
  int comp, nregions, neval, fail;
  cubareal integral[NCOMP], error[NCOMP], prob[NCOMP];

  /* Параметры задачи */
  double R0 = 1;      /* радиус сферы (половина R1) */
  double v0 = 0;
  double t0 = 0;
  double t  = 0;
  double c  = 1;       /* скорость света в м/с */
  double a  = 0.00000001;  /* продольное ускорение*/

  int use_delay                  = 0;
  int use_lorentz_factor         = 0;
  int use_lorentz_general_factor = 1;
  int use_fermi_factor_O         = 1;
  int use_fermi_factor           = 0;
  int use_fermi_general_factor   = 0;
  int use_fast_integrand         = 1;

  integrate(
    R0,    /* Радиус сферы */
    v0,    /* Начальная продольная скорость */
    a,     /* Продольное ускорение */
    c,     /* Скорость света */
    t,     /* Текущее время */
    t0,    /* Начальное время */
    use_delay,
    use_lorentz_factor,
    use_lorentz_general_factor,
    use_fermi_factor_O,
    use_fermi_factor,
    use_fermi_general_factor,
    use_fast_integrand,
    integral, error, prob);

  double v_z = v0 + a * (t - t0);
  double v_c =  v_z / c;  /* отношение скорости заряда к скорости света*/
  printf("v_c = %f\n", v_c);

  printf("VEGAS RESULT:\tneval %d\tfail %d\n",
    neval, fail);
  for( comp = 0; comp < NCOMP; ++comp )
    printf("VEGAS RESULT:\t%.8f +- %.8f\tp = %.3f\n",
      (double)integral[comp], (double)error[comp], (double)prob[comp]);

  return 0;
}
