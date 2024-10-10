#include <stdio.h>
#include <math.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>


// The following program solves the second-order nonlinear pendulum oscillator equation
// U = m * g * R * cos(phi)
// H = R*g*m*cos(phi) + 1/2*p^2/m
// q = phi
// dot_p = R*g*m*sin(phi)
// dot_q = p/m

// y = {p,q}
// y[0] = p
// y[1] = q

// f = {dot_p, dot_q}
// f[0] = R*g*m*sin(phi)
// f[1] = p/m

// f[0] = R*g*m*sin(y[1])
// f[1] = y[0]/m

// J[0,0] = 0
// J[0,1] = R*g*m*cos(phi) = R*g*m*cos(y[1])
// J[1,0] = 1/m
// J[1,1[ = 0

struct pendulum
{
    double m;
    double g;
    double R;
};


int
func (double t, const double y[], double f[],
      void *params)
{
  (void)(t); /* avoid unused parameter warning */
  struct pendulum p = *(struct pendulum *)params;
  f[0] = p.R*p.g*p.m*sin(y[1]);
  f[1] = y[0]/p.m;
  return GSL_SUCCESS;
}

int
jac (double t, const double y[], double *dfdy,
     double dfdt[], void *params)
{
  (void)(t); /* avoid unused parameter warning */
  struct pendulum p = *(struct pendulum *)params;
  gsl_matrix_view dfdy_mat
    = gsl_matrix_view_array (dfdy, 2, 2);
  gsl_matrix * m = &dfdy_mat.matrix;
  gsl_matrix_set (m, 0, 0, 0.0);
  gsl_matrix_set (m, 0, 1, p.R*p.g*p.m*cos(y[1]));
  gsl_matrix_set (m, 1, 0, 1.0/p.m);
  gsl_matrix_set (m, 1, 1, 0.0);
  dfdt[0] = 0.0;
  dfdt[1] = 0.0;
  return GSL_SUCCESS;
}

int
main (void)
{
  struct pendulum p;
  p.m = 1.0;
  p.g = 1.0;
  p.R = 1.0;
  gsl_odeiv2_system sys = {func, jac, 2, &p};

  gsl_odeiv2_driver * d =
    gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
                                  1e-6, 1e-6, 0.0);
  int i;
  double t = 0.0, t1 = 100.0;
  double y[2] = { 0.0001, 0.0 };

  for (i = 1; i <= 100; i++)
    {
      double ti = i * t1 / 10000.0;
      int status = gsl_odeiv2_driver_apply (d, &t, ti, y);

      if (status != GSL_SUCCESS)
        {
          printf ("error, return value=%d\n", status);
          break;
        }

      printf ("%.5e %.5e %.5e\n", t, y[0], y[1]);
    }

  gsl_odeiv2_driver_free (d);
  return 0;
}