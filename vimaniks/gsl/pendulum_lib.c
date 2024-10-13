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
    double y[2];
    double t;
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

static int
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

static gsl_odeiv2_driver * d;
static struct pendulum p;

void init (double m, double g, double R, double p0, double q0, double t0)
{
  p.m = m;
  p.g = g;
  p.R = R;
  p.y[0] = p0;
  p.y[1] = q0;
  p.t = t0;
}

void alloc()
{
  static gsl_odeiv2_system sys;
  sys.function = &func;
  sys.jacobian = &jac;
  sys.dimension = 2;
  sys.params = &p;

  d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd,
                                  1e-6, 1e-6, 0.0);
}

int apply(double ti, double * pt, double *pMomenta, double *pq)
{
  int status = gsl_odeiv2_driver_apply (d, &p.t, ti, p.y);

  if (status != GSL_SUCCESS)
  {
    printf ("error, return value=%d\n", status);
  }

  //printf ("%.5e %.5e %.5e\n", p.t, p.y[0], p.y[1]);
  *pt = p.t;
  *pMomenta = p.y[0];
  *pq = p.y[1];
  return status;
}

void release()
{
  gsl_odeiv2_driver_free (d);
}
