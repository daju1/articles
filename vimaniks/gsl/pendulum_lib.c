#include <stdio.h>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>

#include "pendulum_lib.h"

#define Sq(x) ((x)*(x))
#define Cu(x) ((x)*(x)*(x))

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
    double xc, yc;
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

void init (double m, double g, double R, double p0, double q0, double t0, double xc, double yc)
{
  p.m = m;
  p.g = g;
  p.R = R;
  p.y[0] = p0;
  p.y[1] = q0;
  p.t = t0;
  p.xc = xc;
  p.yc = yc;
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

int logging = 0;

int apply(double ti, double * pt, double *pMomenta, double *pq,
          double * pdot_q, double * pddot_q, double * pdddot_q,
          double * psx, double * psy,
          double * pvx, double * pvy,
          double * pwx, double * pwy,
          double * pdot_wx, double * pdot_wy
         )
{
  if (d->h * (ti - p.t) < 0.0)
  {
    d->h *= -1;
  }

  int status = gsl_odeiv2_driver_apply (d, &p.t, ti, p.y);

  if (status != GSL_SUCCESS)
  {
    printf ("error, return value=%d\n", status);
  }

  //printf ("%.5e %.5e %.5e\n", p.t, p.y[0], p.y[1]);
  if (pt != NULL) *pt = p.t;
  *pMomenta = p.y[0];
  *pq = p.y[1];

  *pdot_q   = *pMomenta / p.m;                  // dphi/dt
  *pddot_q  = p.R * p.g * sin(*pq);             // d2phi/dt2
  *pdddot_q = p.R * p.g * cos(*pq) * (*pdot_q); // d3phi/dt3

  *psx = p.xc + p.R*cos(*pq);
  *psy = p.yc + p.R*sin(*pq);

  *pvx = - p.R*sin(*pq) * (*pdot_q);
  *pvy = + p.R*cos(*pq) * (*pdot_q);

  *pwx = - p.R*cos(*pq) * Sq((*pdot_q)) - p.R*sin(*pq) * (*pddot_q);
  *pwy = - p.R*sin(*pq) * Sq((*pdot_q)) + p.R*cos(*pq) * (*pddot_q);

  *pdot_wx = + p.R*sin(*pq) * Cu((*pdot_q)) - 3 * p.R*cos(*pq) * (*pdot_q) * (*pddot_q) - p.R*sin(*pq) * (*pdddot_q);
  *pdot_wy = - p.R*cos(*pq) * Cu((*pdot_q)) - 3 * p.R*sin(*pq) * (*pdot_q) * (*pddot_q) + p.R*cos(*pq) * (*pdddot_q);

  return status;
}

void release()
{
  gsl_odeiv2_driver_free (d);
}

//sgs

velocity c = 1.0;// = (double)(299792458 * 100);

void cset_c(long double _c)
{
    c = _c;
}

velocity cget_c()
{
    return c;
}

// расчет итерациями запаздывающего момента

static timespan timespan_Epsilon = 1.0e-15;// # погрешность
static distance distance_Epsilon = 1.0e-8;// # погрешность

static long double min_newton_step = 0.1;
static long double newton_step_multiplier = 0.9999;
static int max_steps = 50;
void cset_timespan_Epsilon(long double _eps)
{
    timespan_Epsilon = _eps;
}

void cset_distance_Epsilon(long double _eps)
{
    distance_Epsilon = _eps;
}

void cset_max_steps(int _max_steps)
{
    max_steps = _max_steps;
}

timespan cget_timespan_Epsilon()
{
    return timespan_Epsilon;
}

void cset_min_newton_step(long double min_step)
{
    min_newton_step = min_step;
}

void cset_newton_step_multiplier(long double multiplier)
{
   newton_step_multiplier = multiplier;
}




int no_retardation_test = 0;

void cset_no_retardation_test(int test)
{
    no_retardation_test = test;
}


/*timevalue newton_root_func(coordinate x, coordinate y, coordinate z,
                           timevalue t, timevalue t2,
                           Coordinate sx, Coordinate sy, Coordinate sz,
                           coordinate xc, coordinate yc, coordinate zc,
                           distance R, anglevelocity omega, angle alpha)
{
    long double f = (c*(t-t2)) - sqrtl(Sq(x-sx(t2, xc, yc, zc, R, omega, alpha)) +
                                      Sq(y-sy(t2, xc, yc, zc, R, omega, alpha)) +
                                      Sq(z-sz(t2, xc, yc, zc, R, omega, alpha)));
    return f;
}*/

void newton_root_derivative(coordinate x, coordinate y, coordinate z,
                            timevalue t, timevalue t2,
                            coordinate sx, coordinate sy, coordinate sz,
                            velocity vx, velocity vy, velocity vz,
                            double * cdt, double * d, double * f1,
                            double * dfdt)
{
    distance2 dd =
        Sq(x - sx) +
        Sq(y - sy) +
        Sq(z - sz);

    *cdt = c*(t-t2);
    *d = sqrtl(dd);
    *f1 = (*cdt) - (*d);

    double r_dot_v =
        (x - sx)*vx +
        (y - sy)*vy +
        (z - sz)*vz;

    *dfdt = -c + r_dot_v / (*d);
}

int NewtonIt(double step,
             coordinate x, coordinate y, coordinate z,
             timevalue t, timevalue t2,
             coordinate *psx, coordinate *psy, coordinate *psz,
             velocity *pvx, velocity *pvy, velocity *pvz,
             timevalue * res, double *f1)
{
    int ret = 0;
    double cdt, d, dfdt;
    newton_root_derivative(x, y, z,
                           t, t2,
                           *psx, *psy, *psz,
                           *pvx, *pvy, *pvz,
                           &cdt, &d, f1, &dfdt);

    if(fabsl(*f1) < distance_Epsilon)
    {
        //printf("NewtonIt fabsl(*f) %Le < Epsilon %Le\n", fabsl(*f), Epsilon);
        return +1;
    }

    if (dfdt == 0.0)
    {
        printf("NewtonIt error: derivative is zero\n");
        return -1;
    }

//    if (fabsl(dfdt) > 1.0)
//    {
//        printf("NewtonIt warning: derivative more than one\n");
//        printf("contractive mapping is not contractive\n");
//    }

    double delta      = (*f1)/dfdt;
    double step_delta = step*delta;
    *res                   = t2-step_delta;

    if (*res == t2)
    {
        printf("NewtonIt step_delta = %0.36e delta %0.36e step=%0.36e\n",
            step_delta, delta, step);
        printf("NewtonIt t2 %0.36e *res %0.36e f1 = %0.36e dfdt = %0.36e\n",
            t2, *res, (*f1), dfdt);
        ret = +2;
    }

    if (logging)
    {
        printf("t2 = %0.30e, f1 = %0.30e, dfdt=%0.30e, delta=%0.30e step %0.30e step_delta %0.30e res = %0.30e\n",
               t2, *f1, dfdt, delta, step, step_delta, *res);
    }

    return ret;
}

int find_newton_root(coordinate x, coordinate y, coordinate z, timevalue t, timevalue * pt2,
                           coordinate *psx, coordinate *psy, coordinate *psz,
                           velocity *pvx, velocity *pvy, velocity *pvz)
{
    double step = 1.0;
    double t1;
    double f;
    int ret;

    double Momenta, q;

    double dot_q;   // dphi/dt
    double ddot_q;  // d2phi/dt2
    double dddot_q; // d3phi/dt3

    double wx;
    double wy;

    double dot_wx;
    double dot_wy;

    for (;;)
    {
        if (logging) printf("t2=%0.30e\t", *pt2);
        t1 = *pt2;

        apply(*pt2, NULL, &Momenta, &q,
              &dot_q, &ddot_q, &dddot_q,
              psx, psy,
              pvx, pvy,
              &wx, &wy,
              &dot_wx, &dot_wy
        );

        ret = NewtonIt(step,
                      x, y, z,
                      t, *pt2,
                      psx, psy, psz,
                      pvx, pvy, pvz, pt2, &f);
        if (0 != ret)
        {
            return ret;
        }

        if (*pt2 == t1){
            printf("NewtonIt *pt2 %e == t1 %e step=%e\n", *pt2, t1, step);
            break;
        }

        if (step > min_newton_step)
            step *= newton_step_multiplier;
    }

    return 0;
}


double period_Epsilon = 1.0e-16;

int find_period_by_newton_root(timevalue * pt2, double *pf)
{
    double step = 1.0;
    double t1;
    int ret;

    double Momenta, q;

    double dot_q;   // dphi/dt
    double ddot_q;  // d2phi/dt2
    double dddot_q; // d3phi/dt3

    coordinate sx, sy;
    velocity vx, vy;

    double wx;
    double wy;

    double dot_wx;
    double dot_wy;

    for (;;)
    {
        if (logging) printf("t2=%0.30e\t", *pt2);
        t1 = *pt2;

        apply(*pt2, NULL, &Momenta, &q,
              &dot_q, &ddot_q, &dddot_q,
              &sx, &sy,
              &vx, &vy,
              &wx, &wy,
              &dot_wx, &dot_wy
        );

        /*ret = NewtonIt(step,
                      x, y, z,
                      t, *pt2,
                      psx, psy, psz,
                      pvx, pvy, pvz, pt2, &f);
*/

        {
            int ret = 0;
            double res, dfdt;
            /*newton_root_derivative(x, y, z,
                                t, *pt2,
                                *psx, *psy, *psz,
                                *pvx, *pvy, *pvz,
                                &cdt, &d, &f, &dfdt);*/

            *pf = 2*M_PI - q;
            dfdt = dot_q;

            *pf = sy;
            dfdt = vy;

            if(fabsl(*pf) < period_Epsilon)
            {
                if (logging) printf("NewtonIt fabsl(*f) %Le < Epsilon %e\n",
                    fabsl(*pf), period_Epsilon);
                return +1;
            }

            if (dfdt == 0.0)
            {
                printf("NewtonIt error: derivative is zero\n");
                return -1;
            }

        //    if (fabsl(dfdt) > 1.0)
        //    {
        //        printf("NewtonIt warning: derivative more than one\n");
        //        printf("contractive mapping is not contractive\n");
        //    }

            double delta      = (*pf)/dfdt;
            double step_delta = step*delta;
            res                   = *pt2-step_delta;

            if (res == *pt2)
            {
                printf("NewtonIt step_delta = %0.36e delta %0.36e step=%0.36e\n",
                    step_delta, delta, step);
                printf("NewtonIt t2 %0.36e *res %0.36e f1 = %0.36e dfdt = %0.36e\n",
                    *pt2, res, *pf, dfdt);
                ret = +2;
            }

            if (logging)
            {
                printf("t2 = %0.30e, f1 = %0.30e, dfdt=%0.30e, delta=%0.30e step %0.30e step_delta %0.30e res = %0.30e\n",
                    *pt2, *pf, dfdt, delta, step, step_delta, res);
            }

            *pt2 = res;
        }
        if (*pt2 == t1){
            printf("NewtonIt *pt2 %e == t1 %e step=%e\n", *pt2, t1, step);
            break;
        }

        if (step > min_newton_step)
            step *= newton_step_multiplier;
    }

    return ret;
}

void calc_pendulum_period(double * p_init_T, double max_ti, double dti)
{
    double f;
    int ret;

    double Momenta, q;

    double dot_q;   // dphi/dt
    double ddot_q;  // d2phi/dt2
    double dddot_q; // d3phi/dt3

    coordinate sx, sy;
    velocity vx, vy;

    double wx;
    double wy;

    double dot_wx;
    double dot_wy;

    double pre_ti = 0;
    double pre_sy = 0;
    int start = 1;
    int null_index = 0;

    for (double ti = 0; ti < max_ti; ti += dti)
    {
        apply(ti, NULL, &Momenta, &q,
              &dot_q, &ddot_q, &dddot_q,
              &sx, &sy,
              &vx, &vy,
              &wx, &wy,
              &dot_wx, &dot_wy
        );
        //printf("ti=%f sy=%f\n", ti, sy);

        if (pre_sy * sy < 0.0) {
            if (1 == null_index)
            {
                *p_init_T = (pre_ti+ti)/2;
                printf("init_T=%f\n", *p_init_T);
            }

            printf("(pre_sy=%f,sy=%f)\n", pre_sy,sy);
            printf("(pre_ti=%f,ti=%f)\n", pre_ti,ti);

            null_index += 1;
        }
        if (pre_sy * sy != 0 || start)
        {
            pre_ti = ti;
            pre_sy = sy;
            if (sy != 0.0){
                start = 0;
            }
        }
    }
}

int tlag(coordinate x, coordinate y, coordinate z, timevalue t,
         coordinate *psx, coordinate *psy, coordinate *psz,
         velocity *pvx, velocity *pvy, velocity *pvz,
         acceleration *pwx, acceleration *pwy, acceleration *pwz,
         double * pdot_wx, double * pdot_wy, double * pdot_wz,
         timevalue * pt2, coordinate * rlagerror)
{
    double Momenta, q;

    double dot_q;   // dphi/dt
    double ddot_q;  // d2phi/dt2
    double dddot_q; // d3phi/dt3

    *psz = 0.0;
    *pvz = 0.0;
    *pwz = 0.0;
    *pdot_wz = 0.0;

    if (no_retardation_test){
        *pt2 = t;
        return 0;
    }

    timevalue t1 = t;
    *pt2 = t - 2*timespan_Epsilon;
    double dd, d, cdt;

    int n = 0;

    while (1) {
        t1 = *pt2;

        apply(*pt2, NULL, &Momenta, &q,
              &dot_q, &ddot_q, &dddot_q,
              psx, psy,
              pvx, pvy,
              pwx, pwy,
              pdot_wx, pdot_wy
        );

        dd =
            Sq(x - *psx) +
            Sq(y - *psy) +
            Sq(z - *psz);

        d = sqrtl(dd);
        *pt2 = t - d / c;
        cdt = c*(t-*pt2);

        //printf("d=%0.30Le, c(t-t2)=%0.30Le\n", d, c*(t-t2));
        //printf("d-c(t-t2)=%0.30Le\n", d-c*(t-t2));

        if(d/c < timespan_Epsilon)
            break;

        if(fabsl(cdt - d) < distance_Epsilon)
            break;

        if (*pt2 == t1)
            break;
    }

    int ret = find_newton_root(x, y, z, t, pt2,
                     psx, psy, psz,
                     pvx, pvy, pvz);

    // calc tlagtest

    dd =
        Sq(x - *psx) +
        Sq(y - *psy) +
        Sq(z - *psz);
    d = sqrtl(dd);

    cdt = c*(t-(*pt2));
    *rlagerror = cdt - d;
    //printf("tlag test: d=%Lf cdt=%Lf (cdt - d)=%Le fabsl(cdt - d)=%Le\n", d, cdt, (cdt - d), fabsl(cdt - d));

    if (fabsl(*rlagerror) > distance_Epsilon)
    {
        printf("tlag test error: x=%f y=%f z=%f t=%f t2=%f\n", x, y, z, t, *pt2);
        printf("tlag test error: d=%f cdt=%f (cdt - d)=%e\n", d, cdt, (*rlagerror));
        return -1;
    }

    return 0;
}



/*
Здесь нужно отметить во избежание путаницы, что радиус Лиенара Вихерта $k$ в формуле для электрического и магнитного полей и $k$ в программе отличаются друг от друга тем, что в программе $k$ нормирован на единицу, а в формуле нет. При переходе от формул к програмным кодам по сути произведено преобразование $k \rightarrow k\cdot r$
*/

void calc_k(coordinate x, coordinate y, coordinate z, timevalue t,
            coordinate sx, coordinate sy, coordinate sz,
            velocity vx, velocity vy, velocity vz,
            timevalue t2,
            double * k, distance * r, distance * nx, distance * ny, distance * nz)
{
    (*r) = c * (t - t2);

    //printf("(*r) = %Lf t = %Lf t2 = %Lf (t - t2)=%Lf\n", (*r), t, t2, (t - t2));

    (*nx) = (x - sx)/(*r);
    (*ny) = (y - sy)/(*r);
    (*nz) = (z - sz)/(*r);

    //printf("(*nx) = %Lf (*ny) = %Lf (*nz) = %Lf\n", (*nx), (*ny), (*nz));

    (*k) = 1.0 - ((*nx) * vx +
                  (*ny) * vy +
                  (*nz) * vz) / c;

    //printf("(*r) = %Le (*k) = %Le\n", (*r), (*k));
}

// отношение радиуса Лиенара Вихерта к длине радиус-вектора
int klw(coordinate x, coordinate y, coordinate z, timevalue t, timevalue * pt2,
        double *pk, coordinate * rlagerror)
{
    distance r;
    distance nx;
    distance ny;
    distance nz;

    coordinate sx, sy, sz;
    velocity vx, vy, vz;
    acceleration wx, wy, wz;
    double dot_wx, dot_wy, dot_wz;

    // расчет итерациями запаздывающего момента
    if (0 == tlag(x, y, z, t, &sx, &sy, &sz, &vx, &vy, &vz,
                  &wx, &wy, &wz,
                  &dot_wx, &dot_wy, &dot_wz,
                  pt2, rlagerror)) {

        calc_k(x, y, z, t, sx, sy, sz, vx, vy, vz, *pt2, pk, &r, &nx, &ny, &nz);
        //printf("klw (*r) = %e (*k) = %e t2 = %e\n", r, *pk, t2);
        return 0;
    }
    return -1;
}

// Радиус Лиенара Вихерта
int Rlw(coordinate x, coordinate y, coordinate z, timevalue t,
        double *pRlw, coordinate * rlagerror)
{
    double k;
    distance r;
    distance nx;
    distance ny;
    distance nz;

    coordinate sx, sy, sz;
    velocity vx, vy, vz;
    acceleration wx, wy, wz;
    double dot_wx, dot_wy, dot_wz;

    timevalue t2;
    // расчет итерациями запаздывающего момента
    if (0 == tlag(x, y, z, t, &sx, &sy, &sz, &vx, &vy, &vz,
                  &wx, &wy, &wz,
                  &dot_wx, &dot_wy, &dot_wz,
                  &t2, rlagerror)) {
        calc_k(x, y, z, t, sx, sy, sz, vx, vy, vz, t2, &k, &r, &nx, &ny, &nz);
        //printf("Rlw (*r) = %Le (*k) = %Le t2 = %Le\n", r, k, t2);

        *pRlw = k*r;
        return 0;
    }
    return -1;
}

// phi_lw - скалярный потенциал Лиенара Вихерта
int philw(coordinate x, coordinate y, coordinate z, timevalue t,
          charge q,
          double *pphi, coordinate * rlagerror)
{
    double k;
    distance r;
    distance nx;
    distance ny;
    distance nz;

    coordinate sx, sy, sz;
    velocity vx, vy, vz;
    acceleration wx, wy, wz;
    double dot_wx, dot_wy, dot_wz;

    timevalue t2;
    // расчет итерациями запаздывающего момента
    if (0 == tlag(x, y, z, t, &sx, &sy, &sz, &vx, &vy, &vz,
                  &wx, &wy, &wz,
                  &dot_wx, &dot_wy, &dot_wz,
                  &t2, rlagerror)) {
        calc_k(x, y, z, t, sx, sy, sz, vx, vy, vz, t2, &k, &r, &nx, &ny, &nz);

        *pphi = q/(k*r);

        return 0;
    }

    return -1;
}


// A_lw - векторный потенциал Лиенара Вихерта
int Alw(coordinate x, coordinate y, coordinate z, timevalue t,
        charge q,
        field * A_x, field * A_y, field * A_z, coordinate * rlagerror)
{
    double k;
    distance r;
    distance nx;
    distance ny;
    distance nz;

    coordinate sx, sy, sz;
    velocity vx, vy, vz;
    acceleration wx, wy, wz;
    double dot_wx, dot_wy, dot_wz;

    timevalue t2;
    // расчет итерациями запаздывающего момента
    if (0 == tlag(x, y, z, t, &sx, &sy, &sz, &vx, &vy, &vz,
                  &wx, &wy, &wz,
                  &dot_wx, &dot_wy, &dot_wz,
                  &t2, rlagerror)) {
        calc_k(x, y, z, t, sx, sy, sz, vx, vy, vz, t2, &k, &r, &nx, &ny, &nz);

        (*A_x) = q*vx/(k*r);
        (*A_y) = q*vy/(k*r);
        (*A_z) = q*vz/(k*r);

        return 0;
    }
    return -1;
}


int electr_magnet(coordinate x, coordinate y, coordinate z, timevalue t,
                   charge q,
                   field * E_x, field * E_y, field * E_z,
                   field * B_x, field * B_y, field * B_z,
                   coordinate * rlagerror)
{
    double k;
    distance r;
    distance nx;
    distance ny;
    distance nz;

    coordinate sx, sy, sz;
    velocity vx, vy, vz;
    acceleration wx, wy, wz;
    double dot_wx, dot_wy, dot_wz;

    timevalue t2;
    // расчет итерациями запаздывающего момента
    if (0 == tlag(x, y, z, t, &sx, &sy, &sz, &vx, &vy, &vz,
                  &wx, &wy, &wz,
                  &dot_wx, &dot_wy, &dot_wz,
                  &t2, rlagerror)) {
        calc_k(x, y, z, t, sx, sy, sz, vx, vy, vz, t2, &k, &r, &nx, &ny, &nz);


        double v2_c2 = (Sq(vx) + Sq(vy) + Sq(vz)) / (c*c);
        double ra_c2 = r * (nx*wx + ny*wy + nz*wz) / (c*c);


        (*E_x) = q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nx - vx/c)/(r*r) - (k/r)*wx/(c*c));
        (*E_y) = q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(ny - vy/c)/(r*r) - (k/r)*wy/(c*c));
        (*E_z) = q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nz - vz/c)/(r*r) - (k/r)*wz/(c*c));

        (*B_x) = -q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(ny*vz - nz*vy)/(r*r)/c + (ny*wz - nz*wy)*(k/r)/(c*c));
        (*B_y) = -q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nz*vx - nx*vz)/(r*r)/c + (nz*wx - nx*wz)*(k/r)/(c*c));
        (*B_z) = -q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nx*vy - ny*vx)/(r*r)/c + (nx*wy - ny*wx)*(k/r)/(c*c));

        return 0;
    }
    return -1;
}


int electr_magnet_ex(coordinate x, coordinate y, coordinate z, timevalue t,
                   charge q,
                   field * E_x, field * E_y, field * E_z,
                   field * B_x, field * B_y, field * B_z,
                   field * A_x, field * A_y, field * A_z,
                   coordinate * rlagerror)
{
    double k;
    distance r;
    distance nx;
    distance ny;
    distance nz;

    coordinate sx, sy, sz;
    velocity vx, vy, vz;
    acceleration wx, wy, wz;
    double dot_wx, dot_wy, dot_wz;

    timevalue t2;
    // расчет итерациями запаздывающего момента
    if (0 == tlag(x, y, z, t, &sx, &sy, &sz, &vx, &vy, &vz,
                  &wx, &wy, &wz,
                  &dot_wx, &dot_wy, &dot_wz,
                  &t2, rlagerror)) {

        calc_k(x, y, z, t, sx, sy, sz, vx, vy, vz, t2, &k, &r, &nx, &ny, &nz);


        double v2_c2 = (Sq(vx) + Sq(vy) + Sq(vz)) / (c*c);
        double ra_c2 = r * (nx*wx + ny*wy + nz*wz) / (c*c);

        (*E_x) = q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nx - vx/c)/(r*r) - (k/r)*wx/(c*c));
        (*E_y) = q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(ny - vy/c)/(r*r) - (k/r)*wy/(c*c));
        (*E_z) = q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nz - vz/c)/(r*r) - (k/r)*wz/(c*c));

        (*B_x) = -q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(ny*vz - nz*vy)/(r*r)/c + (ny*wz - nz*wy)*(k/r)/(c*c));
        (*B_y) = -q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nz*vx - nx*vz)/(r*r)/c + (nz*wx - nx*wz)*(k/r)/(c*c));
        (*B_z) = -q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nx*vy - ny*vx)/(r*r)/c + (nx*wy - ny*wx)*(k/r)/(c*c));

        (*A_x) = q*vx/(k*r);
        (*A_y) = q*vy/(k*r);
        (*A_z) = q*vz/(k*r);

        return 0;
    }
    return -1;
}
