#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include "lw.h"
#include "lw_drobyshev_task.h"

#include "stdlib.h"

#define _USE_MATH_DEFINES
#include <math.h>
#define Pi M_PI

#define Sq(x) ((x)*(x))

//sgs

velocity c = 1.0; //(double)(299792458 * 100);

void set_c(double _c)
{
    c = _c;
}

velocity cget_c()
{
    return c;
}

// расчет итерациями запаздывающего момента

static timespan timespan_Epsilon = 1.0e-16;// # погрешность
static distance distance_Epsilon = 1.0e-16;// # погрешность

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


int logging = 0;

int no_retardation_test = 0;

void cset_no_retardation_test(int test)
{
    no_retardation_test = test;
}

typedef timevalue (*Tlag)(coordinate x, coordinate y, coordinate z, timevalue t,
                          Coordinate sx, Coordinate sy, Coordinate sz,
                          Velocity vx, Velocity vy, Velocity vz,
                          coordinate xc, coordinate yc, coordinate zc,
                          distance R, anglevelocity omega, angle alpha);


timevalue newton_root_func2(coordinate x, coordinate y, coordinate z,
                           timevalue t, timevalue t2,
                           Coordinate sx, Coordinate sy, Coordinate sz,
                           coordinate xc, coordinate yc, coordinate zc,
                           distance R, anglevelocity omega, angle alpha)
{
    long double f = Sq(c*(t-t2))
             - Sq(x-sx(t2))
             - Sq(y-sy(t2))
             - Sq(z-sz(t2));
    return f;
}

timevalue newton_root_derivative2(coordinate x, coordinate y, coordinate z,
                                 timevalue t, timevalue t2,
                                 Coordinate sx, Coordinate sy, Coordinate sz,
                                 Velocity vx, Velocity vy, Velocity vz,
                                 coordinate xc, coordinate yc, coordinate zc,
                                 distance R, anglevelocity omega, angle alpha)
{
    long double dfdt2 = - 2 * c*c * (t - t2)
                 + 2*(x - sx(t2))*vx(t2)
                 + 2*(y - sy(t2))*vy(t2)
                 + 2*(z - sz(t2))*vz(t2);
    return dfdt2;
}



timevalue tlag_test(coordinate x, coordinate y, coordinate z, timevalue t1, timevalue t2,
                    Coordinate sx, Coordinate sy, Coordinate sz,
                    coordinate xc, coordinate yc, coordinate zc,
                    distance R, anglevelocity omega, angle alpha, long double * d, long double * cdt)
{
    long double dd =
        Sq(x - sx(t2)) +
        Sq(y - sy(t2)) +
        Sq(z - sz(t2));
    *d = sqrtl(dd);
    *cdt = c*(t1-t2);
    return (*cdt) - (*d);
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
                            Coordinate sx, Coordinate sy, Coordinate sz,
                            Velocity vx, Velocity vy, Velocity vz,
                            long double * cdt, long double * d, long double * f1,
                            long double * dfdt)
{
    long double dd =
        Sq(x - sx(t2)) +
        Sq(y - sy(t2)) +
        Sq(z - sz(t2));

    *cdt = c*(t-t2);
    *d = sqrtl(dd);
    *f1 = (*cdt) - (*d);

    long double r_dot_v =
        (x - sx(t2))*vx(t2) +
        (y - sy(t2))*vy(t2) +
        (z - sz(t2))*vz(t2);

    *dfdt = -c + r_dot_v / (*d);
}

int NewtonIt(long double step,
             coordinate x, coordinate y, coordinate z,
             timevalue t, timevalue t2,
             Coordinate sx, Coordinate sy, Coordinate sz,
             Velocity vx, Velocity vy, Velocity vz,
             timevalue * res, long double *f1)
{
    int ret = 0;
    long double cdt, d, dfdt;
    newton_root_derivative(x, y, z,
                           t, t2,
                           sx, sy, sz,
                           vx, vy, vz,
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

    long double delta      = (*f1)/dfdt;
    long double step_delta = step*delta;
    *res                   = t2-step_delta;

    if (*res == t2)
    {
        printf("NewtonIt step_delta = %0.36Le delta %0.36Le step=%0.36Le\n",
            step_delta, delta, step);
        printf("NewtonIt t2 %0.36Le *res %0.36Le f1 = %0.36Le dfdt = %0.36Le\n",
            t2, *res, (*f1), dfdt);
        ret = +2;
    }

    if (logging)
    {
        printf("t2 = %0.30Le, f1 = %0.30Le, dfdt=%0.30Le, delta=%0.30Le step %0.30Le step_delta %0.30Le res = %0.30Le\n",
               t2, *f1, dfdt, delta, step, step_delta, *res);
    }

    return ret;
}

int find_newton_root(coordinate x, coordinate y, coordinate z, timevalue t, timevalue * pt2,
                           Coordinate sx, Coordinate sy, Coordinate sz,
                           Velocity vx, Velocity vy, Velocity vz)
{
    long double step = 1.0;
    long double t1;
    long double f;
    int ret;

    for (;;)
    {
        if (logging) printf("t2=%0.30Le\t", *pt2);
        t1 = *pt2;
        ret = NewtonIt(step,
                      x, y, z,
                      t, *pt2,
                      sx, sy, sz,
                      vx, vy, vz,
                      pt2, &f);
        if (0 != ret)
        {
            return ret;
        }

        if (*pt2 == t1){
            printf("NewtonIt *pt2 %Le == t1 %Le step=%Le\n", *pt2, t1, step);
            break;
        }

        if (step > min_newton_step)
            step *= newton_step_multiplier;
    }

    return 0;
}

int tlag(coordinate x, coordinate y, coordinate z, timevalue t,
               Coordinate sx, Coordinate sy, Coordinate sz,
               Velocity vx, Velocity vy, Velocity vz,
               timevalue * pt2, coordinate * rlagerror)
{
    if (no_retardation_test){
        *pt2 = t;
        return 0;
    }

    timevalue t1 = t;
    *pt2 = t - 2*timespan_Epsilon;
    long double dd, d, cdt;

    int n = 0;

    while (1) {
        t1 = *pt2;
        dd =
            Sq(x - sx(*pt2)) +
            Sq(y - sy(*pt2)) +
            Sq(z - sz(*pt2));
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
                     sx, sy, sz,
                     vx, vy, vz);

    // calc tlagtest

    dd =
        Sq(x - sx(*pt2)) +
        Sq(y - sy(*pt2)) +
        Sq(z - sz(*pt2));
    d = sqrtl(dd);

    cdt = c*(t-(*pt2));
    *rlagerror = cdt - d;
    //printf("tlag test: d=%Lf cdt=%Lf (cdt - d)=%Le fabsl(cdt - d)=%Le\n", d, cdt, (cdt - d), fabsl(cdt - d));

    if (fabsl(*rlagerror) > distance_Epsilon)
    {
        printf("tlag test error: x=%Lf y=%Lf z=%Lf t=%Lf t2=%Lf\n", x, y, z, t, *pt2);
        printf("tlag test error: d=%Lf cdt=%Lf (cdt - d)=%Le\n", d, cdt, (*rlagerror));
        return -1;
    }

    return 0;
}

/*
Здесь нужно отметить во избежание путаницы, что радиус Лиенара Вихерта $k$ в формуле для электрического и магнитного полей и $k$ в программе отличаются друг от друга тем, что в программе $k$ нормирован на единицу, а в формуле нет. При переходе от формул к програмным кодам по сути произведено преобразование $k \rightarrow k\cdot r$
*/

void calc_k(coordinate x, coordinate y, coordinate z, timevalue t,
            Coordinate sx, Coordinate sy, Coordinate sz,
            Velocity vx, Velocity vy, Velocity vz,
            timevalue t2,
            long double * k, distance * r, distance * nx, distance * ny, distance * nz)
{
    if (no_retardation_test) {
        (*r) = sqrtl((x - sx(t))*(x - sx(t)) +
                     (y - sy(t))*(y - sy(t)) +
                     (z - sz(t))*(z - sz(t)));
    }
    else {
        (*r) = c * (t - t2);
    }

    (*nx) = (x - sx(t2))/(*r);
    (*ny) = (y - sy(t2))/(*r);
    (*nz) = (z - sz(t2))/(*r);

    //printf("(*nx) = %f (*ny) = %f (*nz) = %f\n", (*nx), (*ny), (*nz));

    if (no_retardation_test) {
        (*k) = 1.0;
    }
    else {
        (*k) = 1.0 - ((*nx)*vx(t2) + (*ny) * vy(t2) + (*nz) * vz(t2)) / c;
    }

    //printf("(*r) = %e (*k) = %e\n", (*r), (*k));
}

// отношение радиуса Лиенара Вихерта к радиусу
int klw(coordinate x, coordinate y, coordinate z, timevalue t,
        Coordinate sx, Coordinate sy, Coordinate sz,
        Velocity vx, Velocity vy, Velocity vz,
        long double *pk, coordinate * rlagerror)
{
    distance r;
    distance nx;
    distance ny;
    distance nz;

    long double t2;
    // расчет итерациями запаздывающего момента
    if (0 == tlag(x, y, z, t, sx, sy, sz, vx, vy, vz,
                      &t2, rlagerror)) {

        calc_k(x, y, z, t, sx, sy, sz, vx, vy, vz, t2, pk, &r, &nx, &ny, &nz);
        //printf("klw (*r) = %Le (*k) = %Le t2 = %Le\n", r, k, t2);
        return 0;
    }
    return -1;
}

// Радиус Лиенара Вихерта
int Rlw(coordinate x, coordinate y, coordinate z, timevalue t,
        Coordinate sx, Coordinate sy, Coordinate sz,
        Velocity vx, Velocity vy, Velocity vz,
        long double *pRlw, coordinate * rlagerror)
{
    long double k;
    distance r;
    distance nx;
    distance ny;
    distance nz;

    long double t2;
    // расчет итерациями запаздывающего момента
    if (0 == tlag(x, y, z, t, sx, sy, sz, vx, vy, vz,
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
          Coordinate sx, Coordinate sy, Coordinate sz,
          Velocity vx, Velocity vy, Velocity vz,
          charge q,
          long double *pphi, coordinate * rlagerror)
{
    long double k;
    distance r;
    distance nx;
    distance ny;
    distance nz;

    long double t2;
    // расчет итерациями запаздывающего момента
    if (0 == tlag(x, y, z, t, sx, sy, sz, vx, vy, vz,
                     &t2, rlagerror))
    {
        calc_k(x, y, z, t, sx, sy, sz, vx, vy, vz, t2, &k, &r, &nx, &ny, &nz);

        *pphi = q/(k*r);

        return 0;
    }

    return -1;
}

// A_lw - векторный потенциал Лиенара Вихерта
int Alw(coordinate x, coordinate y, coordinate z, timevalue t,
        Coordinate sx, Coordinate sy, Coordinate sz,
        Velocity vx, Velocity vy, Velocity vz,
        charge q,
        field * A_x, field * A_y, field * A_z,
        coordinate * rlagerror
       )
{
    long double k;
    distance r;
    distance nx;
    distance ny;
    distance nz;

    long double t2;
    if (0 == tlag(x, y, z, t, sx, sy, sz, vx, vy, vz,
                      &t2, rlagerror)){ // расчет итерациями запаздывающего момента
        calc_k(x, y, z, t, sx, sy, sz, vx, vy, vz, t2, &k, &r, &nx, &ny, &nz);

        (*A_x) = q*vx(t2)/(k*r);
        (*A_y) = q*vy(t2)/(k*r);
        (*A_z) = q*vz(t2)/(k*r);
        return 0;
    }
    return -1;
}

int electr_magnet(coordinate x, coordinate y, coordinate z, timevalue t,
                  Coordinate sx, Coordinate sy, Coordinate sz,
                  Velocity vx, Velocity vy, Velocity vz,
                  Acceleration wx, Acceleration wy, Acceleration wz,
                  charge q,
                  field * E_x, field * E_y, field * E_z,
                  field * B_x, field * B_y, field * B_z,
                  coordinate * rlagerror)
{
    long double k;
    distance r;
    distance nx;
    distance ny;
    distance nz;

    // расчет итерациями запаздывающего момента
    long double t2;
    if (0 == tlag(x, y, z, t, sx, sy, sz, vx, vy, vz,
                  &t2, rlagerror)) {

        calc_k(x, y, z, t, sx, sy, sz, vx, vy, vz, t2, &k, &r, &nx, &ny, &nz);

        long double v_x = vx(t2);
        long double v_y = vy(t2);
        long double v_z = vz(t2);

        long double w_x = wx(t2);
        long double w_y = wy(t2);
        long double w_z = wz(t2);

        long double v2_c2 = (Sq(v_x) + Sq(v_y) + Sq(v_z)) / (c*c);
        long double ra_c2 = r * (nx*w_x + ny*w_y + nz*w_z) / (c*c);

        (*E_x) = q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nx - v_x/c)/(r*r) - (k/r)*w_x/(c*c));
        (*E_y) = q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(ny - v_y/c)/(r*r) - (k/r)*w_y/(c*c));
        (*E_z) = q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nz - v_z/c)/(r*r) - (k/r)*w_z/(c*c));

        (*B_x) = -q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(ny*v_z - nz*v_y)/(r*r)/c + (ny*w_z - nz*w_y)*(k/r)/(c*c));
        (*B_y) = -q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nz*v_x - nx*v_z)/(r*r)/c + (nz*w_x - nx*w_z)*(k/r)/(c*c));
        (*B_z) = -q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nx*v_y - ny*v_x)/(r*r)/c + (nx*w_y - ny*w_x)*(k/r)/(c*c));

        return 0;
    }
    return -1;
}

int electr_magnet_ex(coordinate x, coordinate y, coordinate z, timevalue t,
                   Coordinate sx, Coordinate sy, Coordinate sz,
                   Velocity vx, Velocity vy, Velocity vz,
                   Acceleration wx, Acceleration wy, Acceleration wz,
                   charge q,
                   field * E_x, field * E_y, field * E_z,
                   field * B_x, field * B_y, field * B_z,
                   field * A_x, field * A_y, field * A_z,
                   coordinate * rlagerror)
{
    long double k;
    distance r;
    distance nx;
    distance ny;
    distance nz;

    // расчет итерациями запаздывающего момента
    long double t2;
    if (0 == tlag(x, y, z, t, sx, sy, sz, vx, vy, vz,
                  &t2, rlagerror)) {

        calc_k(x, y, z, t, sx, sy, sz, vx, vy, vz, t2, &k, &r, &nx, &ny, &nz);

        long double v_x = vx(t2);
        long double v_y = vy(t2);
        long double v_z = vz(t2);

        long double w_x = wx(t2);
        long double w_y = wy(t2);
        long double w_z = wz(t2);

        long double v2_c2 = (Sq(v_x) + Sq(v_y) + Sq(v_z)) / (c*c);
        long double ra_c2 = r * (nx*w_x + ny*w_y + nz*w_z) / (c*c);

        (*E_x) = q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nx - v_x/c)/(r*r) - (k/r)*w_x/(c*c));
        (*E_y) = q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(ny - v_y/c)/(r*r) - (k/r)*w_y/(c*c));
        (*E_z) = q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nz - v_z/c)/(r*r) - (k/r)*w_z/(c*c));

        (*B_x) = -q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(ny*v_z - nz*v_y)/(r*r)/c + (ny*w_z - nz*w_y)*(k/r)/(c*c));
        (*B_y) = -q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nz*v_x - nx*v_z)/(r*r)/c + (nz*w_x - nx*w_z)*(k/r)/(c*c));
        (*B_z) = -q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nx*v_y - ny*v_x)/(r*r)/c + (nx*w_y - ny*w_x)*(k/r)/(c*c));

        (*A_x) = q*v_x/(k*r);
        (*A_y) = q*v_y/(k*r);
        (*A_z) = q*v_z/(k*r);

        return 0;
    }
    return -1;
}
