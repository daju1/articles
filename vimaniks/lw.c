#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include "lw.h"

#include "stdlib.h"

#define _USE_MATH_DEFINES
#include <math.h>
#define Pi M_PI

#define Sq(x) ((x)*(x))

//sgs 

static velocity c = (double)(299792458 * 100);

void set_c(double _c)
{
    c = _c;
}

velocity cget_c()
{
    return c;
}

typedef coordinate (*Coordinate)(timevalue t_zap);
typedef velocity (*Velocity)(timevalue t_zap);
typedef acceleration (*Acceleration)(timevalue t_zap);

// расчет итерациями запаздывающего момента

static timespan Epsilon = 1.0e-16;// # погрешность

void set_timespan_Epsilon(double _eps)
{
    Epsilon = _eps;
}

timespan cget_timespan_Epsilon()
{
    return Epsilon;
}

_Bool no_retardation_test = 0;

typedef timevalue (*Tlag)(coordinate x, coordinate y, coordinate z, timevalue t,
                          Coordinate sx, Coordinate sy, Coordinate sz);

    
timevalue tlag(coordinate x, coordinate y, coordinate z, timevalue t,
               Coordinate sx, Coordinate sy, Coordinate sz)
{
    if (no_retardation_test){
        return t;
    }

    timevalue t1 = t;
    timevalue t2 = t - 2*Epsilon;
    
    int n = 0;

    while (fabs(t1 - t2) > Epsilon) {
        t1 = t2;
        double dd = (x - sx(t1))*(x - sx(t1)) +
                    (y - sy(t1))*(y - sy(t1)) +
                    (z - sz(t1))*(z - sz(t1));
        double d = sqrt(dd);
        t2 = t - d / c;
        if (++n > 100)
            break;
    }

    return t2;
}

/*
Здесь нужно отметить во избежание путаницы, что радиус Лиенара Вихерта $k$ в формуле для электрического и магнитного полей и $k$ в программе отличаются друг от друга тем, что в программе $k$ нормирован на единицу, а в формуле нет. При переходе от формул к програмным кодам по сути произведено преобразование $k \rightarrow k\cdot r$
*/

void calc_k(coordinate x, coordinate y, coordinate z, timevalue t,
            Coordinate sx, Coordinate sy, Coordinate sz,
            Velocity vx, Velocity vy, Velocity vz,
            timevalue t2,
            double * k, distance * r, distance * nx, distance * ny, distance * nz)
{
    
    if (no_retardation_test) {
        (*r) = sqrt((x - sx(t))*(x - sx(t)) +
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
double klw(coordinate x, coordinate y, coordinate z, timevalue t,
                Coordinate sx, Coordinate sy, Coordinate sz,
                Velocity vx, Velocity vy, Velocity vz,
                Tlag t_lag)
{
    double k;
    distance r;
    distance nx;
    distance ny;
    distance nz;

    double t2 = t_lag(x, y, z, t, sx, sy, sz); // расчет итерациями запаздывающего момента
    
    calc_k(x, y, z, t, sx, sy, sz, vx, vy, vz, t2, &k, &r, &nx, &ny, &nz);
    return k;
}

// Радиус Лиенара Вихерта
double Rlw(coordinate x, coordinate y, coordinate z, timevalue t,
                Coordinate sx, Coordinate sy, Coordinate sz,
                Velocity vx, Velocity vy, Velocity vz,
                Tlag t_lag)
{
    double k;
    distance r;
    distance nx;
    distance ny;
    distance nz;

    double t2 = t_lag(x, y, z, t, sx, sy, sz); // расчет итерациями запаздывающего момента
    calc_k(x, y, z, t, sx, sy, sz, vx, vy, vz, t2, &k, &r, &nx, &ny, &nz);

    return k*r;
}

// phi_lw - скалярный потенциал Лиенара Вихерта
double philw(coordinate x, coordinate y, coordinate z, timevalue t,
           Coordinate sx, Coordinate sy, Coordinate sz,
           Velocity vx, Velocity vy, Velocity vz,
           charge q, Tlag t_lag)
{
    double k;
    distance r;
    distance nx;
    distance ny;
    distance nz;

    double t2 = t_lag(x, y, z, t, sx, sy, sz); // расчет итерациями запаздывающего момента
    calc_k(x, y, z, t, sx, sy, sz, vx, vy, vz, t2, &k, &r, &nx, &ny, &nz);

    return q/(k*r);
}


// A_lw - векторный потенциал Лиенара Вихерта
void Alw(coordinate x, coordinate y, coordinate z, timevalue t,
       Coordinate sx, Coordinate sy, Coordinate sz,
       Velocity vx, Velocity vy, Velocity vz,
       charge q, Tlag t_lag,
       field * A_x, field * A_y, field * A_z
       )
{
    double k;
    distance r;
    distance nx;
    distance ny;
    distance nz;

    double t2 = t_lag(x, y, z, t, sx, sy, sz); // расчет итерациями запаздывающего момента
    calc_k(x, y, z, t, sx, sy, sz, vx, vy, vz, t2, &k, &r, &nx, &ny, &nz);

    (*A_x) = q*vx(t2)/(k*r);
    (*A_y) = q*vy(t2)/(k*r);
    (*A_z) = q*vz(t2)/(k*r);
}


void electr_magnet(coordinate x, coordinate y, coordinate z, timevalue t,
                Coordinate sx, Coordinate sy, Coordinate sz,
                Velocity vx, Velocity vy, Velocity vz,
                Acceleration wx, Acceleration wy, Acceleration wz,
                charge q, Tlag t_lag, 
                field * E_x, field * E_y, field * E_z, field * B_x, field * B_y, field * B_z)
{
    double k;
    distance r;
    distance nx;
    distance ny;
    distance nz;

    double t2 = t_lag(x, y, z, t, sx, sy, sz); // расчет итерациями запаздывающего момента

    calc_k(x, y, z, t, sx, sy, sz, vx, vy, vz, t2, &k, &r, &nx, &ny, &nz);

    double v_x = vx(t2);
    double v_y = vy(t2);
    double v_z = vz(t2);

    double w_x = wx(t2);
    double w_y = wy(t2);
    double w_z = wz(t2);

    double v2_c2 = (v_x*v_x + v_y*v_y + v_z*v_z)) / (c*c);
    double ra_c2 = r * (nx*w_x + ny*w_y + nz*w_z) / (c*c);
    
    (*E_x) = q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nx - v_x/c)/(r*r) - (k/r)*w_x/(c*c));
    (*E_y) = q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(ny - v_y/c)/(r*r) - (k/r)*w_y/(c*c));
    (*E_z) = q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nz - v_z/c)/(r*r) - (k/r)*w_z/(c*c));

    (*B_x) = -q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(ny*v_z - nz*v_y)/(r*r)/c + (ny*w_z - nz*w_y)*(k/r)/(c*c));
    (*B_y) = -q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nz*v_x - nx*v_z)/(r*r)/c + (nz*w_x - nx*w_z)*(k/r)/(c*c));
    (*B_z) = -q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nx*v_y - ny*v_x)/(r*r)/c + (nx*w_y - ny*w_x)*(k/r)/(c*c));
}

void test(double x)
{
    printf("test x = %lf\n", x);
}
