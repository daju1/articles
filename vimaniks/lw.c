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

//sgs 

static velocity c = (double)(299792458 * 100);

void set_c(double _c)
{
    c = _c;
}

const distance R_r = (double)10;
const distance R_l = (double)10;

const double v_c = (double) 0.984;

double get_omega()
{
    return v_c * c / R_r;
}

typedef coordinate (*Coordinate)(timevalue t_zap);
typedef velocity (*Velocity)(timevalue t_zap);
typedef acceleration (*Acceleration)(timevalue t_zap);

// расчет итерациями запаздывающего момента

timespan Epsilon = 1.e-6;// # погрешность
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
    timevalue t2 = t - 0.1;
    //printf("____ x = %f y = %f\n", x, y);
    //printf("____ t1 = %f t2 = %f\n", t1, t2);
    
    int n = 0;

    while (fabs(t1 - t2) > Epsilon) {
        t1 = t2;
        double dd = (x - sx(t1))*(x - sx(t1)) + (y - sy(t1))*(y - sy(t1)) + (z - sz(t1))*(z - sz(t1));
        //printf("dd = %f\n", dd);
        double d = sqrt(dd);
        //printf("d = %f\n", d);
        t2 = t - d / c;
        //printf("t1 = %f t2 = %f\n", t1, t2);
        if (++n > 100)
            break;
    }
    //printf(">>>>> t1 = %f t2 = %f\n", t1, t2);

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
        (*r) = sqrt((x - sx(t))*(x - sx(t)) + (y - sy(t))*(y - sy(t)) + (z - sz(t))*(z - sz(t)));
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

    double v2_c2 = (vx(t2)*vx(t2) + vy(t2)*vy(t2) + vz(t2)*vz(t2)) / (c*c);
    double ra_c2 = r * (nx*wx(t2) + ny*wy(t2) + nz*wz(t2)) / (c*c);
    
    (*E_x) = q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nx - vx(t2)/c)/(r*r) - (k/r)*wx(t2)/(c*c));
    (*E_y) = q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(ny - vy(t2)/c)/(r*r) - (k/r)*wy(t2)/(c*c));
    (*E_z) = q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nz - vz(t2)/c)/(r*r) - (k/r)*wz(t2)/(c*c));

    (*B_x) = -q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(ny*vz(t2) - nz*vy(t2))/(r*r)/c + (ny*wz(t2) - nz*wy(t2))*(k/r)/(c*c));
    (*B_y) = -q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nz*vx(t2) - nx*vz(t2))/(r*r)/c + (nz*wx(t2) - nx*wz(t2))*(k/r)/(c*c));
    (*B_z) = -q*(1.0/(k*k*k))*((1.0 - v2_c2 + ra_c2)*(nx*vy(t2) - ny*vx(t2))/(r*r)/c + (nx*wy(t2) - ny*wx(t2))*(k/r)/(c*c));
}

void test(double x)
{
    printf("test x = %lf\n", x);
}
