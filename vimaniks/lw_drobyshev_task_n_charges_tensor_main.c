#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include "lw.h"
#include "lw_drobyshev_task.h"
#include "lw_drobyshev_task_n_charges_tensor.h"

#include "stdlib.h"

#define _USE_MATH_DEFINES
#include <math.h>
#define Pi M_PI

//const double c = 1.0;
const velocity vk = (const double)0.84;// finish velocity
const acceleration a = (const double)0.3; // acseleration
#define t0 ((vk)/(a)) // time of acseleration

coordinate s1x(timevalue t) {
    coordinate result;
    if (t < 0)
    {
        result = 0;
    }
    else if (t < t0)
    {
        result = a * t*t / 2;
    }
    else
    {
        result = vk * t - a * t0*t0 / 2;
    }

    return result;
}
coordinate s2x(timevalue t) {
    coordinate result;
    if (t < 0)
    {
        result = 0;
    }
    else if (t < t0)
    {
        result = -a * t*t / 2;
    }
    else
    {
        result = -vk * t + a * t0*t0 / 2;
    }

    return result;
}

Coordinate sx[] = {s1x, s2x};

coordinate s1y(timevalue t) {
    coordinate result;
    result = 0;
    return result;
}
coordinate s2y(timevalue t) {
    coordinate result;
    result = 0;
    return result;
}
Coordinate sy[] = {s1y, s2y};

coordinate s1z(timevalue t) {
    coordinate result;
    result = 0;
    return result;
}
coordinate s2z(timevalue t) {
    coordinate result;
    result = 0;
    return result;
}
Coordinate sz[] = {s1z, s2z};

velocity v1x(timevalue t) {
    velocity result;
    if (t < 0)
    {
        result = 0;
    }
    else if (t < t0)
    {
        result = a * t;
    }
    else
    {
        result = vk;
    }

    return result;
}

velocity v2x(timevalue t) {
    velocity result;
    if (t < 0)
    {
        result = 0;
    }
    else if (t < t0)
    {
        result = -a * t;
    }
    else
    {
        result = vk;
    }

    return result;
}

Velocity vx[] = {v1x, v2x};

velocity v1y(timevalue t) {
    velocity result;
    result = 0;
    return result;
}
velocity v2y(timevalue t) {
    velocity result;
    result = 0;
    return result;
}

Velocity vy[] = {v1y, v2y};

velocity v1z(timevalue t) {
    velocity result;
    result = 0;
    return result;
}
velocity v2z(timevalue t) {
    velocity result;
    result = 0;
    return result;
}

Velocity vz[] = {v1z, v2z};


acceleration w1x(timevalue t) {
    acceleration result;
    if (t < 0)
    {
        result = 0;
    }
    else if (t < t0)
    {
        result = a;
    }
    else
    {
        result = 0;
    }

    return result;
}
acceleration w2x(timevalue t) {
    acceleration result;
    if (t < 0)
    {
        result = 0;
    }
    else if (t < t0)
    {
        result = -a;
    }
    else
    {
        result = 0;
    }

    return result;
}

Acceleration wx[] = {w1x, w2x};

acceleration w1y(timevalue t) {
    acceleration result;
    result = 0;
    return result;
}
acceleration w2y(timevalue t) {
    acceleration result;
    result = 0;
    return result;
}
Acceleration wy[] = {w1y, w2y};

acceleration w1z(timevalue t) {
    acceleration result;
    result = 0;
    return result;
}
acceleration w2z(timevalue t) {
    acceleration result;
    result = 0;
    return result;
}
Acceleration wz[] = {w1z, w2z};


dotacceleration dot_w1x(timevalue t) {
    dotacceleration result;
    result = 0;
    return result;
}dotacceleration dot_w2x(timevalue t) {
    dotacceleration result;
    result = 0;
    return result;
}
DotAcceleration dot_wx[] = {dot_w1x, dot_w2x};

dotacceleration dot_w1y(timevalue t) {
    dotacceleration result;
    result = 0;
    return result;
}
dotacceleration dot_w2y(timevalue t) {
    dotacceleration result;
    result = 0;
    return result;
}
DotAcceleration dot_wy[] = {dot_w1y, dot_w2y};

dotacceleration dot_w1z(timevalue t) {
    dotacceleration result;
    result = 0;
    return result;
}
dotacceleration dot_w2z(timevalue t) {
    dotacceleration result;
    result = 0;
    return result;
}
DotAcceleration dot_wz[] = {dot_w1z, dot_w2z};

distance R_lw(coordinate x, coordinate y, coordinate z, timevalue t)
{
    coordinate Xa = x;
    coordinate Ya = y;
    coordinate Za = z;
    timevalue ta = t;

    //set_c(1.0);
    distance r_lw;
    coordinate rlagerror;
    Rlw(Xa, Ya, Za, ta, s1x, s1y, s1z, v1x, v1y, v1z, &r_lw, &rlagerror);

    return r_lw;
}


int n_charges = 2;
charge q[] = {(long double)1.0, (long double)-1.0};

long double spherical_x_calc_En_R_t (long double xc, long double theta, long double varphi, long double t)
{
    long double _xc = xc;
    long double _theta = theta;
    long double _varphi = varphi;
    long double _t = t;

    long double Txn;
    long double Tyn;
    long double Tzn;
    long double Nx;
    long double Ny;
    long double Nz;
    long double Sn;
    long double En;
    long double Hn;
    long double An;
    long double jn;

    coordinate sum_rlagerror_square;

    int ret = spherical_x_ccalc_Maxwells_stress_tensor_R_t(_xc, _theta, _varphi, _t,
        n_charges, q,
        sx, sy, sz, vx, vy, vz, wx, wy, wz,
        dot_wx, dot_wy, dot_wz,
        &Txn,
        &Tyn,
        &Tzn,
        &Nx,
        &Ny,
        &Nz,
        &Sn,
        &En,
        &Hn,
        &An,
        &jn,
        &sum_rlagerror_square);

    printf("sum_rlagerror_square = %Le\n", sum_rlagerror_square);

    return jn;
}

void main()
{
    double x = 5;
    double y = 5;
    double z = 0;
    double t = 7.5;

    double R = R_lw(x, y, z, t);

    long double xc = s1x(t);
    long double theta = M_PI / 4;
    long double varphi = 0.1;

    cset_sphere_R(1.0);

    long double En = spherical_x_calc_En_R_t(xc, theta, varphi, t);

    printf("R = %f\n", R);
    printf("xc = %Lf\n", xc);
    printf("En = %Lf\n", En);
}
