#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include "lw.h"
#include "lw_drobyshev_task.h"
#include "lw_drobyshev_task_tensor.h"

#include "stdlib.h"

#define _USE_MATH_DEFINES
#include <math.h>
#define Pi M_PI

//const double c = 1.0;
const velocity vk = (const double)0.84;// finish velocity
const acceleration a = (const double)0.3; // acseleration
#define t0 ((vk)/(a)) // time of acseleration

coordinate sx(timevalue t) {
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

coordinate sy(timevalue t) {
    coordinate result;
    result = 0;
    return result;
}

coordinate sz(timevalue t) {
    coordinate result;
    result = 0;
    return result;
}

velocity vx(timevalue t) {
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

velocity vy(timevalue t) {
    velocity result;
    result = 0;
    return result;
}

velocity vz(timevalue t) {
    velocity result;
    result = 0;
    return result;
}


acceleration wx(timevalue t) {
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

acceleration wy(timevalue t) {
    acceleration result;
    result = 0;
    return result;
}

acceleration wz(timevalue t) {
    acceleration result;
    result = 0;
    return result;
}


dotacceleration dot_wx(timevalue t) {
    dotacceleration result;
    result = 0;
    return result;
}

dotacceleration dot_wy(timevalue t) {
    dotacceleration result;
    result = 0;
    return result;
}

dotacceleration dot_wz(timevalue t) {
    dotacceleration result;
    result = 0;
    return result;
}

distance R_lw(coordinate x, coordinate y, coordinate z, timevalue t)
{
    coordinate Xa = x;
    coordinate Ya = y;
    coordinate Za = z;
    timevalue ta = t;

    //set_c(1.0);
    distance r_lw;
    coordinate rlagerror;
    Rlw(Xa, Ya, Za, ta, sx, sy, sz, vx, vy, vz, &r_lw, &rlagerror);

    return r_lw;
}


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

    // printf("sum_rlagerror_square = %Le\n", sum_rlagerror_square);
    printf("t = %Le Tn = %Le %Le %Le N = %Le %Le %Le SEHAj =  %Le %Le %Le %Le %Le\n",
        t,
        Txn, Tyn, Tzn,
        Nx, Ny, Nz,
        Sn, En, Hn, An, jn);

    return jn;
}

void main()
{
    for ( int i = 0; i < 10; ++i)
    {
        double t = 1+0.1*i;

        long double xc = sx(t);
        long double theta = M_PI / 4;
        long double varphi = 0.1;

        cset_sphere_R(1.0);

        long double En = spherical_x_calc_En_R_t(xc, theta, varphi, t);

        printf("xc = %Lf\n", xc);
        // printf("En = %Lf\n", En);
    }
}
