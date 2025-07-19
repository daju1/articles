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

const long double L = 2;     // dipole length
const long double X = 5;     // dipole coordinate
const long double d = 0.1;   // dipole friction aplitude
const long double omega = 10;

coordinate s1x(timevalue t) {
    coordinate s1 = d*sin(omega*t) - L/2 + X;
    return s1;
}
coordinate s2x(timevalue t) {
    coordinate s2 = -d*sin(omega*t) + L/2 + X;
    return s2;
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
    velocity v1 = d*omega*cos(omega*t);
    return v1;
}

velocity v2x(timevalue t) {
    velocity v2 = -d*omega*cos(omega*t);
    return v2;
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
    acceleration w1 = -d*omega*omega*sin(omega*t);
    return w1;
}
acceleration w2x(timevalue t) {
    acceleration w2 = d*omega*omega*sin(omega*t);
    return w2;
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
    dotacceleration dot_w1 = -d*omega*omega*omega*cos(omega*t);
    return dot_w1;
}dotacceleration dot_w2x(timevalue t) {
    dotacceleration dot_w2 = d*omega*omega*omega*cos(omega*t);
    return dot_w2;
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

        long double xc = s1x(t);
        long double theta = M_PI / 4;
        long double varphi = 0.1;

        cset_sphere_R(1.0);

        long double En = spherical_x_calc_En_R_t(xc, theta, varphi, t);

        printf("xc = %Lf\n", xc);
        // printf("En = %Lf\n", En);
    }
}
