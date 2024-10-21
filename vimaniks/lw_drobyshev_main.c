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

distance R_lw(coordinate x, coordinate y, coordinate z, timevalue t)
{
    coordinate Xa = x;
    coordinate Ya = y;
    coordinate Za = z;
    timevalue ta = t;

    //set_c(1.0);
    distance r_lw = Rlw(Xa, Ya, Za, ta, sx, sy, sz, vx, vy, vz);

    return r_lw;
}

void main()
{
    double x = 5;
    double y = 5;
    double z = 0;
    double t = 7.5;

    double R = R_lw(x, y, z, t);

    printf("R = %f\n", R);
}