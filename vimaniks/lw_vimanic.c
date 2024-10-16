#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include "lw.h"
#include "lw_rotate.h"
#include "lw_vimanic.h"

#include "stdlib.h"


#define _USE_MATH_DEFINES
#include <math.h>
#define Pi M_PI

#define Sq(x) ((x)*(x))

long double R_r = 1;
long double R_l = 1;
static long double S = 0.05;
static long double v_c = 0.8;

void cset_vc(long double vc)
{
    v_c = vc;
}

long double cget_vc()
{
    return v_c;
}

long double cget_omega() {
    return v_c * cget_c() / cget_R_r();
}

//# centers of circles
#ifdef USE_LEFT_CENTERED_CHARGE_ONLY
long double xc_l() { return 0.0; }
#else
long double xc_l() { return -S/2 - R_l; }
#endif

long double xc_r() { return  S/2 + R_r; }


long double yc_l() { return 0; }
long double yc_r() { return 0; }

long double zc_l() { return 0; }
long double zc_r() { return 0; }

extern velocity c;

long double sx(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha)
{
    long double result;
    long double current_angle = omega * t + alpha;
    result = xc + R*cosl(current_angle);
    return result;
}

long double sy(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha)
{
#ifdef TEST_WITH_X_DIRECTED_DIPOLE
    return yc;
#endif
    long double result;
    long double current_angle = omega * t + alpha;
    result = yc + R*sinl(current_angle);
    return result;
}

long double sz(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha)
{
    long double result;
    result = zc;
    return result;
}

long double vx(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha)
{
    long double result;
    long double current_angle = omega * t + alpha;
    result = -omega*R*sinl(current_angle);
    return result;
}

long double vy(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha)
{
#ifdef TEST_WITH_X_DIRECTED_DIPOLE
    return (long double)(0.0);
#endif
    long double result;
    long double current_angle = omega * t + alpha;
    result = omega*R*cosl(current_angle);
    return result;
}

long double vz(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha)
{
    long double result;
    result = 0;
    return result;
}

long double wx(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha)
{
    long double result;
    long double current_angle = omega * t + alpha;
    result = -omega*omega*R*cosl(current_angle);
    return result;
}
long double wy(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha)
{
#ifdef TEST_WITH_X_DIRECTED_DIPOLE
    return (long double)(0.0);
#endif
    long double result;
    long double current_angle = omega * t + alpha;
    result = -omega*omega*R*sinl(current_angle);
    return result;
}
long double wz(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha)
{
    long double result;
    result = 0;
    return result;
}

long double dot_wx(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha)
{
    long double result;
    long double current_angle = omega * t + alpha;
    result = omega*omega*omega*R*sinl(current_angle);
    return result;
}
long double dot_wy(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha)
{
#ifdef TEST_WITH_X_DIRECTED_DIPOLE
    return (long double)(0.0);
#endif
    long double result;
    long double current_angle = omega * t + alpha;
    result = -omega*omega*omega*R*cosl(current_angle);
    return result;
}
long double dot_wz(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha)
{
    long double result;
    result = 0;
    return result;
}


long double cget_S()  { return S;}
long double cget_R_l(){ return R_l; }
long double cget_R_r(){ return R_r; }

void cset_S(long double s)  { S = s;}
void cset_R_l(long double r){ R_l = r; }
void cset_R_r(long double r){ R_r = r; }

long double cget_xc_l(){ return xc_l(); }
long double cget_xc_r(){ return xc_r(); }

long double cget_yc_l(){ return yc_l(); }
long double cget_yc_r(){ return yc_r(); }

long double cget_zc_l(){ return zc_l(); }
long double cget_zc_r(){ return zc_r(); }
