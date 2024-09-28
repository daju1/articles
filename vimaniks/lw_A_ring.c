#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include "lw.h"
#include "lw_rotate.h"
#include "lw_A_ring.h"

#include "stdlib.h"

#define _USE_MATH_DEFINES
#include <math.h>
#define Pi M_PI

#define Sq(x) (x)*(x)

static long double R_ring = 1;
static long double v_c = 0.8;

void cset_vc(long double vc)
{
    v_c = vc;
}

long double cget_vc()
{
    return v_c;
}

// center of ring
long double xc_ring = 0.0;
long double yc_ring = 0.0;
long double zc_ring = 0.0;

long double sx(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha)
{
    double result;
    double current_angle = omega * t + alpha;
    result = xc + R*cos(current_angle);
    return result;
}

long double sy(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha)
{
    long double result;
    long double current_angle = omega * t + alpha;
    result = yc + R*sin(current_angle);
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
    result = -omega*R*sin(current_angle);
    return result;
}

long double vy(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha)
{
    long double result;
    long double current_angle = omega * t + alpha;
    result = omega*R*cos(current_angle);
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
    result = -omega*omega*R*cos(current_angle);
    return result;
}
long double wy(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha)
{
    double result;
    double current_angle = omega * t + alpha;
    result = -omega*omega*R*sin(current_angle);
    return result;
}
long double wz(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha)
{
    long double result;
    result = 0;
    return result;
}

int ccalc_A_ring(long double X_a, long double Y_a, long double Z_a, long double t,
                 charge q, long double alpha0, long double * A_phi, coordinate * rlagerror)
{
    long double omega = cget_omega();

    long double A_x, A_y, A_z;

    // A_lw - векторный потенциал Лиенара Вихерта
    int ret = Alw(X_a, Y_a, Z_a, t,
        sx, sy, sz, vx, vy, vz,
        q,
        &A_x, &A_y, &A_z, rlagerror,
        xc_ring, yc_ring, zc_ring,
        R_ring, omega, alpha0);
    if (0 != ret)
    {
        return ret;
    }

    // угловая координата точки наблюдения в полярных координатах с центром в центре токового кольца
    long double phi_a = atan2(Y_a-yc_ring, X_a-xc_ring);

    // направление фитой оси полярной системы координат в точке наблюдения
    long double phi_phi = phi_a + M_PI/2;
    if (phi_phi > 2*M_PI)
        phi_phi -= 2*M_PI;
    if (phi_phi < -2*M_PI)
        phi_phi += 2*M_PI;

    long double phi_A = atan2(A_y, A_x);

    // angle between vector A and polar axis with center of rotated charge
    long double beta  =  phi_phi - phi_A;
    *A_phi = sqrt(Sq(A_x) + Sq(A_y)) * cos(beta);

    return ret;
}

long double cget_R_ring(){ return R_ring; }
void cset_R_ring(long double r){ R_ring = r; }

long double cget_xc_ring(){ return xc_ring; }
long double cget_yc_ring(){ return yc_ring; }
long double cget_zc_ring(){ return zc_ring; }

void cset_xc_ring(long double x){ return xc_ring = x; }
void cset_yc_ring(long double y){ return yc_ring = y; }
void cset_zc_ring(long double z){ return zc_ring = z; }


long double cget_omega() {
    return v_c * cget_c() / cget_R_ring();
}
