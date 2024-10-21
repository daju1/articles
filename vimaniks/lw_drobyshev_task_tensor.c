#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include "lw.h"
#include "lw_drobyshev_task.h"

#define _USE_MATH_DEFINES
#include <math.h>
#define Pi M_PI

extern  velocity c;

int ccalc_Maxwells_stress_tensor(
    long double X_a, long double Y_a, long double Z_a, long double t_i,
    Coordinate sx, Coordinate sy, Coordinate sz,
    Velocity vx, Velocity vy, Velocity vz,
    Acceleration wx, Acceleration wy, Acceleration wz,
    long double cos_nx, long double cos_ny, long double cos_nz,
    long double alpha0_l,
    long double alpha0_r,
    long double * T_xn,
    long double * T_yn,
    long double * T_zn,
    long double * N_x,
    long double * N_y,
    long double * N_z,
    long double * S_n,
    long double * E_n,
    long double * H_n,
    long double * A_n
    )
{
    long double _q = 1.0;

    long double Ex = 0;
    long double Ey = 0;
    long double Ez = 0;

    long double Hx = 0;
    long double Hy = 0;
    long double Hz = 0;

    long double Ax = 0;
    long double Ay = 0;
    long double Az = 0;

    long double E_x, E_y, E_z, B_x, B_y, B_z, A_x, A_y, A_z;

    // поле в точке наблюдения создаваемое левым зарядом
    electr_magnet_ex(X_a, Y_a, Z_a, t_i,
            sx, sy, sz, vx, vy, vz, wx, wy, wz,
            _q,
            &E_x, &E_y, &E_z,
            &B_x, &B_y, &B_z,
            &A_x, &A_y, &A_z);

    Ex += E_x;
    Ey += E_y;
    Ez += E_z;

    Hx += B_x;
    Hy += B_y;
    Hz += B_z;

    Ax += A_x;
    Ay += A_y;
    Az += A_z;


    // ЛЛ2 (33,3)

    long double sigma_xx = ((long double)(1.0))/(8*M_PI)*( Ex*Ex + Hx*Hx - Ey*Ey - Ez*Ez - Hy*Hy - Hz*Hz);
    long double sigma_yy = ((long double)(1.0))/(8*M_PI)*( Ey*Ey + Hy*Hy - Ez*Ez - Ex*Ex - Hz*Hz - Hx*Hx);
    long double sigma_zz = ((long double)(1.0))/(8*M_PI)*( Ez*Ez + Hz*Hz - Ex*Ex - Ey*Ey - Hx*Hx - Hy*Hy);

    long double sigma_xy = ((long double)(1.0))/(4*M_PI)*( Ex*Ey + Hx*Hy );
    long double sigma_xz = ((long double)(1.0))/(4*M_PI)*( Ex*Ez + Hx*Hz );
    long double sigma_yz = ((long double)(1.0))/(4*M_PI)*( Ey*Ez + Hy*Hz );

    long double sigma_yx = ((long double)(1.0))/(4*M_PI)*( Ey*Ex + Hy*Hx );
    long double sigma_zx = ((long double)(1.0))/(4*M_PI)*( Ez*Ex + Hz*Hx );
    long double sigma_zy = ((long double)(1.0))/(4*M_PI)*( Ez*Ey + Hz*Hy );

    //T = [[sigma_xx, sigma_xy, sigma_xz],
    //     [sigma_yx, sigma_yy, sigma_yz],
    //     [sigma_zx, sigma_zy, sigma_zz]]

    // Тамм параграф 33 формула (33.5)
    // сила натяжения действующая на площадку поверхности интегрирования
    // со стороны поля создаваемого вращающимися зарядами
    // Интегральная величина количества имульса электромагнитного поля,
    // вытекающего в единицу времени из замкнутого объёма через площадку ЛЛ2 32.14
    *T_xn = (sigma_xx * cos_nx + sigma_xy * cos_ny + sigma_xz * cos_nz);
    *T_yn = (sigma_yx * cos_nx + sigma_yy * cos_ny + sigma_yz * cos_nz);
    *T_zn = (sigma_zx * cos_nx + sigma_zy * cos_ny + sigma_zz * cos_nz);

    // Момент сил натяжения, приложенных к поверхности $S$ объема $V$
    *N_x = Y_a * (*T_zn) - Z_a * (*T_yn);
    *N_y = Z_a * (*T_xn) - X_a * (*T_zn);
    *N_z = X_a * (*T_yn) - Y_a * (*T_xn);

    // ЛЛ2 (31,2)

    long double Sx = (c)/(4*M_PI)*(Ey * Hz - Ez * Hy);
    long double Sy = (c)/(4*M_PI)*(Ez * Hx - Ex * Hz);
    long double Sz = (c)/(4*M_PI)*(Ex * Hy - Ey * Hx);

    // складова вектора Пойнтінга,
    // перпендикулярна до поверхні
    *S_n = (Sx * cos_nx + Sy * cos_ny + Sz * cos_nz);
    *E_n = (Ex * cos_nx + Ey * cos_ny + Ez * cos_nz);
    *H_n = (Hx * cos_nx + Hy * cos_ny + Hz * cos_nz);
    *A_n = (Ax * cos_nx + Ay * cos_ny + Az * cos_nz);

    return 0;
}


// Интегрируем в сферической системе координат,
// которая однако в соотвествие с принятыми в задаче
// наименованиями осей повернута так,
// что главная ось (зет в сферической системе)
// направлена вдоль оси игрек декартовой системы

// направление векторов нормали к сферической воображаемой поверхности инвертировано - снаружи вовнутрь
int spherical_y_ccalc_Maxwells_stress_tensor(
    long double r, long double theta, long double varphi, long double t,
    Coordinate sx, Coordinate sy, Coordinate sz,
    Velocity vx, Velocity vy, Velocity vz,
    Acceleration wx, Acceleration wy, Acceleration wz,
    long double * Txn, long double * Tyn, long double * Tzn,
    long double * Nx, long double * Ny, long double * Nz,
    long double * Sn,
    long double * En,
    long double * Hn,
    long double * An)
{
    return ccalc_Maxwells_stress_tensor(
        r*sinl(theta)*sinl(varphi), // X_a - в декартовой -> y в сферической
        r*cosl(theta),              // Y_a - в декартовой -> z в сферической
        r*sinl(theta)*cosl(varphi), // Z_a - в декартовой -> x в сферической
        t,
        sx, sy, sz, vx, vy, vz, wx, wy, wz,
        - sinl(theta)*sinl(varphi), // cos_nx - i в декартовой -> j в сферической
        - cosl(theta),              // cos_ny - j в декартовой -> k в сферической
        - sinl(theta)*cosl(varphi), // cos_nz - k в декартовой -> i в сферической
        0,
        0,
        Txn, Tyn, Tzn,
        Nx, Ny, Nz,
        Sn,
        En,
        Hn,
        An);
}

// Интегрируем в сферической системе координат,
// которая однако в соотвествие с принятыми в задаче
// наименованиями осей повернута так,
// что главная ось (зет в сферической системе)
// направлена вдоль оси икс декартовой системы

// направление векторов нормали к сферической воображаемой поверхности инвертировано - снаружи вовнутрь
int spherical_x_ccalc_Maxwells_stress_tensor(long double xc,
    long double r, long double theta, long double varphi, long double t,
    Coordinate sx, Coordinate sy, Coordinate sz,
    Velocity vx, Velocity vy, Velocity vz,
    Acceleration wx, Acceleration wy, Acceleration wz,
    long double * Txn, long double * Tyn, long double * Tzn,
    long double * Nx, long double * Ny, long double * Nz,
    long double * Sn,
    long double * En,
    long double * Hn,
    long double * An)
{
    return ccalc_Maxwells_stress_tensor(
        xc + r*cosl(theta),         // X_a - в декартовой -> z в сферической
        r*sinl(theta)*cosl(varphi), // Y_a - в декартовой -> x в сферической
        r*sinl(theta)*sinl(varphi), // Z_a - в декартовой -> y в сферической
        t,
        sx, sy, sz, vx, vy, vz, wx, wy, wz,
        - cosl(theta),              // cos_nx - i в декартовой -> k в сферической
        - sinl(theta)*cosl(varphi), // cos_ny - j в декартовой -> i в сферической
        - sinl(theta)*sinl(varphi), // cos_nz - k в декартовой -> j в сферической
        0,
        0,
        Txn, Tyn, Tzn,
        Nx, Ny, Nz,
        Sn,
        En,
        Hn,
        An);
}
long double sphere_R;
void cset_sphere_R(long double R)
{
    sphere_R = R;
}

long double cget_sphere_R()
{
    return sphere_R;
}

int spherical_y_ccalc_Maxwells_stress_tensor_R_t(
    long double theta, long double varphi, long double t,
    Coordinate sx, Coordinate sy, Coordinate sz,
    Velocity vx, Velocity vy, Velocity vz,
    Acceleration wx, Acceleration wy, Acceleration wz,
    long double * pTxn, long double * pTyn, long double * pTzn,
    long double * pNx, long double * pNy, long double * pNz,
    long double * pSn,
    long double * pEn,
    long double * pHn,
    long double * pAn)
{
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

    int ret = spherical_y_ccalc_Maxwells_stress_tensor(sphere_R, theta, varphi, t,
        sx, sy, sz, vx, vy, vz, wx, wy, wz,
        &Txn, &Tyn, &Tzn,
        &Nx, &Ny, &Nz,
        &Sn,
        &En,
        &Hn,
        &An);

    *pTxn = sphere_R * sphere_R * sinl(theta) * Txn;
    *pTyn = sphere_R * sphere_R * sinl(theta) * Tyn;
    *pTzn = sphere_R * sphere_R * sinl(theta) * Tzn;

    *pNx = sphere_R * sphere_R * sinl(theta) * Nx;
    *pNy = sphere_R * sphere_R * sinl(theta) * Ny;
    *pNz = sphere_R * sphere_R * sinl(theta) * Nz;

    *pSn  = - sphere_R * sphere_R * sinl(theta) * Sn; // берём количество энергии излучения
    // протекающей через поверхность воображаемой сферы со знаком минус, потому что
    // направление векторов нормали к поверхности внутри функции spherical_ccalc_Maxwells_stress_tensor
    // инвертировано

    *pEn  = - sphere_R * sphere_R * sinl(theta) * En;
    *pHn  = - sphere_R * sphere_R * sinl(theta) * Hn;
    *pAn  = - sphere_R * sphere_R * sinl(theta) * An;

    return ret;
}

int spherical_x_ccalc_Maxwells_stress_tensor_R_t(long double xc,
    long double theta, long double varphi, long double t,
    Coordinate sx, Coordinate sy, Coordinate sz,
    Velocity vx, Velocity vy, Velocity vz,
    Acceleration wx, Acceleration wy, Acceleration wz,
    long double * pTxn, long double * pTyn, long double * pTzn,
    long double * pNx, long double * pNy, long double * pNz,
    long double * pSn,
    long double * pEn,
    long double * pHn,
    long double * pAn)
{
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

    int ret = spherical_x_ccalc_Maxwells_stress_tensor(xc, sphere_R, theta, varphi, t,
        sx, sy, sz, vx, vy, vz, wx, wy, wz,
        &Txn, &Tyn, &Tzn,
        &Nx, &Ny, &Nz,
        &Sn,
        &En,
        &Hn,
        &An);

    *pTxn = sphere_R * sphere_R * sinl(theta) * Txn;
    *pTyn = sphere_R * sphere_R * sinl(theta) * Tyn;
    *pTzn = sphere_R * sphere_R * sinl(theta) * Tzn;

    *pNx = sphere_R * sphere_R * sinl(theta) * Nx;
    *pNy = sphere_R * sphere_R * sinl(theta) * Ny;
    *pNz = sphere_R * sphere_R * sinl(theta) * Nz;

    *pSn  = - sphere_R * sphere_R * sinl(theta) * Sn; // берём количество энергии излучения
    // протекающей через поверхность воображаемой сферы со знаком минус, потому что
    // направление векторов нормали к поверхности внутри функции spherical_ccalc_Maxwells_stress_tensor
    // инвертировано
    *pEn  = - sphere_R * sphere_R * sinl(theta) * En;
    *pHn  = - sphere_R * sphere_R * sinl(theta) * Hn;
    *pAn  = - sphere_R * sphere_R * sinl(theta) * An;

    return ret;
}
