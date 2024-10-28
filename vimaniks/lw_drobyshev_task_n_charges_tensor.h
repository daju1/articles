#include "lw.h"
#include "lw_drobyshev_task.h"

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
    long double * A_n,
    long double * j_n,
    long double * sum_rlagerror_square
    );

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
    long double * An,
    long double * jn,
    long double * sum_rlagerror_square);

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
    long double * An,
    long double * jn,
    long double * sum_rlagerror_square);

void cset_sphere_R(long double R);
long double cget_sphere_R();

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
    long double * pAn,
    long double * pjn,
    long double * sum_rlagerror_square);

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
    long double * pAn,
    long double * pjn,
    long double * sum_rlagerror_square);
