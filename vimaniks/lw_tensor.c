#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include "lw.h"
#include "lw_rotate.h"
#include "lw_vimanic.h"
#include "lw_tensor.h"

#include "stdlib.h"

#define _USE_MATH_DEFINES
#include <math.h>
#define Pi M_PI

#define Sq(x) ((x)*(x))

extern velocity c;

extern long double R_r;
extern long double R_l;

long double get_Omega_l() {return + cget_omega();}
long double get_Omega_r() {return - cget_omega();}

int ccalc_sum_F_t(int N, long double t_i,
                  long double alpha0_l,
                  long double alpha0_r,
                  long double * Fx,
                  long double * Fy,
                  long double * Fz,
                  long double * F_alpha_l,
                  long double * F_alpha_r,
                  long double * sum_rlagerror_square,
                  _Bool to_log)
{
    long double sign_a;
    long double sign_q;

    long double Alpha_l;
    long double Alpha_r;

    long double current_angle_l;
    long double current_angle_r;

    long double Omega_l = get_Omega_l();
    long double Omega_r = get_Omega_r();

    long double Xa, Ya, Za; // X_l, Y_l, Z_l
    long double Xq, Yq, Zq; // X_r, Y_r, Z_r

    long double vx_l, vy_l, vz_l;
    long double vx_r, vy_r, vz_r;

    long double Fx_l = (long double)(0.0);
    long double Fy_l = (long double)(0.0);
    long double Fz_l = (long double)(0.0);
    long double Fx_r = (long double)(0.0);
    long double Fy_r = (long double)(0.0);
    long double Fz_r = (long double)(0.0);

    if (F_alpha_l) {
        *F_alpha_l = (long double)(0.0);
    }

    if (F_alpha_r) {
        *F_alpha_r = (long double)(0.0);
    }

    long double
        fx_l = (long double)(0.0),
        fx_r = (long double)(0.0),
        fy_l = (long double)(0.0),
        fy_r = (long double)(0.0),
        fz_l = (long double)(0.0),
        fz_r = (long double)(0.0);
    long double f_alpha_l = (long double)(0.0), f_alpha_r = (long double)(0.0);

#ifndef USE_LEFT_CHARGE_ONLY
    long double E_x, E_y, E_z, B_x, B_y, B_z;
#endif

    long double fx_3l, fy_3l;
    long double fx_3r, fy_3r;

    coordinate rlagerror;
    *sum_rlagerror_square = (long double)(0.0);

    for(int i_a = 0; i_a < N; ++i_a) {
        for (int i_q = 0; i_q < N; ++i_q) {

            //i_l = i_a
            //i_r = i_q

            sign_a = (i_a%2)*2-1;
            sign_q = -((i_q%2)*2-1);

            Alpha_l = M_PI - (i_a * (2*M_PI)/N + alpha0_l);
            Alpha_r = i_q * (2*M_PI)/N + alpha0_r;

            Xa = sx(t_i, xc_l(), yc_l(), zc_l(), R_l, Omega_l, Alpha_l);
            Xq = sx(t_i, xc_r(), yc_r(), zc_r(), R_r, Omega_r, Alpha_r);

            Ya = sy(t_i, xc_l(), yc_l(), zc_l(), R_l, Omega_l, Alpha_l);
            Yq = sy(t_i, xc_r(), yc_r(), zc_r(), R_r, Omega_r, Alpha_r);

            Za = sz(t_i, xc_l(), yc_l(), zc_l(), R_l, Omega_l, Alpha_l);
            Zq = sz(t_i, xc_r(), yc_r(), zc_r(), R_r, Omega_r, Alpha_r);

            vx_l = vx(t_i, xc_l(), yc_l(), zc_l(), R_l, Omega_l, Alpha_l);
            vx_r = vx(t_i, xc_r(), yc_r(), zc_r(), R_r, Omega_r, Alpha_r);

            vy_l = vy(t_i, xc_l(), yc_l(), zc_l(), R_l, Omega_l, Alpha_l);
            vy_r = vy(t_i, xc_r(), yc_r(), zc_r(), R_r, Omega_r, Alpha_r);

            vz_l = vz(t_i, xc_l(), yc_l(), zc_l(), R_l, Omega_l, Alpha_l);
            vz_r = vz(t_i, xc_r(), yc_r(), zc_r(), R_r, Omega_r, Alpha_r);

            // торможение излучением (разложение функции Лагранжа
            // до членов третьего порядка)
            // ЛЛ2 (75,8)
            fx_3r =  2 * sign_q * sign_q / (3*c*c*c) * dot_wx(t_i, xc_r(), yc_r(), zc_r(), R_r, Omega_r, Alpha_r);
            fy_3r =  2 * sign_q * sign_q / (3*c*c*c) * dot_wy(t_i, xc_r(), yc_r(), zc_r(), R_r, Omega_r, Alpha_r);
            fx_3l =  2 * sign_a * sign_a / (3*c*c*c) * dot_wx(t_i, xc_l(), yc_l(), zc_l(), R_l, Omega_l, Alpha_l);
            fy_3l =  2 * sign_a * sign_a / (3*c*c*c) * dot_wy(t_i, xc_l(), yc_l(), zc_l(), R_l, Omega_l, Alpha_l);

            // релятивистское преобразование электрического поля
            // при переходе из системы координат,
            // связанной с зарядом в лабораторную систему координат
            // потому что согласно ЛЛ2 параграф 73
            // для точечного (не распределённого в объеме) заряда
            // формула (75,8) является точным выражением для обратного действия излучения
            // в той системе отсчёта, в которой заряд покоится
            // а интегрирование производится в лабораторной системе координат
            // поэтому необходимо в соответствии с преобразованиями Лоренца
            // формулы (26.35) (26.38) Фейнмановских лекций по физике, том 6
            // произвести преобразование формулы ЛЛ2 (75,8) из системы координат,
            // связанной с зарядом в лабораторную систему координат
            fx_3r /= sqrtl((long double)(1.0) - cget_vc()*cget_vc());
            fy_3r /= sqrtl((long double)(1.0) - cget_vc()*cget_vc());
            fx_3l /= sqrtl((long double)(1.0) - cget_vc()*cget_vc());
            fy_3l /= sqrtl((long double)(1.0) - cget_vc()*cget_vc());
#ifndef USE_LEFT_CHARGE_ONLY
            // поле создаваемое правым вращающимся зарядом
            // в области левого вращающегося заряда
            if (0 != electr_magnet(Xa, Ya, Za, t_i,
                  //sx, sy, sz, vx, vy, vz, wx, wy, wz,
                  sign_q,
                  &E_x, &E_y, &E_z, &B_x, &B_y, &B_z, &rlagerror, xc_r(), yc_r(), zc_r(), R_r, Omega_r, Alpha_r))
            {
                return -1;
            }

            *sum_rlagerror_square += Sq(rlagerror);

            // сила действующая на левый заряд
            // со стороны поля правого заряда
            fx_l = sign_a*(E_x + (vy_l*B_z-vz_l*B_y)/cget_c());
            fy_l = sign_a*(E_y + (vz_l*B_x-vx_l*B_z)/cget_c());
            fz_l = sign_a*(E_z + (vx_l*B_y-vy_l*B_x)/cget_c());
#endif
            // добавка обусловленная торможением излучением
            fx_l += fx_3l;
            fy_l += fy_3l;

            Fx_l += fx_l;
            Fy_l += fy_l;
            Fz_l += fz_l;

            if (F_alpha_l) {

                current_angle_l = Omega_l * t_i + Alpha_l;

                // Расчёт углового усилия действующего на левый заряд
                // со стороны поля правого заряда
                f_alpha_l  = fy_l * cosl(current_angle_l) - fx_l * sinl(current_angle_l);
                *F_alpha_l += f_alpha_l;
            }
#ifndef USE_LEFT_CHARGE_ONLY
            // поле создаваемое левым вращающимся зарядом
            // в области правого вращающегося заряда
            if (0 != electr_magnet(Xq, Yq, Zq, t_i,
                  //sx, sy, sz, vx, vy, vz, wx, wy, wz,
                  sign_a,
                  &E_x, &E_y, &E_z, &B_x, &B_y, &B_z, &rlagerror, xc_l(), yc_l(), zc_l(), R_l, Omega_l, Alpha_l))
            {
                return -1;
            }

            *sum_rlagerror_square += Sq(rlagerror);

            // сила действующая на правый заряд
            // со стороны поля левого заряда
            fx_r = sign_q*(E_x + (vy_r*B_z-vz_r*B_y)/cget_c());
            fy_r = sign_q*(E_y + (vz_r*B_x-vx_r*B_z)/cget_c());
            fz_r = sign_q*(E_z + (vx_r*B_y-vy_r*B_x)/cget_c());

            // добавка обусловленная торможением излучением
            fx_r += fx_3r;
            fy_r += fy_3r;

            Fx_r += fx_r;
            Fy_r += fy_r;
            Fz_r += fz_r;

            if (F_alpha_r) {
                current_angle_r = Omega_r * t_i + Alpha_r;

                // Расчёт углового усилия действующего на правый заряд
                // со стороны поля левого заряда
                f_alpha_r  = fy_r * cosl(current_angle_r) - fx_r * sinl(current_angle_r);
                *F_alpha_r += f_alpha_r;
            }
#endif
            if (to_log){
                printf("fx_l=%Lf fy_l=%Lf fx_r=%Lf fy_r=%Lf",fx_l, fy_l, fx_r, fy_r);
            }
        }
    }

    if (Fx)
    {
        *Fx = Fx_l + Fx_r;
    }

    if (Fy)
    {
        // Интегральная величина тяги в направлении оси y
        // printf("Fy_l = %Lf Fy_r = %Lf\n", Fy_l, Fy_r);
        *Fy = Fy_l + Fy_r;
    }

    if (Fz)
    {
        *Fz = Fz_l + Fz_r;
    }

    return 0;
}

int ccalc_Maxwells_stress_tensor(long double X_a, long double Y_a, long double Z_a, long double t_i,
                                 int N,
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
                               )
{
    long double sign_a;
    long double sign_q;

    long double Alpha_l;
    long double Alpha_r;

    long double Omega_l = get_Omega_l();
    long double Omega_r = get_Omega_r();

    //long double Xa, Ya, Za; // X_l, Y_l, Z_l
    //long double Xq, Yq, Zq; // X_r, Y_r, Z_r

    long double Ex = 0;
    long double Ey = 0;
    long double Ez = 0;

    long double Hx = 0;
    long double Hy = 0;
    long double Hz = 0;

    long double Ax = 0;
    long double Ay = 0;
    long double Az = 0;

    long double jx = 0;
    long double jy = 0;
    long double jz = 0;

    long double E_x, E_y, E_z, B_x, B_y, B_z, A_x, A_y, A_z;
    long double j_x, j_y, j_z;

    coordinate rlagerror;
    *sum_rlagerror_square = 0.0;

    for(int i_a = 0; i_a < N; ++i_a) {
        for (int i_q = 0; i_q < N; ++i_q) {
            //i_l = i_a
            //i_r = i_q

            sign_a = (i_a%2)*2-1;
            sign_q = -((i_q%2)*2-1);

            Alpha_l = M_PI - (i_a * (2*M_PI)/N + alpha0_l);
            Alpha_r = i_q * (2*M_PI)/N + alpha0_r;

            //printf("sign_a = %Lf sign_q = %Lf\n", sign_a, sign_q);
            //printf("Alpha_l = %Lf Alpha_r = %Lf\n", Alpha_l, Alpha_r);

#ifndef USE_LEFT_CHARGE_ONLY
            // поле в точке наблюдения создаваемое правым зарядом
            if (0 != electr_magnet_ex(X_a, Y_a, Z_a, t_i,
                  //sx, sy, sz, vx, vy, vz, wx, wy, wz,
                  sign_q,
                  &E_x, &E_y, &E_z,
                  &B_x, &B_y, &B_z,
                  &A_x, &A_y, &A_z,
                  &j_x, &j_y, &j_z,
                  &rlagerror, xc_r(), yc_r(), zc_r(), R_r, Omega_r, Alpha_r))
            {
                return -1;
            }

            *sum_rlagerror_square += Sq(rlagerror);

            Ex += E_x;
            Ey += E_y;
            Ez += E_z;

            Hx += B_x;
            Hy += B_y;
            Hz += B_z;

            Ax += A_x;
            Ay += A_y;
            Az += A_z;

    jx += j_x;
    jy += j_y;
    jz += j_z;

#endif

            // поле в точке наблюдения создаваемое левым зарядом
            if (0 != electr_magnet_ex(X_a, Y_a, Z_a, t_i,
                  //sx, sy, sz, vx, vy, vz, wx, wy, wz,
                  sign_a,
                  &E_x, &E_y, &E_z,
                  &B_x, &B_y, &B_z,
                  &A_x, &A_y, &A_z,
                  &j_x, &j_y, &j_z,
                  &rlagerror, xc_l(), yc_l(), zc_l(), R_l, Omega_l, Alpha_l))
            {
                return -1;
            }

            *sum_rlagerror_square += Sq(rlagerror);

            Ex += E_x;
            Ey += E_y;
            Ez += E_z;

            Hx += B_x;
            Hy += B_y;
            Hz += B_z;

            Ax += A_x;
            Ay += A_y;
            Az += A_z;

    jx += j_x;
    jy += j_y;
    jz += j_z;

        }
    }

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
    *j_n = (jx * cos_nx + jy * cos_ny + jz * cos_nz);

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
    long double * Txn, long double * Tyn, long double * Tzn,
    long double * Nx, long double * Ny, long double * Nz,
    long double * Sn,
    long double * En,
    long double * Hn,
    long double * An,
    long double * jn,
    long double * sum_rlagerror_square)
{
    return ccalc_Maxwells_stress_tensor(
        r*sinl(theta)*sinl(varphi), // X_a - в декартовой -> y в сферической
        r*cosl(theta),              // Y_a - в декартовой -> z в сферической
        r*sinl(theta)*cosl(varphi), // Z_a - в декартовой -> x в сферической
        t,
        1,                          // N
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
        An,
        jn,
        sum_rlagerror_square);
}

// Интегрируем в сферической системе координат,
// которая однако в соотвествие с принятыми в задаче
// наименованиями осей повернута так,
// что главная ось (зет в сферической системе)
// направлена вдоль оси икс декартовой системы

// направление векторов нормали к сферической воображаемой поверхности инвертировано - снаружи вовнутрь
int spherical_x_ccalc_Maxwells_stress_tensor(
    long double r, long double theta, long double varphi, long double t,
    long double * Txn, long double * Tyn, long double * Tzn,
    long double * Nx, long double * Ny, long double * Nz,
    long double * Sn,
    long double * En,
    long double * Hn,
    long double * An,
    long double * jn,
    long double * sum_rlagerror_square)
{
    return ccalc_Maxwells_stress_tensor(
        r*cosl(theta),              // X_a - в декартовой -> z в сферической
        r*sinl(theta)*cosl(varphi), // Y_a - в декартовой -> x в сферической
        r*sinl(theta)*sinl(varphi), // Z_a - в декартовой -> y в сферической
        t,
        1,                          // N
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
        An,
        jn,
        sum_rlagerror_square);
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
    long double * pTxn, long double * pTyn, long double * pTzn,
    long double * pNx, long double * pNy, long double * pNz,
    long double * pSn,
    long double * pEn,
    long double * pHn,
    long double * pAn,
    long double * pjn,
    long double * sum_rlagerror_square)
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
    long double jn;

    int ret = spherical_y_ccalc_Maxwells_stress_tensor(sphere_R, theta, varphi, t,
                                                     &Txn, &Tyn, &Tzn,
                                                     &Nx, &Ny, &Nz,
                                                     &Sn,
                                                     &En,
                                                     &Hn,
                                                     &An,
                                                     &jn,
                                                     sum_rlagerror_square);
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
    *pjn  = - sphere_R * sphere_R * sinl(theta) * jn;

    return ret;
}

int spherical_x_ccalc_Maxwells_stress_tensor_R_t(
    long double theta, long double varphi, long double t,
    long double * pTxn, long double * pTyn, long double * pTzn,
    long double * pNx, long double * pNy, long double * pNz,
    long double * pSn,
    long double * pEn,
    long double * pHn,
    long double * pAn,
    long double * pjn,
    long double * sum_rlagerror_square)
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
    long double jn;

    int ret = spherical_x_ccalc_Maxwells_stress_tensor(sphere_R, theta, varphi, t,
                                                     &Txn, &Tyn, &Tzn,
                                                     &Nx, &Ny, &Nz,
                                                     &Sn,
                                                     &En,
                                                     &Hn,
                                                     &An,
                                                     &jn,
                                                     sum_rlagerror_square);
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
    *pjn  = - sphere_R * sphere_R * sinl(theta) * jn;

    return ret;
}
