#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include <gsl/gsl_math.h>

#include "pendulum_lib.h"

#include "stdlib.h"

#define _USE_MATH_DEFINES
#include <math.h>
#define Pi M_PI

#define Sq(x) ((x)*(x))

extern velocity c;

extern gsl_odeiv2_driver * driver_left;
extern gsl_odeiv2_driver * driver_right;

extern struct pendulum pendulum_left;
extern struct pendulum pendulum_right;


int ccalc_sum_F_t(timevalue t_i,
                  force * Fx,
                  force * Fy,
                  force * Fz,
                  force * F_alpha_l,
                  force * F_alpha_r,
                  double * sum_rlagerror_sqare,
                  _Bool to_log)
{
    charge sign_a;
    charge sign_q;

    angle current_angle_l;
    angle current_angle_r;

    coordinate X_l, Y_l, Z_l; //Xa, Ya, Za
    coordinate X_r, Y_r, Z_r; //Xq, Yq, Zq

    velocity vx_l, vy_l, vz_l;
    velocity vx_r, vy_r, vz_r;

    acceleration wx_l, wy_l, wz_l;
    acceleration wx_r, wy_r, wz_r;

    double dot_wx_l, dot_wy_l, dot_wz_l;
    double dot_wx_r, dot_wy_r, dot_wz_r;

    force Fx_l = (force)(0.0);
    force Fy_l = (force)(0.0);
    force Fz_l = (force)(0.0);
    force Fx_r = (force)(0.0);
    force Fy_r = (force)(0.0);
    force Fz_r = (force)(0.0);

    if (F_alpha_l) {
        *F_alpha_l = (force)(0.0);
    }

    if (F_alpha_r) {
        *F_alpha_r = (force)(0.0);
    }

    double
        fx_l = (force)(0.0),
        fx_r = (force)(0.0),
        fy_l = (force)(0.0),
        fy_r = (force)(0.0),
        fz_l = (force)(0.0),
        fz_r = (force)(0.0);
    force f_alpha_l = (force)(0.0), f_alpha_r = (force)(0.0);

#ifndef USE_LEFT_CHARGE_ONLY
    field E_x, E_y, E_z, B_x, B_y, B_z;
#endif

    double fx_3l, fy_3l;
    double fx_3r, fy_3r;

    coordinate rlagerror;
    *sum_rlagerror_sqare = (long double)(0.0);

    //i_l = i_a
    //i_r = i_q

    charge q_left = -1.0;
    charge q_right = +1.0;

    sign_a = q_left;
    sign_q = q_right;

    // определяем текущее состояние левого и правого маятника в момент t_i

    {
        double Momenta, q;

        double dot_q;   // dphi/dt
        double ddot_q;  // d2phi/dt2
        double dddot_q; // d3phi/dt3

        apply(&pendulum_left, driver_left,
            t_i, NULL, &Momenta, &q,
            &dot_q, &ddot_q, &dddot_q,
            &X_l, &Y_l,
            &vx_l, &vy_l,
            &wx_l, &wy_l,
            &dot_wx_l, &dot_wy_l
        );
        current_angle_l = q;
    }

    {
        double Momenta, q;

        double dot_q;   // dphi/dt
        double ddot_q;  // d2phi/dt2
        double dddot_q; // d3phi/dt3

        apply(&pendulum_right, driver_right,
            t_i, NULL, &Momenta, &q,
            &dot_q, &ddot_q, &dddot_q,
            &X_r, &Y_r,
            &vx_r, &vy_r,
            &wx_r, &wy_r,
            &dot_wx_r, &dot_wy_r
        );
        current_angle_r = q;
    }

    // торможение излучением (разложение функции Лагранжа
    // до членов третьего порядка)
    // ЛЛ2 (75,8)
    fx_3r =  2 * sign_q * sign_q / (3*c*c*c) * dot_wx_r;
    fy_3r =  2 * sign_q * sign_q / (3*c*c*c) * dot_wy_r;
    fx_3l =  2 * sign_a * sign_a / (3*c*c*c) * dot_wx_l;
    fy_3l =  2 * sign_a * sign_a / (3*c*c*c) * dot_wy_l;

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

    // но тут возникает вопрос как правильно выбрать v/c
#if 1
    // или основываясь на общей скорости
    double vc_l = sqrt(Sq(vx_l) + Sq(vy_l))/c;
    double vc_r = sqrt(Sq(vx_r) + Sq(vy_r))/c;

    fx_3r /= sqrtl((long double)(1.0) - vc_r*vc_r);
    fy_3r /= sqrtl((long double)(1.0) - vc_r*vc_r);
    fx_3l /= sqrtl((long double)(1.0) - vc_l*vc_l);
    fy_3l /= sqrtl((long double)(1.0) - vc_l*vc_l);
#else
    // или основываясь на соответсвующей проекции скорости
    double vxc_l = vx_l/c;
    double vyc_l = vy_l/c;

    double vxc_r = vx_r/c;
    double vyc_r = vy_r/c;

    fx_3r /= sqrtl((long double)(1.0) - vxc_r*vxc_r);
    fy_3r /= sqrtl((long double)(1.0) - vyc_r*vyc_r);
    fx_3l /= sqrtl((long double)(1.0) - vxc_l*vxc_l);
    fy_3l /= sqrtl((long double)(1.0) - vyc_l*vyc_l);
#endif

#ifndef USE_LEFT_CHARGE_ONLY
            // поле создаваемое правым вращающимся зарядом
            // в области левого вращающегося заряда
            if (0 != electr_magnet(&pendulum_right, driver_right,
                  X_l, Y_l, Z_l, t_i, q_right,
                  &E_x, &E_y, &E_z, &B_x, &B_y, &B_z, &rlagerror))
            {
                return -1;
            }

            *sum_rlagerror_sqare += Sq(rlagerror);

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

                //current_angle_l = Omega_l * t_i + Alpha_l;

                // Расчёт углового усилия действующего на левый заряд
                // со стороны поля правого заряда
                f_alpha_l  = fy_l * cosl(current_angle_l) - fx_l * sinl(current_angle_l);
                *F_alpha_l += f_alpha_l;
            }
#ifndef USE_LEFT_CHARGE_ONLY
            // поле создаваемое левым вращающимся зарядом
            // в области правого вращающегося заряда
            if (0 != electr_magnet(&pendulum_left, driver_left,
                  X_r, Y_r, Z_r, t_i, q_left,
                  &E_x, &E_y, &E_z, &B_x, &B_y, &B_z, &rlagerror))
            {
                return -1;
            }

            *sum_rlagerror_sqare += Sq(rlagerror);

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
                //current_angle_r = Omega_r * t_i + Alpha_r;

                // Расчёт углового усилия действующего на правый заряд
                // со стороны поля левого заряда
                f_alpha_r  = fy_r * cosl(current_angle_r) - fx_r * sinl(current_angle_r);
                *F_alpha_r += f_alpha_r;
            }
#endif
            if (to_log){
                printf("fx_l=%f fy_l=%f fx_r=%f fy_r=%f",fx_l, fy_l, fx_r, fy_r);
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
                                 long double * sum_rlagerror_sqare
                               )
{
    long double sign_a;
    long double sign_q;

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

    field E_x, E_y, E_z, B_x, B_y, B_z, A_x, A_y, A_z;

    coordinate rlagerror;
    *sum_rlagerror_sqare = 0.0;

    //i_l = i_a
    //i_r = i_q

    double q_left = -1.0;
    double q_right = +1.0;

    sign_a = q_left;
    sign_q = q_right;


            //printf("sign_a = %Lf sign_q = %Lf\n", sign_a, sign_q);
            //printf("Alpha_l = %Lf Alpha_r = %Lf\n", Alpha_l, Alpha_r);

#ifndef USE_LEFT_CHARGE_ONLY
            // поле в точке наблюдения создаваемое правым зарядом
            if (0 != electr_magnet_ex(&pendulum_right, driver_right,
                  X_a, Y_a, Z_a, t_i,
                  sign_q,
                  &E_x, &E_y, &E_z,
                  &B_x, &B_y, &B_z,
                  &A_x, &A_y, &A_z,
                  &rlagerror))
            {
                return -1;
            }

            *sum_rlagerror_sqare += Sq(rlagerror);

            Ex += E_x;
            Ey += E_y;
            Ez += E_z;

            Hx += B_x;
            Hy += B_y;
            Hz += B_z;

            Ax += A_x;
            Ay += A_y;
            Az += A_z;
#endif

            // поле в точке наблюдения создаваемое левым зарядом
            if (0 != electr_magnet_ex(&pendulum_left, driver_left,
                  X_a, Y_a, Z_a, t_i,
                  sign_a,
                  &E_x, &E_y, &E_z,
                  &B_x, &B_y, &B_z,
                  &A_x, &A_y, &A_z,
                  &rlagerror))
            {
                return -1;
            }

            *sum_rlagerror_sqare += Sq(rlagerror);

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
    long double * Txn, long double * Tyn, long double * Tzn,
    long double * Nx, long double * Ny, long double * Nz,
    long double * Sn,
    long double * En,
    long double * Hn,
    long double * An,
    long double * sum_rlagerror_sqare)
{
    return ccalc_Maxwells_stress_tensor(
        r*sinl(theta)*sinl(varphi), // X_a - в декартовой -> y в сферической
        r*cosl(theta),              // Y_a - в декартовой -> z в сферической
        r*sinl(theta)*cosl(varphi), // Z_a - в декартовой -> x в сферической
        t,
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
        sum_rlagerror_sqare);
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
    long double * sum_rlagerror_sqare)
{
    return ccalc_Maxwells_stress_tensor(
        r*cosl(theta),              // X_a - в декартовой -> z в сферической
        r*sinl(theta)*cosl(varphi), // Y_a - в декартовой -> x в сферической
        r*sinl(theta)*sinl(varphi), // Z_a - в декартовой -> y в сферической
        t,
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
        sum_rlagerror_sqare);
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
    long double * sum_rlagerror_sqare)
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
                                                     &Txn, &Tyn, &Tzn,
                                                     &Nx, &Ny, &Nz,
                                                     &Sn,
                                                     &En,
                                                     &Hn,
                                                     &An,
                                                     sum_rlagerror_sqare);
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

int spherical_x_ccalc_Maxwells_stress_tensor_R_t(
    long double theta, long double varphi, long double t,
    long double * pTxn, long double * pTyn, long double * pTzn,
    long double * pNx, long double * pNy, long double * pNz,
    long double * pSn,
    long double * pEn,
    long double * pHn,
    long double * pAn,
    long double * sum_rlagerror_sqare)
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

    int ret = spherical_x_ccalc_Maxwells_stress_tensor(sphere_R, theta, varphi, t,
                                                     &Txn, &Tyn, &Tzn,
                                                     &Nx, &Ny, &Nz,
                                                     &Sn,
                                                     &En,
                                                     &Hn,
                                                     &An,
                                                    sum_rlagerror_sqare);
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