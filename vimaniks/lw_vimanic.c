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

#define Sq(x) (x)*(x)

static long double R_r = 1;
static long double R_l = 1;
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

long double xc_l() { return -S/2 - R_l; }
long double xc_r() { return  S/2 + R_r; }

long double yc_l() { return 0; }
long double yc_r() { return 0; }

long double zc_l() { return 0; }
long double zc_r() { return 0; }

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

void ccalc_sum_Fy_t(int N, long double t_i,
                  long double alpha0_l,
                  long double alpha0_r,
                  long double * Fy,
                  long double * F_alpha_l,
                  long double * F_alpha_r,
                  _Bool to_log)
{
    long double omega = cget_omega();

    long double sign_a;
    long double sign_q;
    
    long double Alpha_l;
    long double Alpha_r;
    
    long double current_angle_l;
    long double current_angle_r;
    
    long double Omega_l = + omega;
    long double Omega_r = - omega;
    
    long double Xa, Ya, Za;
    long double Xq, Yq, Zq;

    //# current positions of rotated masses
    /*
    sign_r = []
    sign_l = []
    alpha_r = []
    alpha_l = []

    # n - number of charges per circle
    for i in range(n):
        sign_r += [lambda i=i : -((i%2)*2-1)]
        sign_l += [lambda i=i :   (i%2)*2-1]

        alpha_r += [lambda i=i : i * np.float128(2*np.pi)/n + alpha0_r]
        alpha_l += [lambda i=i : i * np.float128(2*np.pi)/n + alpha0_l]
        
        */

    long double Fx_l = 0.0;
    long double Fy_l = 0.0;
    long double Fx_r = 0.0;
    long double Fy_r = 0.0;

    if (F_alpha_l) {
        *F_alpha_l = 0.0;
    }

    if (F_alpha_r) {
        *F_alpha_r = 0.0;
    }
    
    long double fx_l, fx_r, fy_l, fy_r;
    long double f_alpha_l, f_alpha_r;
    
    long double E_x, E_y, E_z, B_x, B_y, B_z;
    
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
            
            
            // поле создаваемое правым вращающимся зарядом в области левого вращающегося заряда
            //(E_x, E_y, E_z, B_x, B_y, B_z) = EB_lw(Xa, Ya, Za, t_i, sign_q, xc_r, yc_r, zc_r, R_r, Omega_r, Alpha_r)
                
            electr_magnet(Xa, Ya, Za, t_i,
                  sx, sy, sz, vx, vy, vz, wx, wy, wz,
                  sign_q,
                  &E_x, &E_y, &E_z, &B_x, &B_y, &B_z, xc_r(), yc_r(), zc_r(), R_r, Omega_r, Alpha_r);

            // сила действующая на левый заряд со стороны поля правого заряда
            fx_l = (E_x)*sign_a;
            fy_l = (E_y)*sign_a;

            Fx_l += fx_l;
            Fy_l += fy_l;
            
            if (F_alpha_l) {
            
                current_angle_l = Omega_l * t_i + Alpha_l;

                // Расчёт углового усилия действующего на левый заряд со стороны поля правого заряда
                f_alpha_l  = fy_l * cos(current_angle_l) - fx_l * sin(current_angle_l);
                *F_alpha_l += f_alpha_l;
            }
            
            // поле создаваемое левым вращающимся зарядом в области правого вращающегося заряда
            //(E_x, E_y, E_z, B_x, B_y, B_z) = EB_lw(Xq, Yq, Zq, t_i, sign_a, xc_l, yc_l, zc_l, R_l, Omega_l, Alpha_l)
                
            electr_magnet(Xq, Yq, Zq, t_i,
                  sx, sy, sz, vx, vy, vz, wx, wy, wz,
                  sign_a,
                  &E_x, &E_y, &E_z, &B_x, &B_y, &B_z, xc_l(), yc_l(), zc_l(), R_l, Omega_l, Alpha_l);

            // сила действующая на правый заряд со стороны поля левого заряда
            fx_r = (E_x)*sign_q;
            fy_r = (E_y)*sign_q;
            
            Fx_r += fx_r;
            Fy_r += fy_r;
            
            if (F_alpha_r) {
                current_angle_r = Omega_r * t_i + Alpha_r;

                // Расчёт углового усилия действующего на правый заряд со стороны поля левого заряда
                f_alpha_r  = fy_r * cos(current_angle_r) - fx_r * sin(current_angle_r);
                *F_alpha_r += f_alpha_r;
            }

            if (to_log){
                printf("fx_l=%Lf fy_l=%Lf fx_r=%Lf fy_r=%Lf",fx_l, fy_l, fx_r, fy_r);
            }
        }
    }

    if (Fy)
    {
        // Интегральная величина тяги в направлении оси y
        *Fy = Fy_l + Fy_r;
    }
}

long double ccalc_Maxwells_stress_tensor(long double X_a, long double Y_a, long double Z_a, long double t_i,
                                 int N,
                                long double cos_nx, long double cos_ny, long double cos_nz,
                                long double alpha0_l,
                                long double alpha0_r
                               )
{
    long double omega = cget_omega();
    
    long double sign_a;
    long double sign_q;
    
    long double Alpha_l;
    long double Alpha_r;
    
    long double Omega_l = + omega;
    long double Omega_r = - omega;
    
    long double Xa, Ya, Za;
    long double Xq, Yq, Zq;


    long double Ex = 0;
    long double Ey = 0;
    long double Ez = 0;
    
    long double Hx = 0;
    long double Hy = 0;
    long double Hz = 0;


    long double E_x, E_y, E_z, B_x, B_y, B_z;
    
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

            
            //(E_x, E_y, E_z, B_x, B_y, B_z) = EB_lw(X_a, Y_a, Z_a, t_i, sign_q, xc_r, yc_r, zc_r, R_r, Omega_r, Alpha_r)
            electr_magnet(X_a, Y_a, Z_a, t_i,
                  sx, sy, sz, vx, vy, vz, wx, wy, wz,
                  sign_q,
                  &E_x, &E_y, &E_z, &B_x, &B_y, &B_z, xc_r(), yc_r(), zc_r(), R_r, Omega_r, Alpha_r);

            Ex += E_x;
            Ey += E_y;
            Ez += E_z;
            
            Hx += B_x;
            Hy += B_y;
            Hz += B_z;
            

            //(E_x, E_y, E_z, B_x, B_y, B_z) = EB_lw(X_a, Y_a, Z_a, t_i, sign_a, xc_l, yc_l, zc_l, R_l, Omega_l, Alpha_l)
            electr_magnet(X_a, Y_a, Z_a, t_i,
                  sx, sy, sz, vx, vy, vz, wx, wy, wz,
                  sign_a,
                  &E_x, &E_y, &E_z, &B_x, &B_y, &B_z, xc_l(), yc_l(), zc_l(), R_l, Omega_l, Alpha_l);

            Ex += E_x;
            Ey += E_y;
            Ez += E_z;
            
            Hx += B_x;
            Hy += B_y;
            Hz += B_z;
        }
    }
            
    long double sigma_xx = ((long double)(1.0))/(4*M_PI)*( - Ex*Ex - Hx*Hx + Ey*Ey + Ez*Ez + Hy*Hy + Hz*Hz);
    long double sigma_yy = ((long double)(1.0))/(4*M_PI)*( - Ey*Ey - Hy*Hy + Ez*Ez + Ex*Ex + Hz*Hz + Hx*Hx);
    long double sigma_zz = ((long double)(1.0))/(4*M_PI)*( - Ez*Ez - Hz*Hz + Ex*Ex + Ey*Ey + Hx*Hx + Hy*Hy);
    
    long double sigma_xy = ((long double)(1.0))/(4*M_PI)*( - Ex*Ey - Hx*Hy );
    long double sigma_xz = ((long double)(1.0))/(4*M_PI)*( - Ex*Ez - Hx*Hz );
    long double sigma_yz = ((long double)(1.0))/(4*M_PI)*( - Ey*Ez - Hy*Hz );

    long double sigma_yx = ((long double)(1.0))/(4*M_PI)*( - Ey*Ex - Hy*Hx );
    long double sigma_zx = ((long double)(1.0))/(4*M_PI)*( - Ez*Ex - Hz*Hx );
    long double sigma_zy = ((long double)(1.0))/(4*M_PI)*( - Ez*Ey - Hz*Hy );
    
    //#T = [[sigma_xx, sigma_xy, sigma_xz],
    //#     [sigma_yx, sigma_yy, sigma_yz],
    //#     [sigma_zx, sigma_zy, sigma_zz]]
    
    //# Тамм параграф 33 формула (33.5)
    //# сила натяжения действующая на площадку поверхности интегрирования
    //# со стороны поля создаваемого вращающимися зарядами
    long double py = (sigma_yx * cos_nx + sigma_yy * cos_ny + sigma_yz * cos_nz);
    
    long double S_x = ((long double)(1.0))/(4*M_PI)*(Ey * Hz - Ez * Hy);
    long double S_y = ((long double)(1.0))/(4*M_PI)*(Ez * Hx - Ex * Hz);
    long double S_z = ((long double)(1.0))/(4*M_PI)*(Ex * Hy - Ey * Hx);

    //# Интегральная величина количества имульса электромагнитного поля, 
    //# вытекающего в единицу времени из замкнутого обьёма через площадку ЛЛ2 32.14
    return py;
}


// Интегрируем в сферической системе координат,
// у которой однако в соотвествие с принятыми в задаче
// наименованиями осей главная ось игрек вместо зет

// направление векторов нормали к сферической воображаемой поверхности инвертировано - снаружи вовнутрь
long double spherical_ccalc_Maxwells_stress_tensor(
    long double r, long double theta, long double varphi, long double t)
{
    return ccalc_Maxwells_stress_tensor(
        r*sin(theta)*cos(varphi),
        r*cos(theta),
        r*sin(theta)*sin(varphi),
        t,
        1,
        - sin(theta)*cos(varphi),
        - cos(theta),
        - sin(theta)*sin(varphi),
        0,
        0);
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

long double spherical_ccalc_Maxwells_stress_tensor_R_t(
    long double theta, long double varphi, long double t)
{
    return sphere_R * sphere_R * spherical_ccalc_Maxwells_stress_tensor(sphere_R, theta, varphi, t);
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
