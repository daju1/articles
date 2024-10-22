#include <stdio.h>
#include "pendulum_tensor.h"

extern gsl_odeiv2_driver * driver_left;
extern gsl_odeiv2_driver * driver_right;

extern struct pendulum pendulum_left;
extern struct pendulum pendulum_right;

int main (void)
{
    double m = 1.0;
    double g = 0.1;
    double R = 1.0;
    double p0 = 0.05;
    double q0 = 0.0;
    double t0 = 0.0;
    double xc = 0.0;
    double yc = 0.0;

    cset_timespan_Epsilon(1.e-15);
    cset_distance_Epsilon(1.e-8);
    cset_min_newton_step(0.1);
    cset_newton_step_multiplier(0.9999);

    init_left (m, -g, R, -p0, q0+M_PI, t0, xc, yc);
    init_right (m, g, R, p0, q0, t0, xc, yc);
    alloc_left();
    alloc_right();

    double init_T;
    double max_ti = 50;
    double dti = 0.01;

    calc_pendulum_period(&pendulum_right, driver_right, &init_T, max_ti, dti);

    timevalue T = init_T;
    double fT;
    find_period_by_newton_root(&pendulum_right, driver_right, &T, &fT);

    printf("period T = %f\n", T);
    printf("fT = %0.30e\n", fT);


    int To_log = 0;

    double sum_rlagerror_sqare;

    force Fx;
    force Fy;
    force Fz;
    force F_alpha_l;
    force F_alpha_r;

    timevalue t_i = T/3;
    if (0 != ccalc_sum_F_t(t_i, &Fx, &Fy, &Fz, &F_alpha_l, &F_alpha_r, &sum_rlagerror_sqare, To_log))
    {
        printf("ccalc_sum_Fy_t error\n");
    }

    printf("Fx=%f, Fy=%f, Fz=%f, F_alpha_l=%f, F_alpha_r=%f\n\n", Fx, Fy, Fz, F_alpha_l, F_alpha_r);

    return 0;

    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double t = 7.5;
    double t2_left, t2_right;
    double k_left, k_right, rlagerror;

    double x0 = -7.0;
    double dx = 0.1;
    double y0 = -7.0;
    double dy = 0.1;

    int ret;

    for (int i = 0; i < 140; ++i){
        for (int j = 0; j < 140; ++j){
            x = x0 + i * dx;
            y = y0 + j * dy;
            // отношение радиуса Лиенара Вихерта к длине радиус-вектора
            ret = klw(&pendulum_left, driver_left, x, y, z, t, &t2_left,
                    &k_left, &rlagerror);
            ret = klw(&pendulum_right, driver_right, x, y, z, t, &t2_right,
                    &k_right, &rlagerror);

            printf("x = %f y = %f t2 =%f %f k = %0.12f %0.12f rlagerror = %0.12e\n",
                x, y, t2_left, t2_right, k_left, k_right, rlagerror);
        }
    }

    return 0;
    
}