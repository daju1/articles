#include <stdio.h>
#include "pendulum_lib.h"

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

    init (m, g, R, p0, q0, t0, xc, yc);
    alloc();

    double init_T;
    double max_ti = 50;
    double dti = 0.01;

    calc_pendulum_period(&init_T, max_ti, dti);

    timevalue T = init_T;
    double fT;
    find_period_by_newton_root(&T, &fT);

    printf("period T = %f\n", T);
    printf("fT = %0.30e\n", fT);

    return 0;


    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
    double t = 7.5;
    double t2;
    double k, rlagerror;

    double x0 = -7.0;
    double dx = 0.1;
    double y0 = -7.0;
    double dy = 0.1;

    for (int i = 0; i < 140; ++i){
        for (int j = 0; j < 140; ++j){
            x = x0 + i * dx;
            y = y0 + j * dy;
            // отношение радиуса Лиенара Вихерта к длине радиус-вектора
            int ret = klw(x, y, z, t, &t2,
                    &k, &rlagerror);

            printf("x = %f y = %f t2 =%f k = %0.12f rlagerror = %0.12e\n", x, y, t2, k, rlagerror);
        }
    }



    return 0;
    
}