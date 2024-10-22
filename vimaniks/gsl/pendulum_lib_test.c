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
    {
    double f;
    int ret;

    double Momenta, q;

    double dot_q;   // dphi/dt
    double ddot_q;  // d2phi/dt2
    double dddot_q; // d3phi/dt3

    coordinate sx, sy;
    velocity vx, vy;

    double wx;
    double wy;

    double dot_wx;
    double dot_wy;

    double pre_ti = 0;
    double pre_sy = 0;
    int start = 1;
    int null_index = 0;

    for (double ti = 0; ti < 50; ti += 0.01)
    {
        apply(ti, NULL, &Momenta, &q,
              &dot_q, &ddot_q, &dddot_q,
              &sx, &sy,
              &vx, &vy,
              &wx, &wy,
              &dot_wx, &dot_wy
        );
        //printf("ti=%f sy=%f\n", ti, sy);

        if (pre_sy * sy < 0.0) {
            if (1 == null_index)
            {
                init_T = (pre_ti+ti)/2;
                printf("init_T=%f\n", init_T);
            }

            printf("(pre_sy=%f,sy=%f)\n", pre_sy,sy);
            printf("(pre_ti=%f,ti=%f)\n", pre_ti,ti);

            null_index += 1;
        }
        if (pre_sy * sy != 0 || start)
        {
            pre_ti = ti;
            pre_sy = sy;
            if (sy != 0.0){
                start = 0;
            }
        }       
    }
    }


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