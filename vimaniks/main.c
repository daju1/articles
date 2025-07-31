#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include "lw.h"
#include "lw_rotate.h"

#include "stdlib.h"

#define _USE_MATH_DEFINES
#include <math.h>
#define Pi M_PI


#include "lw.h"
#include "lw_rotate.h"
#include "lw_vimanic.h"
#include "lw_tensor.h"

extern int logging;

int main()
{
    long double xcl = cget_xc_l();
    printf("xcl = %Lf\n", xcl);
    printf("sizeof(double) = %ld\n", sizeof(double));
    printf("sizeof(long double) = %ld\n", sizeof(long double));
    
    cset_c(1.0);
    cset_timespan_Epsilon(1.e-15);
    cset_distance_Epsilon(1.e-8);
    cset_vc(0.8);
    cset_no_retardation_test(0);

    long double omega = cget_omega();
    long double T = (2*M_PI)/cget_omega(); // период вращения
    int time_steps_number = 36000;                      // разбиваем период на шаги
    long double dt = T / time_steps_number;             // длительность шага

    int n = 1;
    long double t_i = T/3;

    long double Alpha0_l = 0;
    long double Alpha0_r = 0;
    int To_log = 0;

    long double sum_rlagerror_sqare;

    cset_max_steps(50);

    long double Fx;
    long double Fy;
    long double Fz;
    long double F_alpha_l;
    long double F_alpha_r;
    if (0 != ccalc_sum_F_t(n, t_i, Alpha0_l, Alpha0_r, &Fx, &Fy, &Fz, &F_alpha_l, &F_alpha_r, &sum_rlagerror_sqare, To_log))
    {
        printf("ccalc_sum_Fy_t error\n");
    }

    long double R;
    // радиус сферы интегрирования
    R = 4 * cget_R_l() + 2 * cget_S();
    R *= 1.5;

    long double theta = M_PI / 4;
    long double varphi = 0;

    long double Txn;
    long double Tyn;
    long double Tzn;
    long double Nx;
    long double Ny;
    long double Nz;
    long double Sn;
    long double E1n;
    long double E2n;
    long double En;
    long double Hn;
    long double An;
    long double jn;
    int ret = spherical_y_ccalc_Maxwells_stress_tensor(R, theta, varphi, t_i,
                                                     &Txn, &Tyn, &Tzn,
                                                     &Nx, &Ny, &Nz,
                                                     &Sn,
                                                     &E1n,
                                                     &E2n,
                                                     &En,
                                                     &Hn,
                                                     &An,
                                                     &jn,
                                                     &sum_rlagerror_sqare);
    if (0 != ret)
    {
        printf("spherical_y_ccalc_Maxwells_stress_tensor error\n");
    }
    printf("Txn = %0.36Le Tyn = %0.36Lf Tzn = %0.36Lf\n", Txn, Tyn, Tzn);
    printf("Nx = %0.36Le Ny = %0.36Lf Nz = %0.36Lf\n", Nx, Ny, Nz);
    printf("Sn = %0.36Lf Fy = %0.36Lf\n", Sn, Fy);
    printf("En = %0.36Lf jn = %0.36Lf\n", En, jn);
    printf("E1n = %0.36Lf\n", E1n);
    printf("E2n = %0.36Lf\n", E2n);
    printf("E1n+E2n = %0.36Lf En = %0.36Lf\n", E1n+E2n, En);

    ret = spherical_x_ccalc_Maxwells_stress_tensor(R, theta, varphi, t_i,
                                                     &Txn, &Tyn, &Tzn,
                                                     &Nx, &Ny, &Nz,
                                                     &Sn,
                                                     &E1n,
                                                     &E2n,
                                                     &En,
                                                     &Hn,
                                                     &An,
                                                     &jn,
                                                     &sum_rlagerror_sqare);
    if (0 != ret)
    {
        printf("spherical_x_ccalc_Maxwells_stress_tensor error\n");
    }
    printf("Txn = %0.36Le Tyn = %0.36Lf Tzn = %0.36Lf\n", Txn, Tyn, Tzn);
    printf("Nx = %0.36Le Ny = %0.36Lf Nz = %0.36Lf\n", Nx, Ny, Nz);
    printf("Sn = %0.36Lf Fy = %0.36Lf\n", Sn, Fy);
    printf("En = %0.36Lf jn = %0.36Lf\n", En, jn);
    printf("E1n = %0.36Lf\n", E1n);
    printf("E2n = %0.36Lf\n", E2n);
    printf("E1n+E2n = %0.36Lf En = %0.36Lf\n", E1n+E2n, En);

    return 0;
}
