#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include "lw.h"
#include "lw_rotate.h"
#include "lw_vimanic.h"
#include "lw_tensor.h"

double sum_Fy_Integrand(double t, void *user_data)
{
    int n = 1;

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
    if (0 != ccalc_sum_F_t(n, t, Alpha0_l, Alpha0_r, &Fx, &Fy, &Fz, &F_alpha_l, &F_alpha_r, &sum_rlagerror_sqare, To_log))
    {
        printf("ccalc_sum_Fy_t error\n");
    }

    return (double)Fy;
}

// scipy.LowLevelCallable interface
double Maxwells_stress_tensor_Integrand(int n, double *xx, void *user_data){
    // theta, varphi, t

    #define theta  xx[0]
    #define varphi xx[1]
    #define t      xx[2]

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

    long double sum_rlagerror_sqare;
    int ret = spherical_y_ccalc_Maxwells_stress_tensor_R_t(theta, varphi, (long double)t,
                                                     &Txn, &Tyn, &Tzn,
                                                     &Nx, &Ny, &Nz,
                                                     &Sn,
                                                     &En,
                                                     &Hn,
                                                     &An,
                                                     &sum_rlagerror_sqare);
    if (0 != ret)
    {
        printf("spherical_y_ccalc_Maxwells_stress_tensor error\n");
    }

    return (double)Tyn;
}
