#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include "lw.h"
#include "lw_rotate.h"
#include "lw_vimanic.h"

double sum_Fy_Integrand(double t, void *user_data)
{
    long double res = ccalc_sum_Fy_t(1, (long double) t,
                  (long double) 0.0,
                  (long double) 0.0, 0);

    return (double)res;
}

// scipy.LowLevelCallable interface
double Maxwells_stress_tensor_Integrand(int n, double *xx, void *user_data){
    // theta, varphi, t

    #define theta  xx[0]
    #define varphi xx[1]
    #define t      xx[2]

    long double res = spherical_ccalc_Maxwells_stress_tensor_R_t((long double)theta, (long double)varphi, (long double)t);

    return (double)res;
}
