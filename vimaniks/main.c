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

int main()
{
	long double xcl = cget_xc_l();
    printf("xcl = %Lf\n", xcl);
    
    cset_c(1.0);
    cset_timespan_Epsilon(1.e-15);
    cset_vc(0.8);
    
    long double omega = cget_omega();
    long double T = (2*M_PI)/cget_omega(); // период вращения
    int time_steps_number = 36000;                      // разбиваем период на шаги
    long double dt = T / time_steps_number;             // длительность шага
    
    int n = 1;
    long double t_i = T/3;
    
    long double Alpha0_l = 0;
    long double Alpha0_r = 0;
    int To_log = 0;
    

    cset_max_steps(50);

    long double Fy = ccalc_sum_Fy_t(n, t_i, Alpha0_l, Alpha0_r, To_log);

    long double R;
    // радиус сферы интегрирования
    R = 4 * cget_R_l() + 2 * cget_S();
    R *= 1.5;

    long double theta = M_PI / 2;
    long double varphi = 0;


    long double p = spherical_ccalc_Maxwells_stress_tensor(
        R, theta, varphi, t_i);
    printf("p = %0.36Le Fy = %0.36Lf\n", p, Fy);

    
    return 0;
}