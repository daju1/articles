#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define Pi M_PI
#include "tzap.h"
#include <assert.h>
//#define USE_DEBUG
#include "dbg_info.h"
#include "integrate.h"
#include "calc.h"

extern velocity g_c;

int test_v0()
{
	int error = 0;
	coordinate r_zap;
	distance R_zap, R_lw_zap;
	potential phi_lw, A_lw;
	field E_minus_grad_phi_R0, E_minus_1_c_dA_dt_R0;

	/* Текущий момент */
	timevalue t = 5;
	/* расстояние от центра сферы к точке наблюдения. Точка наблюдения расположена на оси z в сферической системе координат */
	coordinate R0 = 2;
	/* начальный радиус заряженной сферы в момент t=t_start */
	coordinate r0 = 1;
	/* скорость в момент t_start*/
	velocity v0 = 0.0;
	/* ускорение */
	acceleration a0 = 0.0;
	/* минимально возможный радиус заряженной сферы (из соображений упругости) в момент более ранний чем t=t_start */
	coordinate r_min = 0.1;
	/* угловая координата заряда на заряженной сфере в сферической системе координат */
	angle theta = Pi/2;
	/* Заряд сферы */
	charge q = 1.0;
	timevalue t_zap;
	coordinate r;

	error = integral_phi(/*q*/-1, /*t*/0, /*R0*/2.0, /*r0*/2.0, /*v0*/g_c/3.0, /*a0*/0.0, /*r_min*/0.1, &phi_lw);
	printf("\nphi_lw = %0.10Lf error = %d r = %Lf\n", phi_lw, error, r);
	#ifdef CALC_LW_WITHOUT_LAGGING
	printf ("should be -0.4315231087\n\n");
	#else
	printf ("should be -0.5198603854\n\n");
	#endif

	error = integral_phi(/*q*/-1, /*t*/0, /*R0*/2.0, /*r0*/2.0, /*v0*/0.0, /*a0*/0.0, /*r_min*/0.1, &phi_lw);
	printf("\nphi_lw = %0.10Lf error = %d r = %Lf\n", phi_lw, error, r);
	#ifdef CALC_LW_WITHOUT_LAGGING
	printf ("should be -0.5\n\n");
	#else
	printf ("should be -0.5\n\n");
	#endif

	error = calc_tzap(q, /*t*/0.0, /*R0*/0.0, /*r0*/2.0, /*v0*/1.0, /*a0*/0.0, /*theta*/0.0, r_min, &t_zap, &r_zap);
	printf("t_zap = %Lf\n", t_zap);

	error = integral_phi(/*q*/-1, /*t*/0, /*R0*/0.0, /*r0*/2.0, /*v0*/g_c/3.0, /*a0*/0.0, /*r_min*/0.1, &phi_lw);
	printf("\nphi_lw = %0.10Lf error = %d r = %Lf\n", phi_lw, error, r);
	#ifdef CALC_LW_WITHOUT_LAGGING
	printf ("should be -0.375\n\n");
	#else
	printf ("should be -0.5\n\n");
	#endif

	error = integral_phi(/*q*/-1, /*t*/0, /*R0*/1.0, /*r0*/2.0, /*v0*/g_c/3.0, /*a0*/0.0, /*r_min*/0.1, &phi_lw);
	printf("\nphi_lw = %0.10Lf error = %d r = %Lf\n", phi_lw, error, r);
	#ifdef CALC_LW_WITHOUT_LAGGING
	printf ("should be -0.38348\n\n");
	#else
	printf ("should be -0.5047083549\n\n");
	#endif

	error = integral_phi(/*q*/-1, /*t*/0, /*R0*/3.0, /*r0*/2.0, /*v0*/g_c/3.0, /*a0*/0.0, /*r_min*/0.1, &phi_lw);
	printf("\nphi_lw = %0.10Lf error = %d r = %Lf\n", phi_lw, error, r);
	#ifdef CALC_LW_WITHOUT_LAGGING
	printf ("should be -0.31720\n\n");
	#else
	printf ("should be -0.3465735903\n\n");
	#endif


	error = integral_phi(/*q*/-1, /*t*/0, /*R0*/0.0, /*r0*/2.0, /*v0*/1.0*g_c/3.0, /*a0*/0.0, /*r_min*/0.1, &phi_lw);
	printf("\nphi_lw = %0.10Lf error = %d r = %Lf\n", phi_lw, error, r);
#ifndef CALC_LW_WITHOUT_LAGGING
	printf("should be -0.5\n\n");
#endif

	error = integral_phi(/*q*/-1, /*t*/0, /*R0*/0.0, /*r0*/2.0*g_c/3.0, /*v0*/0.0, /*a0*/0.1*g_c/3.0, /*r_min*/0.1, &phi_lw);
	printf("\nphi_lw = %0.10Lf error = %d r = %Lf\n", phi_lw, error, r);
	#ifndef CALC_LW_WITHOUT_LAGGING
	printf ("should be -0.50576 tzap should be -0.67424\n\n");
	#endif

	error = integral_phi(/*q*/-1, /*t*/0, /*R0*/0.0, /*r0*/1.0*g_c / 3.0, /*v0*/0.0, /*a0*/0.1*g_c / 3.0, /*r_min*/0.1, &phi_lw);
	printf("\nphi_lw = %0.10Lf error = %d r = %Lf\n", phi_lw, error, r);
#ifndef CALC_LW_WITHOUT_LAGGING
	printf("should be -1.0057 tzap should be -0.33521\n\n");
#endif

	error = integral_phi(/*q*/-1, /*t*/0, /*R0*/1.0*g_c/3.0, /*r0*/1.0*g_c / 3.0, /*v0*/0.0, /*a0*/0.1*g_c / 3.0, /*r_min*/0.1, &phi_lw);
	printf("\nphi_lw = %0.10Lf error = %d r = %Lf\n", phi_lw, error, r);
#ifndef CALC_LW_WITHOUT_LAGGING
	printf("should be -0.98892 tzap should be 2.*sqrt(445.+5.*cos(theta)-5.*sqrt(cos(theta)^2+180.*cos(theta)+7919.))\n\n");
#endif

	error = integral_phi(/*q*/-1, /*t*/0, /*R0*/2.0*g_c / 3.0, /*r0*/1.0*g_c / 3.0, /*v0*/0.0, /*a0*/0.1*g_c / 3.0, /*r_min*/0.1, &phi_lw);
	printf("\nphi_lw = %0.10Lf error = %d r = %Lf\n", phi_lw, error, r);
#ifndef CALC_LW_WITHOUT_LAGGING
	printf("should be -0.49443 tzap should be 2.*sqrt(445.+10.*cos(theta)-10.*sqrt(cos(theta)^2+90.*cos(theta)+1979.))\n\n");
#endif


return 0;
	/*
        cdef res = integral_phi(q, t, R0, r0, v0, a0)
        if res < -500:
            print (res, q, t, R0, r0, v0, a0)

	(-6749.344473145355, -1, 5.5, 2.8421709430404007e-13, 2, 0, -0.150000000000000)
	(-2674.2184374462895, -1, 6.0, 2.8421709430404007e-13, 2, 0, -0.150000000000000)
	(-1650.9694948050178, -1, 6.5, 2.8421709430404007e-13, 2, 0, -0.150000000000000)
	(-1186.1298919436645, -1, 7.0, 2.8421709430404007e-13, 2, 0, -0.150000000000000)
	(-921.0452564374108, -1, 7.5, 2.8421709430404007e-13, 2, 0, -0.150000000000000)
	*/
	error = integral_phi(-1, 5.5, 2.8421709430404007e-13, 2, 0, -0.150000000000000, r_min, &phi_lw);
	printf("phi_lw = %Lf error = %d\n", phi_lw, error);

	for (long double t = 0; t < 8.0; t += 0.1)
	{
		error = integral_phi(-1, t, 2.8421709430404007e-13, 2, 0, -0.150000000000000, r_min, &phi_lw);
		printf("t = %Lf phi_lw = %Lf error = %d\n", t, phi_lw, error);
	}

	error = integral_phi(-1, 6.0, 2.8421709430404007e-13, 2, 0, -0.150000000000000, r_min, &phi_lw);
	printf("phi_lw = %Lf\n", phi_lw);
	error = integral_phi(-1, 6.5, 2.8421709430404007e-13, 2, 0, -0.150000000000000, r_min, &phi_lw);
	printf("phi_lw = %Lf\n", phi_lw);
	error = integral_phi(-1, 7.0, 2.8421709430404007e-13, 2, 0, -0.150000000000000, r_min, &phi_lw);
	printf("phi_lw = %Lf\n", phi_lw);
	error = integral_phi(-1, 7.5, 2.8421709430404007e-13, 2, 0, -0.150000000000000, r_min, &phi_lw);
	printf("phi_lw = %Lf\n", phi_lw);

	get_R(1.5, 1.499999999999999, 0.0);
	get_R(1.5, 1.4999999979313956, 0.0);

	error = calc_R_lw(q, t, R0, r0, v0, a0, theta, &t_zap, &r_zap, &R_zap, r_min, &R_lw_zap);

	/* скалярный потенциал Лиенара Вихерта зарядов равномерно распределённых по сферической поверхности радиуса r0 и движущихся из центра.  */
	/* varphi := proc (q, t, r__0, v__0, a__0, R__0) options operator, arrow;
	int(2*Pi*r(t, r__0, v__0, a__0)^2*sin(theta)*sigma(q, r__0)/K__zap(tzap(t, r__0, v__0, a__0, R__0, theta), r__0, v__0, a__0, R__0, theta), theta = 0 .. Pi) end proc;*/
	error = integral_phi(q, t, R0, r0, v0, a0, r_min, &phi_lw);
	printf("phi_lw = %Lf\n", phi_lw);
	error = integral_phi_and_E(q, t, R0, r0, v0, a0, &E_minus_grad_phi_R0, &E_minus_1_c_dA_dt_R0, r_min, &phi_lw, &A_lw);
	printf("phi_lw = %Lf E1=%Lf E2 = %Lf\n", phi_lw, E_minus_grad_phi_R0, E_minus_1_c_dA_dt_R0);

	/*infinity error result*/
	error = integral_phi(q, t, -1.0000000000000142, 1, 0, 0, r_min, &phi_lw);
	printf("phi_lw = %Lf\n", phi_lw);
	phi_lw = integral_phi_and_E(q, t, -1.0000000000000142, 1, 0, 0, &E_minus_grad_phi_R0, &E_minus_1_c_dA_dt_R0, r_min, &phi_lw, &A_lw);
	printf("phi_lw = %Lf E1=%Lf E2 = %Lf\n", phi_lw, E_minus_grad_phi_R0, E_minus_1_c_dA_dt_R0);

	/*hung
	calc_tzap(t=5.000000, R0=-1.500000, r0=2.000000, v0=0.0, a0=1.000000, theta=0.000000
	*/
	error = calc_tzap(q, /*t*/5.0, /*R0*/-1.5, /*r0*/2.0, /*v0*/0.0, /*a0*/1.0, /*theta*/0.0, r_min, &t_zap, &r_zap);
	printf("t_zap = %Lf\n", t_zap);

	/*hung

	integral_phi(q=-1.000000 t=5.000000, R0=-1.500000, r0=2.000000, v0=0.0, a0=1.000000)

	t2=1.000667 t1=0.999333 t=5.000000 R=3.999333 dR=1.333710e-03 dR_pre=-1.333710e-03
	t2=0.999333 t1=1.000667 t=5.000000 R=4.000667 dR=-1.333710e-03 dR_pre=1.333710e-03
	t2=1.000667 t1=0.999333 t=5.000000 R=3.999333 dR=1.333710e-03 dR_pre=-1.333710e-03
	t2=0.999333 t1=1.000667 t=5.000000 R=4.000667 dR=-1.333710e-03 dR_pre=1.333710e-03
	*/

	error = integral_phi(-q, 5.0, -1.5, 2.0, 0.0, 1.0, r_min, &phi_lw);
	printf("phi_lw = %Lf\n", phi_lw);
	error = integral_phi_and_E(-q, 5.0, -1.5, 2.0, 0.0, 1.0, &E_minus_grad_phi_R0, &E_minus_1_c_dA_dt_R0, r_min, &phi_lw, &A_lw);
	printf("phi_lw = %Lf E1=%Lf E2 = %Lf\n", phi_lw, E_minus_grad_phi_R0, E_minus_1_c_dA_dt_R0);

	{
		field E1_p, E1_n, E2_p, E2_n;
		potential phi_p, phi_n;
		potential A_p, A_n;
		velocity v0_p = 0.0;
		velocity v0_n = 0.0;
		acceleration a0_p = 0.0001;
		acceleration a0_n = 0.01;
		field E_p, E_n, E1, E2, E;
		charge q_add;

		for (double ti = 0.0; ti <= 15.0; ti += 0.1)
		//for (double R0_i = r0+0.5; R0_i < 20.0; R0_i += 0.5)
		{
			for (double R0_i = r0+0.5; R0_i < 20.0; R0_i += 0.5)
			//for (double ti = 0.0; ti <= 15.0; ti += 0.1)
			{
				printf("ti=%03.1f R0_i=%03.1f ", ti, R0_i);

				error = integral_phi_and_E(+q, ti, R0_i, r0, v0_p, a0_p, &E1_p, &E2_p, r_min, &phi_p, &A_p);
				error = integral_phi_and_E(-q, ti, R0_i, r0, v0_n, a0_n, &E1_n, &E2_n, r_min, &phi_n, &A_n);
				//printf("phi_p = %Lf ", phi_p);
				//printf("phi_n = %Lf ", phi_n);

				printf("E1_p=%Le E2_p=%Le ", E1_p, E2_p);
				printf("E1_n=%Le E2_n=%Le ", E1_n, E2_n);

				E_p = E1_p + E2_p;
				E_n = E1_n + E2_n;

				printf("E_p=%Le ", E_p);
				printf("E_n=%Le ", E_n);

				E1 = E1_p + E1_n;
				E2 = E2_p + E2_n;

				printf("E1=%Le ", E1);
				printf("E2=%Le ", E2);

				E = E_p + E_n;

				printf("E=%Lf\n", E);

				q_add = E * 4 * Pi*R0_i*R0_i;
				printf("q_add=%Lf\n", q_add);
			}
			printf("\n");
		}
		error = integral_phi_and_E(q, t, R0, r0, 0.0, 0.0, &E_minus_grad_phi_R0, &E_minus_1_c_dA_dt_R0, r_min, &phi_lw, &A_lw);
		printf("phi_lw = %Lf E1=%Lf E2 = %Lf\n", phi_lw, E_minus_grad_phi_R0, E_minus_1_c_dA_dt_R0);

		error = integral_phi_and_E(-q, t, R0, r0, 0.0, 0.0, &E_minus_grad_phi_R0, &E_minus_1_c_dA_dt_R0, r_min, &phi_lw, &A_lw);
		printf("phi_lw = %Lf E1=%Lf E2 = %Lf\n", phi_lw, E_minus_grad_phi_R0, E_minus_1_c_dA_dt_R0);

	}

	return 0;
}
