#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define Pi M_PI
#include "tzap.h"
#include <assert.h>
//#define USE_DEBUG
#include "dbg_info.h"
#include "integrate.h"

#define DBG_INFO printf

extern velocity g_c;
extern timevalue * v_t;
extern long double multiplier_E;

long double get_sigma(charge q, coordinate r)
{
	long double sigma = q / (4*Pi*r*r);
	return sigma;
}

long double get_dS_dtheta(coordinate r, angle theta)
{
	long double dS_dtheta = 2*Pi*r*r*sin(theta);
	return dS_dtheta;
}

/* Радиус Лиенара Вихерта */
int calc_R_lw(charge q, timevalue t, coordinate R0, coordinate r0, velocity v0, acceleration a0, angle theta, timevalue * pt_zap, coordinate * pr_zap, distance * pR_zap, coordinate r_min, distance * R_lw_zap)
{
//#define DBG_INFO printf
	int err, error = 0;
	velocity v;
	distance R_zap;

	/* численный расчёта запаздывающего момента */
	coordinate r1_zap;
	err = calc_tzap(q, t, R0, r0, v0, a0, theta, r_min, pt_zap, &r1_zap);
	if (0 != err)
	{
		error += 1;
	}
	DBG_INFO("theta = %Lf t_zap = %Lf ", theta, *pt_zap);
	/* Запаздывающий радиус в зависимости от текущего момента */
#ifdef ALGORITHM_VERSION_0
	err = get_r(q, *pt_zap, r0, v0, a0, r_min, pr_zap); /* расстояние от заряда до центра сферы в запаздывающий момент времени */
#endif
#ifdef ALGORITHM_VERSION_1
	err = get_r_ex1(q, *pt_zap, r0, v0, r_min, pr_zap, 0);
#endif
#ifdef ALGORITHM_VERSION_2
	err = get_r_ex2(q, *pt_zap, r0, v0, r_min, pr_zap);
#endif
	if (0 != err)
	{
		error += 1;
	}
	assert(*pr_zap >= 0);
	*pR_zap = get_R(R0, *pr_zap, theta); /* расстояние от заряда до точки наблюдения в запаздывающий момент времени */
	DBG_INFO("r_zap = %Lf ", *pr_zap);
	DBG_INFO("R_zap = %Lf ", *pR_zap);
	R_zap = g_c * (t-*pt_zap);

	if (fabs(R_zap - *pR_zap) >= 1e-8)
	{
		printf("r_zap = %Lf\n", *pr_zap);
		printf("*pR_zap = %Lf\n", *pR_zap);
		printf("R_zap = %Lf\n", R_zap);
		printf("dR_zap = %e\n", fabs(R_zap - *pR_zap));
		//assert(0);
		//int * p = 0;
		//*p += 1;
	}

#ifdef ALGORITHM_VERSION_0
	v = get_v(q, *pt_zap, v0, a0);
#endif
#ifdef ALGORITHM_VERSION_1
	v = get_v_ex1(*pt_zap, v0, q);
#endif
#ifdef ALGORITHM_VERSION_2
	v = get_v_ex1(*pt_zap, v0, q);
#endif
	/* Радиус Лиенара Вихерта */
	*R_lw_zap = get_R(R0, *pr_zap, theta) - (v / g_c) * (R0 * cos(theta) - *pr_zap);
	DBG_INFO("R_lw_zap = %Lf ", *R_lw_zap);

//#define DBG_INFO
	return error;
}

/* Для расчёта радиальной компоненты электрического поля в точке наблюдения введём вспомогательную величину -
косинус угла между запаздывающим радиус-вектором (вектор из запаздывающего положения заряда в точку наблюдения)
и радиус-вектором из центра сферы в точку наблюдения
*/
/*
cos_alpha__zap := proc (t__zap, r__0, v__0, a__0, R__0, theta) options operator, arrow;
(R__0-r(t__zap, r__0, v__0, a__0)*cos(theta))/R__zap(t__zap, r__0, v__0, a__0, R__0, theta) end proc;
*/

long double get_cos_alpha_zap(coordinate R0, angle theta, coordinate r_zap, distance R_zap)
{
	long double cos_alpha_zap = (R0 - r_zap*cos(theta)) / R_zap;
	return cos_alpha_zap;
}

/* Скалярное произведение ускорения частицы в запаздывающий момент времени на запаздывающий радиус-вектор (вектор из запаздывающего положения заряда в точку наблюдения) */
/*
aR__zap := proc (t__zap, r__0, v__0, a__0) options operator, arrow;
a__r(t__zap, r__0, v__0, a__0)*(R__0*cos(theta)-r(t__zap, r__0, v__0, a__0)) end proc;
*/

long double get_aR_zap(coordinate R0, angle theta, coordinate r_zap, acceleration a_zap)
{
	long double aR_zap = a_zap*(R0*cos(theta) - r_zap);
	return aR_zap;
}

/* Первое слагаемое радиальной компоненты электрического поля - минус градиент скалярного потенциала */
/*
E_minus_grad_varphi__R__0 := proc (q, t__zap, r__0, v__0, a__0, R__0, theta) options operator, arrow;
sigma(q, r__0)*
(
R__zap(t__zap, r__0, v__0, a__0, R__0, theta)*cos_alpha__zap(t__zap, r__0, v__0, a__0, R__0, theta)*(1+aR__zap(t__zap, r__0, v__0, a__0)/c^2-v__r(t__zap, r__0, v__0, a__0)^2/c^2)
/
K__zap(t__zap, r__0, v__0, a__0, R__0, theta)
- v__r(t__zap, r__0, v__0, a__0)*cos(theta)/c
)
/
K__zap(t__zap, r__0, v__0, a__0, R__0, theta)^2
end proc;
*/

field get_E_minus_grad_phi_R0(angle theta, velocity v_zap, distance R_zap, long double aR_zap, distance R_lw_zap, long double cos_alpha_zap)
{
	field E_minus_grad_phi_R0 =
		(
			(cos_alpha_zap * R_zap / R_lw_zap) * (1.0 + (aR_zap - v_zap * v_zap) / (g_c * g_c) )
			- v_zap*cos(theta) / g_c
		)
		/ (R_lw_zap * R_lw_zap);
#ifndef USE_NORM
	DBG_INFO("E_minus_grad_phi_R0 = %0.25Lf\n", E_minus_grad_phi_R0);
	E_minus_grad_phi_R0 *= multiplier_E;
#endif
	DBG_INFO("E_minus_grad_phi_R0 = %0.25Le\n", E_minus_grad_phi_R0);
	return E_minus_grad_phi_R0;
}

/*double calc_E_minus_grad_varphi_R0(double t, double R0, double r0, double a0, double theta, double * pt_zap, double * pr_zap, double * pR_zap)
{

}*/


/* Второе слагаемое компоненты электрического поля */
/*
E_minus_1_c_dA_dt__R__0 := proc (q, t__zap, r__0, v__0, a__0, R__0, theta) options operator, arrow;
cos(theta)*sigma(q, r__0)*
(g_c
v__r(t__zap, r__0, v__0, a__0)*(R__zap(t__zap, r__0, v__0, a__0, R__0, theta)*(v__r(t__zap, r__0, v__0, a__0)^2/c-aR__zap(t__zap, r__0, v__0, a__0)/c-c)/K__zap(t__zap, r__0, v__0, a__0, R__0, theta)+c)/c^2
- a__r(t__zap, r__0, v__0, a__0)*R__zap(t__zap, r__0, v__0, a__0, R__0, theta)/c^2
)
/
K__zap(t__zap, r__0, v__0, a__0, R__0, theta)^2
end proc;
*/

field get_E_minus_1_c_dA_dt_R0(angle theta, velocity v_zap, acceleration a_zap, distance R_zap, long double aR_zap, distance R_lw_zap)
{
	field E_minus_1_c_dA_dt_R0 =
		cos(theta) *
		(
			(v_zap / (g_c * g_c)) * ( (R_zap / R_lw_zap) * ( (v_zap * v_zap - aR_zap) / g_c - g_c) + g_c)
			- a_zap * R_zap / (g_c * g_c)
		)
		/
		(R_lw_zap * R_lw_zap);
#ifndef USE_NORM
	DBG_INFO("E_minus_1_c_dA_dt_R0 = %0.25Lf\n", E_minus_1_c_dA_dt_R0);
	E_minus_1_c_dA_dt_R0 *= multiplier_E;
#endif
	DBG_INFO("E_minus_1_c_dA_dt_R0 = %0.25Lf\n", E_minus_1_c_dA_dt_R0);
	return E_minus_1_c_dA_dt_R0;
}


field get__E(angle theta, velocity v_zap, acceleration a_zap, distance R_zap, long double aR_zap, distance R_lw_zap, long double cos_alpha_zap)
{
	field E =
		(
			((cos_alpha_zap - v_zap * cos(theta) / g_c) * R_zap) * (1.0 + (aR_zap - v_zap * v_zap) / (g_c * g_c))
			- a_zap*cos(theta) * R_zap * R_lw_zap / (g_c * g_c)
		)
		/ (R_lw_zap * R_lw_zap * R_lw_zap);
#ifndef USE_NORM
	DBG_INFO("E = %0.25Lf\n", E);
	E *= multiplier_E;
#endif
	DBG_INFO("E = %0.25Le\n", E);
	return E;
}

field get_E(angle theta, timevalue t, timevalue t_zap, coordinate R0, coordinate r_zap, velocity v_zap, acceleration a_zap)
{
	//Радиус Лиенара Вихерта в сферической системе координат через запаздывающий момент вычисляется как

	//	\[{ {R}^ {*}} = c\left(t - t' \right)-\frac{v}{c}\left( {{R}_{0}}\cos \theta -r \right)\]
	distance R_lw_zap = g_c * (t - t_zap) - (v_zap / g_c) * (R0 * cos(theta) - r_zap);
	distance R_lw_zap_2 = R_lw_zap * R_lw_zap;
	distance R_lw_zap_3 = R_lw_zap_2 * R_lw_zap;
	DBG_INFO("R_lw_zap = %0.25Lf\n", R_lw_zap);


	// $$\frac{dE}{dq}=\frac{\left( {{R}_{0}}-\left( r+\left( t-{t}' \right)v \right)\cos \left( \theta  \right) \right)}{{{R}^{*}}^{3}}\left( 1+\frac{a\left( {{R}_{0}}\cos \theta -r \right)}{{{c}^{2}}}-\frac{{{v}^{2}}}{{{c}^{2}}} \right)-a\cos \left( \theta  \right)\frac{\left( t-{t}' \right)}{c{{R}^{*}}^{2}}$$

	field E = (R0 - (r_zap + (t - t_zap) * v_zap) * cos(theta))
		/ (R_lw_zap_3)
		* (1.0 + (a_zap * (R0*cos(theta) - r_zap) - v_zap * v_zap) / (g_c * g_c))
		- a_zap * cos(theta) * (t - t_zap) / (g_c * R_lw_zap_2);
#ifndef USE_NORM
	DBG_INFO("E = %0.25Lf\n", E);
	E *= multiplier_E;
#endif
	DBG_INFO("E = %0.25Lf\n", E);
	return E;
}

//#define OLD_DS_THETA_ALG



int integrand_phi_and_E(charge q, timevalue t, coordinate R0, coordinate r0, velocity v0, acceleration a0, angle theta, field * pE_minus_grad_phi_R0, field *pE_minus_1_c_dA_dt_R0, field *pE, coordinate r_min, potential *phi, potential *A)
{
	int err, error = 0;

	timevalue t_zap;
	coordinate r_zap;
	distance R_zap;
	distance R_lw_zap;

	velocity v_zap;
	acceleration a_zap;
	long double aR_zap;
	long double cos_alpha_zap;
	field E_minus_grad_varphi_R0;
	field E_minus_1_c_dA_dt_R0;
	field E;

	DBG_INFO("integrand_phi_and_E(q=%Lf t=%Lf, R0=%0.20Lf, r0=%0.20Lf, v0=%0.10Lf, a0=%0.10Lf)\n", q, t, R0, r0, v0, a0);

	*phi = 0.0;
	*A = 0.0;
	*pE_minus_grad_phi_R0 = 0.0;
	*pE_minus_1_c_dA_dt_R0 = 0.0;
	*pE = 0.0;


		err = calc_R_lw(q, t, R0, r0, v0, a0, theta, &t_zap, &r_zap, &R_zap, r_min, &R_lw_zap);
		if (0 != err)
		{
			error += 1;
		}
#ifdef ALGORITHM_VERSION_0
		v_zap = get_v(q, t_zap, v0, a0);
		a_zap = get_a(q, t_zap, a0);
#endif
#ifdef ALGORITHM_VERSION_1
		v_zap = get_v_ex1(t_zap, v0, q);
		a_zap = get_a_ex1(t_zap, q);
#endif
#ifdef ALGORITHM_VERSION_2
		v_zap = get_v_ex2(t_zap, v0, q);
		a_zap = get_a_ex2(t_zap, q);
#endif
		DBG_INFO("v_zap = %Lf ", v_zap);
		DBG_INFO("a_zap = %Lf ", a_zap);

		aR_zap = get_aR_zap(R0, theta, r_zap, a_zap);
		DBG_INFO("aR_zap = %Lf ", aR_zap);
		cos_alpha_zap = get_cos_alpha_zap(R0, theta, r_zap, R_zap);
		DBG_INFO("cos_alpha_zap = %Lf ", cos_alpha_zap);
		E_minus_grad_varphi_R0 = get_E_minus_grad_phi_R0(theta, v_zap, R_zap, aR_zap, R_lw_zap, cos_alpha_zap);
		E_minus_1_c_dA_dt_R0 = get_E_minus_1_c_dA_dt_R0(theta, v_zap, a_zap, R_zap, aR_zap, R_lw_zap);
		E = get_E(theta, t, t_zap, R0, r_zap, v_zap, a_zap);


		if (fabs(E_minus_grad_varphi_R0 + E_minus_1_c_dA_dt_R0 - E) > 1e-6)
		{
			printf("fabs(E_minus_grad_varphi_R0 + E_minus_1_c_dA_dt_R0 - E) = %e\n", fabs(E_minus_grad_varphi_R0 + E_minus_1_c_dA_dt_R0 - E));
			printf("E1 = %Le E2 = %Le E = %Le\n", E_minus_grad_varphi_R0, E_minus_1_c_dA_dt_R0, E);
			//int * p = 0;
			//p += 1;
		}


		DBG_INFO("theta = %Lf "
			"r_zap = %0.6Le "
			"R_zap %0.6Le "
			"R_lw_zap %0.6Le "
			"v_zap = %0.6Le "
			"t_zap = %0.6Le "
			"a_zap = %0.6Le "
			"aR_zap = %0.6Le "
			"E1 %Lf "
			"E2 %0.20Lf "
			"E %0.20Le "
			"\n"
			, theta
			, r_zap
			, R_zap
			, R_lw_zap
			, v_zap
			, t_zap
			, a_zap
			, aR_zap
			, E_minus_grad_varphi_R0
			, E_minus_1_c_dA_dt_R0
			, E
		);

		if (0.0 != R_lw_zap) {
			*phi = 1.0 / R_lw_zap;
			*A = 1.0 * cos(theta) * v_zap / (g_c * R_lw_zap);
			*pE_minus_grad_phi_R0 = E_minus_grad_varphi_R0;
			*pE_minus_1_c_dA_dt_R0 = E_minus_1_c_dA_dt_R0;
			*pE = E;
			DBG_INFO("phi = %Lf ", *phi);
			DBG_INFO("E1 = %Lf ", *pE_minus_grad_phi_R0);
			DBG_INFO("E2 = %Lf ", *pE_minus_1_c_dA_dt_R0);
			DBG_INFO("E2 = %Lf ", *pE);
		}
		else
		{
			DBG_INFO("0.0 != R_lw_zap ");
		}

	DBG_INFO("phi = %Lf E1 = %Lf E2 = %Lf\n", *phi, *pE_minus_grad_phi_R0, *pE_minus_1_c_dA_dt_R0);


	printf("t = %Lf theta = %Lf t_zap = %Le r_zap = %Le v_zap = %Le a_zap = %Le phi = %Le E1 = %Le E2 = %Le E = %Le\n", t, theta, t_zap, r_zap, v_zap, a_zap, *phi, *pE_minus_grad_phi_R0, *pE_minus_1_c_dA_dt_R0, *pE);
	return error;
}


/* скалярный потенциал Лиенара Вихерта зарядов равномерно распределённых по сферической поверхности радиуса r0 и движущихся из центра. */
/*
varphi := proc (q, t, r__0, v__0, a__0, R__0) options operator, arrow;
int(2*Pi*r(t, r__0, v__0, a__0)^2*sin(theta)*sigma(q, r__0)/K__zap(tzap(t, r__0, v__0, a__0, R__0, theta), r__0, v__0, a__0, R__0, theta), theta = 0 .. Pi) end proc;
*/
int integral_phi(charge q, timevalue t, coordinate R0, coordinate r0, velocity v0, acceleration a0, coordinate r_min, potential * result)
{
	int err, error = 0;
	int i;
	angle theta;
	timevalue t_zap;
	coordinate r_zap;
	distance R_zap;
	distance R_lw_zap;
	long double dS_dtheta;

	int N = 1000;
	angle dtheta = Pi / N;
	*result = 0.0;
	long double sigma0 = get_sigma(q, r0);
	long double ommited_S = 0.0;
	long double S0 = 4*Pi*r0*r0;
	long double S = 0.0;
	DBG_INFO("integral_phi(q=%Lf t=%Lf, R0=%Lf, r0=%Lf, v0=%Lf, a0=%Lf)\n", q, t, R0, r0, v0, a0);


	//printf("r = %Lf err = %d ", *r, err);

	for (i = 0; i <= N; ++i)
	{
		theta = (i * dtheta);

		err = calc_R_lw(q, t, R0, r0, v0, a0, theta, &t_zap, &r_zap, &R_zap, r_min, &R_lw_zap);
		if (0 != err)
		{
			error += 1;
		}
		if (i % 100 == 0)
			printf("%d %Lf R_lw_zap = %Lf t_zap = %Le\n", i, theta, R_lw_zap, t_zap);

#ifdef OLD_DS_THETA_ALG
		dS_dtheta = get_dS_dtheta(*r, theta);
		//printf("dS_dtheta = %Lf ", dS_dtheta);
#else
		dS_dtheta = get_dS_dtheta(r0, theta);
		//printf("r0 = %Lf dS_dtheta = %Lf \n", r0, dS_dtheta);
#endif

		S += dS_dtheta * dtheta;
		//printf("S = %Lf S0 = %Lf\n", S, S0);
		DBG_INFO("dS_dtheta = %Lf ", dS_dtheta);
		if (0.0 != R_lw_zap){
			*result += dS_dtheta / R_lw_zap * dtheta;
			DBG_INFO("result = %Lf ", *result);
		}
		else
		{
			ommited_S += dS_dtheta * dtheta;
			DBG_INFO("ommited_S = %Le ", ommited_S);
		}

		DBG_INFO("\n");
	}
	DBG_INFO("result = %Lf\n", *result);
	//printf("S = %Lf ommited_S = %Lf S0 = %Lf\n", S, ommited_S, S0);
	if (0.0 != ommited_S)
	{
		S -= ommited_S;
		DBG_INFO("corrected S = %Lf\n", S);
	}
	long double sigma = q / S;
	DBG_INFO("sigma0 = %Lf\n", sigma0);
	DBG_INFO("sigma = %Lf\n", sigma);
	*result *= sigma;
	DBG_INFO("result = %Lf\n", *result);
	return error;
}

int integral_phi_and_E(charge q, timevalue t, coordinate R0, coordinate r0, velocity v0, acceleration a0, field * pE_minus_grad_phi_R0, field *pE_minus_1_c_dA_dt_R0, field *pE, coordinate r_min, potential *phi, potential *A)
{
	int err, error = 0;
	int i;
	angle theta;
	timevalue t_zap;
	coordinate r_zap;
	distance R_zap;
	distance R_lw_zap;
	long double dS_dtheta;

	velocity v_zap;
	acceleration a_zap;
	long double aR_zap;
	long double cos_alpha_zap;
	field E_minus_grad_varphi_R0;
	field E_minus_1_c_dA_dt_R0;
	field E;

	int N = 10000;
	angle dtheta = Pi / N;
	*phi = 0.0;
	*A = 0.0;
	long double sigma0 = get_sigma(q, r0);
	long double ommited_S = 0.0;
	long double S0 = 4*Pi*r0*r0;
	long double S = 0.0;
	DBG_INFO("integral_phi_and_E(q=%Lf t=%Lf, R0=%0.20Lf, r0=%0.20Lf, v0=%0.10Lf, a0=%0.10Lf)\n", q, t, R0, r0, v0, a0);

	*pE_minus_grad_phi_R0  = 0.0;
	*pE_minus_1_c_dA_dt_R0 = 0.0;
	*pE                    = 0.0;

	for (i = 1; i <= N; ++i)
	{
		theta = (i * dtheta);

		err = calc_R_lw(q, t, R0, r0, v0, a0, theta, &t_zap, &r_zap, &R_zap, r_min, &R_lw_zap);
		if (0 != err)
		{
			error += 1;
		}
#ifdef ALGORITHM_VERSION_0
		v_zap = get_v(q, t_zap, v0, a0);
		a_zap = get_a(q, t_zap, a0);
#endif
#ifdef ALGORITHM_VERSION_1
		v_zap = get_v_ex1(t_zap, v0, q);
		a_zap = get_a_ex1(t_zap, q);
#endif
#ifdef ALGORITHM_VERSION_2
		v_zap = get_v_ex2(t_zap, v0, q);
		a_zap = get_a_ex2(t_zap, q);
#endif
		DBG_INFO("v_zap = %Lf ", v_zap);
		DBG_INFO("a_zap = %Lf ", a_zap);

		aR_zap = get_aR_zap(R0, theta, r_zap, a_zap);
		DBG_INFO("aR_zap = %Lf ", aR_zap);
		cos_alpha_zap = get_cos_alpha_zap(R0, theta, r_zap, R_zap);
		DBG_INFO("cos_alpha_zap = %Lf ", cos_alpha_zap);
		E_minus_grad_varphi_R0 = get_E_minus_grad_phi_R0 (theta, v_zap, R_zap, aR_zap, R_lw_zap, cos_alpha_zap);
		E_minus_1_c_dA_dt_R0 = get_E_minus_1_c_dA_dt_R0(theta, v_zap, a_zap, R_zap, aR_zap, R_lw_zap);
		E = get_E(theta, t, t_zap, R0, r_zap, v_zap, a_zap);
		E = get__E(theta, v_zap, a_zap, R_zap, aR_zap, R_lw_zap, cos_alpha_zap);

		// 2*Pi*r^2*sin(theta)
#ifdef OLD_DS_THETA_ALG
		dS_dtheta = get_dS_dtheta(*r, theta);
		//printf("dS_dtheta = %Lf ", dS_dtheta);
#else
		dS_dtheta = get_dS_dtheta(r0, theta);
		//printf("r0 = %Lf dS_dtheta = %Lf ", r0, dS_dtheta);
#endif
		S += dS_dtheta * dtheta;
		//printf("S = %Lf S0 = %Lf\n", S, S0);
		DBG_INFO("dS_dtheta = %Le ", dS_dtheta);

		if (i % 1000 == 0)
		if (0.0 != E_minus_1_c_dA_dt_R0)
			DBG_INFO("theta = %Lf "
				"r_zap = %0.6Le "
				"R_zap %0.6Le "
				"R_lw_zap %0.6Le "
				"v_zap = %0.6Le "
				"t_zap = %0.6Le "
				"a_zap = %0.6Le "
				"aR_zap = %0.6Le "
				"E1 %Lf "
				"E2 %0.20Lf "
				"e2 %0.20Le "
				"\n"
				, theta
				, r_zap
				, R_zap
				, R_lw_zap
				, v_zap
				, t_zap
				, a_zap
				, aR_zap
				, E_minus_grad_varphi_R0
				, E_minus_1_c_dA_dt_R0
				, E_minus_1_c_dA_dt_R0 * sigma0 * dS_dtheta * dtheta
			);

		if (0.0 != R_lw_zap){
			*phi                   += dS_dtheta * dtheta / R_lw_zap ;
			*A                     += dS_dtheta * dtheta * cos(theta) * v_zap / (g_c * R_lw_zap);
			*pE_minus_grad_phi_R0  += dS_dtheta * dtheta * E_minus_grad_varphi_R0;
			*pE_minus_1_c_dA_dt_R0 += dS_dtheta * dtheta * E_minus_1_c_dA_dt_R0;
			*pE                    += dS_dtheta * dtheta * E;
			DBG_INFO("phi = %Lf ", *phi);
			DBG_INFO("E1 = %Lf ", *pE_minus_grad_phi_R0);
			DBG_INFO("E2 = %Lf ", *pE_minus_1_c_dA_dt_R0);
			DBG_INFO("E2 = %Lf ", *pE);
		}
		else
		{
			ommited_S += dS_dtheta * dtheta;
			DBG_INFO("ommited_S = %Le ", ommited_S);
		}

		DBG_INFO("\n");
	}
	DBG_INFO("phi = %Lf E1 = %Lf E2 = %Lf\n", *phi, *pE_minus_grad_phi_R0, *pE_minus_1_c_dA_dt_R0);
	DBG_INFO("sigma0 = %Lf\n", sigma0);
	//printf("S = %Lf ommited_S = %Lf S0 = %Lf\n", S, ommited_S, S0);
	if (0.0 != ommited_S)
	{
		S -= ommited_S;
		DBG_INFO("corrected S = %Lf\n", S);
	}
	long double sigma = q / S;
	DBG_INFO("sigma = %Lf\n", sigma);
	*phi                   *= sigma;
	*A                     *= sigma;
	*pE_minus_grad_phi_R0  *= sigma;
	*pE_minus_1_c_dA_dt_R0 *= sigma;
	*pE                    *= sigma;

	printf("phi = %Lf E1 = %Lf E2 = %Lf E = %Lf\n", *phi, *pE_minus_grad_phi_R0, *pE_minus_1_c_dA_dt_R0, *pE);
	return error;
}


int dbl_integral_phi_and_E(charge q, timevalue t, coordinate R0, coordinate r0_min, velocity v0_min, acceleration a0_min, coordinate r0_max, velocity v0_max, acceleration a0_max, field * pE_minus_grad_phi_R0, field *pE_minus_1_c_dA_dt_R0, field *pE, coordinate r_min, potential *phi, potential *A)
{
	int error = 0;
	int j;

	int N_r0 = 1000;
	distance dr0 = (r0_max - r0_min) / N_r0;
	velocity dv0 = (v0_max - v0_min) / N_r0;
	acceleration da0 = (a0_max - a0_min) / N_r0;

	*phi = 0.0;
	*A = 0.0;

	*pE_minus_grad_phi_R0  = 0.0;
	*pE_minus_1_c_dA_dt_R0 = 0.0;
	*pE                    = 0.0;

	coordinate r0_pre = r0_min;
	long double V_of_q = 4.0 / 3.0 * Pi * ( (r0_max * r0_max * r0_max) - (r0_min * r0_min * r0_min) );

	for (j = 1; j <= N_r0; ++j)
	{
		coordinate r0 = r0_min + j * dr0;
		coordinate v0 = v0_min + j * dv0;
		coordinate a0 = a0_min + j * da0;

		charge dV_of_q = 4.0 / 3.0 * Pi * ( (r0 * r0 * r0) - (r0_pre * r0_pre * r0_pre) );
		charge dq = q * dV_of_q / V_of_q;

		field dE_minus_grad_phi_R0  = 0.0;
		field dE_minus_1_c_dA_dt_R0 = 0.0;
		field dE                    = 0.0;
		potential dphi              = 0.0;
		potential dA                = 0.0;


		integral_phi_and_E(dq, t, R0, r0, v0, a0, &dE_minus_grad_phi_R0, &dE_minus_1_c_dA_dt_R0, &dE, r_min, &dphi, &dA);

		printf("dq = %Le\n", dq);
		printf("q = %Lf\n", q);

		printf("dV_of_q = %Le V_of_q = %Le\n", dV_of_q, V_of_q);
		printf("r0 = %Lf v0 = %Lf a0 = %Lf\n", r0, v0, a0);

		printf("dphi = %Lf dA = %Le dE1 = %Lf dE2 = %Lf dE = %Lf\n", dphi, dA, dE_minus_grad_phi_R0, dE_minus_1_c_dA_dt_R0, dE);

		*phi                   += dphi;
		*A                     += dA;
		*pE_minus_grad_phi_R0  += dE_minus_grad_phi_R0;
		*pE_minus_1_c_dA_dt_R0 += dE_minus_1_c_dA_dt_R0;
		*pE                    += dE;
		printf("phi = %Lf A = %Le E1 = %Lf E2 = %Lf E = %Lf\n", *phi, *A, *pE_minus_grad_phi_R0, *pE_minus_1_c_dA_dt_R0, *pE);

		r0_pre = r0;
	}

	//printf("phi = %Lf E1 = %Lf E2 = %Lf\n", *phi, *pE_minus_grad_phi_R0, *pE_minus_1_c_dA_dt_R0);
	return error;
}

/* радиальная компонента векторного потенциала Лиенара Вихерта зарядов равномерно распределённых по сферической поверхности радиуса r0 и движущихся из центра. */

/*
A__R__0 := proc (q, t, r__0, v__0, a__0, R__0) options operator, arrow;
int(2*Pi*r(t, r__0, v__0, a__0)^2*sin(theta)*sigma(q, r__0)*v__r(t__zap, r__0, v__0, a__0)*cos(theta)
/
K__zap(tzap(t, r__0, v__0, a__0, R__0, theta), r__0, v__0, a__0, R__0, theta), theta = 0 .. Pi) end proc;
*/


/* Для расчёта радиальной компоненты электрического поля в точке наблюдения введём вспомогательную величину -
косинус угла между запаздывающим радиус-вектором (вектор из запаздывающего положения заряда в точку наблюдения)
и радиус-вектором из центра сферы в точку наблюдения
*/
/*
cos_alpha__zap := proc (t__zap, r__0, v__0, a__0, R__0, theta) options operator, arrow;
(R__0-r(t__zap, r__0, v__0, a__0)*cos(theta))/R__zap(t__zap, r__0, v__0, a__0, R__0, theta) end proc;
*/

/* Скалярное произведение ускорения частицы в запаздывающий момент времени на запаздывающий радиус-вектором (вектор из запаздывающего положения заряда в точку наблюдения) */
/*
aR__zap := proc (t__zap, r__0, v__0, a__0) options operator, arrow;
a__r(t__zap, r__0, v__0, a__0)*(R__0*cos(theta)-r(t__zap, r__0, v__0, a__0)) end proc;
*/

/* Первое слагаемое радиальной компоненты электрического поля - минус градиент скалярного потенциала */
/*
E_minus_grad_varphi__R__0 := proc (q, t__zap, r__0, v__0, a__0, R__0, theta) options operator, arrow;
sigma(q, r__0)*
(
R__zap(t__zap, r__0, v__0, a__0, R__0, theta)*cos_alpha__zap(t__zap, r__0, v__0, a__0, R__0, theta)*(1+aR__zap(t__zap, r__0, v__0, a__0)/c^2-v__r(t__zap, r__0, v__0, a__0)^2/c^2)
/
K__zap(t__zap, r__0, v__0, a__0, R__0, theta)
- v__r(t__zap, r__0, v__0, a__0)*cos(theta)/c
)
/
K__zap(t__zap, r__0, v__0, a__0, R__0, theta)^2
end proc;
*/

/* Интегрируя по поверхности сферы */
/*
E_minus_grad_&varphi;_integral__R__0 := proc (q, t, r__0, v__0, a__0, R__0) options operator, arrow;
int(2*Pi*r(t, r__0, v__0, a__0)^2*sin(theta)*E_minus_grad_varphi__R__0(q, tzap(t, r__0, v__0, a__0, R__0, theta), r__0, v__0, a__0, R__0, theta), theta = 0 .. Pi) end proc;
*/

/* Второе слагаемое компоненты электрического поля */
/*
E_minus_1_c_dA_dt__R__0 := proc (q, t__zap, r__0, v__0, a__0, R__0, theta) options operator, arrow;
cos(theta)*sigma(q, r__0)*
(
v__r(t__zap, r__0, v__0, a__0)*(R__zap(t__zap, r__0, v__0, a__0, R__0, theta)*(v__r(t__zap, r__0, v__0, a__0)^2/c-aR__zap(t__zap, r__0, v__0, a__0)/c-c)/K__zap(t__zap, r__0, v__0, a__0, R__0, theta)+c)/c^2
- a__r(t__zap, r__0, v__0, a__0)*R__zap(t__zap, r__0, v__0, a__0, R__0, theta)/c^2
)
/
K__zap(t__zap, r__0, v__0, a__0, R__0, theta)^2
end proc;
*/

/* Интегрируя по поверхности сферы */
/*
E_minus_1_c_dA_dt_integral__R__0 := proc (q, t, r__0, v__0, a__0, R__0) options operator, arrow;
int(2*Pi*r(t, r__0, v__0, a__0)^2*sin(theta)*E_minus_1_c_dA_dt__R__0(q, tzap(t, r__0, v__0, a__0, R__0, theta), r__0, v__0, a__0, R__0, theta), theta = 0 .. Pi) end proc;


*/
