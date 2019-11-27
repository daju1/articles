#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define Pi M_PI
#include "tzap.h"
#include <assert.h>
//#define USE_DEBUG
#include "dbg_info.h"
#include "testtzap.h"

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

	/* численный расчёта запаздывающего момента */
	err = calc_tzap(q, t, R0, r0, v0, a0, theta, r_min, pt_zap);
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
#ifdef ALGORITHM_VERSION_0
	v = get_v(*pt_zap, v0, a0);
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
	DBG_INFO("E_minus_grad_phi_R0 = %0.25e\n", E_minus_grad_phi_R0);
	return E_minus_grad_phi_R0;
}

/*double calc_E_minus_grad_varphi_R0(double t, double R0, double r0, double a0, double theta, double * pt_zap, double * pr_zap, double * pR_zap)
{

}*/


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

field get_E_minus_1_c_dA_dt_R0(angle theta, velocity v_zap, acceleration a_zap, distance R_zap, long double aR_zap, distance R_lw_zap)
{
	field E_minus_1_c_dA_dt_R0 =
		cos(theta) *
		(
			(v_zap / (g_c * g_c)) * ( (R_zap / R_lw_zap) * ( (v_zap * v_zap - aR_zap) / g_c - g_c)  + g_c)
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

//#define OLD_DS_THETA_ALG

/* скалярный потенциал Лиенара Вихерта зарядов равномерно распределённых по сферической поверхности радиуса r0 и движущихся из центра. */
/*
varphi := proc (q, t, r__0, v__0, a__0, R__0) options operator, arrow;
int(2*Pi*r(t, r__0, v__0, a__0)^2*sin(theta)*sigma(q, r__0)/K__zap(tzap(t, r__0, v__0, a__0, R__0, theta), r__0, v__0, a__0, R__0, theta), theta = 0 .. Pi) end proc;
*/
int integral_phi(charge q, timevalue t, coordinate R0, coordinate r0, velocity v0, double a0, coordinate r_min, potential * result)
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

int integral_phi_and_E(charge q, timevalue t, coordinate R0, coordinate r0, velocity v0, acceleration a0, field * pE_minus_grad_phi_R0, field *pE_minus_1_c_dA_dt_R0, coordinate r_min, potential *phi)
{
	int err, error = 0;
	int i;
	angle theta;
	timevalue t_zap;
	coordinate r_zap;
	distance R_zap;
	distance R_lw_zap;
	long double dS_dtheta;

	double v_zap;
	acceleration a_zap;
	long double aR_zap;
	long double cos_alpha_zap;
	double E_minus_grad_varphi_R0;
	double E_minus_1_c_dA_dt_R0;

	int N = 10000;
	angle dtheta = Pi / N;
	*phi = 0.0;
	long double sigma0 = get_sigma(q, r0);
	long double ommited_S = 0.0;
	long double S0 = 4*Pi*r0*r0;
	long double S = 0.0;
	DBG_INFO("integral_phi_and_E(q=%Lf t=%Lf, R0=%0.20Lf, r0=%0.20Lf, v0=%0.10Lf, a0=%0.10Lf)\n", q, t, R0, r0, v0, a0);

	*pE_minus_grad_phi_R0 = 0.0;
	*pE_minus_1_c_dA_dt_R0 = 0.0;

	for (i = 1; i <= N; ++i)
	{
		theta = (i * dtheta);

		err = calc_R_lw(q, t, R0, r0, v0, a0, theta, &t_zap, &r_zap, &R_zap, r_min, &R_lw_zap);
		if (0 != err)
		{
			error += 1;
		}
#ifdef ALGORITHM_VERSION_0
		v_zap = get_v(t_zap, v0, a0);
		a_zap = get_a(t_zap, a0);
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
		E_minus_1_c_dA_dt_R0   = get_E_minus_1_c_dA_dt_R0(theta, v_zap, a_zap, R_zap, aR_zap, R_lw_zap);

#ifdef OLD_DS_THETA_ALG
		dS_dtheta = get_dS_dtheta(*r, theta);
		//printf("dS_dtheta = %Lf ", dS_dtheta);
#else
		dS_dtheta = get_dS_dtheta(r0, theta);
		//printf("r0 = %Lf dS_dtheta = %Lf ", r0, dS_dtheta);
#endif
		S += dS_dtheta * dtheta;
		//printf("S = %Lf S0 = %Lf\n", S, S0);
		DBG_INFO("dS_dtheta = %Lf ", dS_dtheta);

		//if (i % 100 == 0)
		if (0.0 != E_minus_1_c_dA_dt_R0)
			DBG_INFO("theta = %Lf "
				"r_zap = %0.6e "
				"R_zap %0.6e "
				"R_lw_zap %0.6e "
				"v_zap = %0.6e "
				"t_zap = %0.6e "
				"a_zap = %0.6e "
				"aR_zap = %0.6e "
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
			*pE_minus_grad_phi_R0  += dS_dtheta * dtheta * E_minus_grad_varphi_R0;
			*pE_minus_1_c_dA_dt_R0 += dS_dtheta * dtheta * E_minus_1_c_dA_dt_R0;
			DBG_INFO("phi = %Lf ", *phi);
			DBG_INFO("E1 = %Lf ", *pE_minus_grad_phi_R0);
			DBG_INFO("E2 = %Lf ", *pE_minus_1_c_dA_dt_R0);
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
	*pE_minus_grad_phi_R0  *= sigma;
	*pE_minus_1_c_dA_dt_R0 *= sigma;

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

int test_v0()
{
	int error = 0;
	coordinate r_zap;
	distance R_zap, R_lw_zap;
	potential phi_lw;
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

	error = calc_tzap(q, /*t*/0.0, /*R0*/0.0, /*r0*/2.0, /*v0*/1.0, /*a0*/0.0, /*theta*/0.0, r_min, &t_zap);
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
	error = integral_phi_and_E(q, t, R0, r0, v0, a0, &E_minus_grad_phi_R0, &E_minus_1_c_dA_dt_R0, r_min, &phi_lw);
	printf("phi_lw = %Lf E1=%Lf E2 = %Lf\n", phi_lw, E_minus_grad_phi_R0, E_minus_1_c_dA_dt_R0);

	/*infinity error result*/
	error = integral_phi(q, t, -1.0000000000000142, 1, 0, 0, r_min, &phi_lw);
	printf("phi_lw = %Lf\n", phi_lw);
	phi_lw = integral_phi_and_E(q, t, -1.0000000000000142, 1, 0, 0, &E_minus_grad_phi_R0, &E_minus_1_c_dA_dt_R0, r_min, &phi_lw);
	printf("phi_lw = %Lf E1=%Lf E2 = %Lf\n", phi_lw, E_minus_grad_phi_R0, E_minus_1_c_dA_dt_R0);

	/*hung
	calc_tzap(t=5.000000, R0=-1.500000, r0=2.000000, v0=0.0, a0=1.000000, theta=0.000000
	*/
	error = calc_tzap(q, /*t*/5.0, /*R0*/-1.5, /*r0*/2.0, /*v0*/0.0, /*a0*/1.0, /*theta*/0.0, r_min, &t_zap);
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
	error = integral_phi_and_E(-q, 5.0, -1.5, 2.0, 0.0, 1.0, &E_minus_grad_phi_R0, &E_minus_1_c_dA_dt_R0, r_min, &phi_lw);
	printf("phi_lw = %Lf E1=%Lf E2 = %Lf\n", phi_lw, E_minus_grad_phi_R0, E_minus_1_c_dA_dt_R0);

	{
		field E1_p, E1_n, E2_p, E2_n;
		potential phi_p, phi_n;
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

				error = integral_phi_and_E(+q, ti, R0_i, r0, v0_p, a0_p, &E1_p, &E2_p, r_min, &phi_p);
				error = integral_phi_and_E(-q, ti, R0_i, r0, v0_n, a0_n, &E1_n, &E2_n, r_min, &phi_n);
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
		error = integral_phi_and_E(q, t, R0, r0, 0.0, 0.0, &E_minus_grad_phi_R0, &E_minus_1_c_dA_dt_R0, r_min, &phi_lw);
		printf("phi_lw = %Lf E1=%Lf E2 = %Lf\n", phi_lw, E_minus_grad_phi_R0, E_minus_1_c_dA_dt_R0);

		error = integral_phi_and_E(-q, t, R0, r0, 0.0, 0.0, &E_minus_grad_phi_R0, &E_minus_1_c_dA_dt_R0, r_min, &phi_lw);
		printf("phi_lw = %Lf E1=%Lf E2 = %Lf\n", phi_lw, E_minus_grad_phi_R0, E_minus_1_c_dA_dt_R0);

	}

	return 0;
}

extern int v_n_t; // итератор полноты заполения двумерных массивов по координате времени

int test_v1()
{
	/* начальный радиус заряженной сферы в момент t=t_start */
	coordinate r0_pos = 0.1;
	coordinate r0_neg = 0.1;

	/* минимально возможный радиус заряженной сферы (из соображений упругости) в момент более ранний чем t=t_start */
	coordinate r_min_pos = 0.01;
	coordinate r_min_neg = 0.01;

	/* скорость в момент t_start*/
	velocity v0_pos = 0.0;
	velocity v0_neg = 0.0;

	/* ускорение вызванное причинами неэлектрического характера, например вледствие подвода энергии извне */
	acceleration a0_pos = 0.1*g_c;
	acceleration a0_neg = 1.0*g_c;

	mass m_pos = 10.0;
	mass m_neg = 1.0;

	/* Заряд сферы */
	charge q = 1.0;

	/* время действия ускорения вызванного причинами неэлектрического характера, например вледствие подвода энергии извне */
	double t_a0 = 1.0;

	return do_v1_calc(q, m_pos, m_neg, r0_pos, r0_neg, v0_pos, v0_neg, a0_pos, a0_neg, t_a0);
}

/*
sage: attach("copper_explosion_lw.sage")
M__Cu_sgs =  0.00140115032350105
should be 0.001401150324
a(2567) =  1.06780733067000
should be 1.067807331
c__p =  [(300, 385.000000000000), (400, 397.700000000000), (500, 408.000000000000), (600, 416.900000000000), (700, 425.100000000000), (800, 432.900000000000), (900, 441.700000000000), (1000, 451.400000000000), (1100, 464.300000000000), (1200, 480.800000000000), (1300, 506.600000000000), (1357.60000000000, 525.200000000000)]
Q__nagrev_do_T_pl =  466024.203987
should be 466024.2040
c__p(1357)=  525.001217859
should be 525.001217859185
c__p_exptrapotation=  33.3617273900798
should be 33.3617273900798
0.652969964175658
should be .6529699644
V__Cu_sgs__crit =  0.000718538627436435
should be 0.7185386277e-3
R__i =  0.000555631898340631
should be 0.5556318983e-3
Q__pl =  0.286911910827956
should be .2869119109
Q__isp =  6.71624317090642
should be 6.716243174
Q__nagrev =  4.99696901906082
should be 4.99696902105744
Q__ionization =  16.4267930476864
should be 16.42679305
V__Cu_1atm =  0.0000130035233983192
should be 0.1300352340e-4
V__Cu_450atm =  2.88967186629317e-8
should be 2.889671867*10^(-8)
V__Cu__crit =  7.18538627436435e-10
should be 7.185386277*10^(-10)
A =  12.7480190448392
should be 12.74801905
Delta_E =  93.8250638066792
should be 93.8250637789426
Delta_t =  0.000347500236321034
should be 0.347500236218306e-3
Delta_T =  127548.050829933
should be 1.27548050741264*10^5
T__e =  134641.050829933
should be 1.34641050741264*10^5
11.6024640799239
should be 11.6024640725366
7.18538627436435e-10
should be 7.185386277*10^(-10)
7.18538627436435e-10
should be 7.185386273*10^(-10)
n__i =  1.84797988858465e28
463689.588837354
should be 4.636895888*10^5
2.02023137307735e6
should be 2.02023137241212*10^6
1357.55450752652
should be 1357.554508
5914.67713054400
should be 5914.67712890453
(2.36433873744671e-12)*r__0*sqrt(T/m)
(2.41093551366813e-7)*r__0/(R_init*m*sqrt(T/m))
0.000433908766013665*r0/(m*sqrt(T/m))
should be 0.433908765869441e-3*r0/(sqrt(T/m)*m)
Compiling ./tzap.spyx...
q =  2.12744210611535 кулон
a0_pos =  3.81193295102697e-7
a0_neg =  0.000130201300422018
nr = 10000 nt = 10000
dr = 0.00000027781594917032 dt = 0.00000017375011816052
t_a0 = 0.00034750023632103418
m_pos = 0.00000141116074503615 m_neg = 0.00000000001209584880
r0_pos = 0.00027781594917031547 r0_neg = 0.00027781594917031547
v0_pos = 0.00000000000000000000 v0_neg = 0.00000000000000000000
a0_pos = 0.00000038119329510270 a0_neg = 0.00013020130042201811
sage:
*/

int copper_explosion_lw_v1()
{
	charge q = 2.12744210611535; // кулон
	acceleration a0_pos = 3.81193295102697e-7;
	acceleration a0_neg = 0.000130201300422018;
	int nr = 10000, nt = 10000;
	timespan t_a0 = 0.00034750023632103418;
	mass m_pos = 0.00000141116074503615, m_neg = 0.00000000001209584880;
	coordinate r0_pos = 0.00027781594917031547, r0_neg = 0.00027781594917031547;
	velocity v0_pos = 0.00000000000000000000, v0_neg = 0.00000000000000000000;
	return do_v1_calc(q, m_pos, m_neg, r0_pos, r0_neg, v0_pos, v0_neg, a0_pos, a0_neg, t_a0);
}

int do_v1_calc(charge q, mass m_pos, mass m_neg, coordinate r0_pos, coordinate r0_neg, velocity v0_pos, velocity v0_neg, acceleration a0_pos, acceleration a0_neg, timespan t_a0)
{
	printf("do_v1_calc\n");
	printf("q = %0.20Le\n", q);
	printf("t_a0 = %0.20Le\n", t_a0);
	printf("m_pos = %0.20Le m_neg = %0.20Le\n", m_pos, m_neg);
	printf("r0_pos = %0.20Le r0_neg = %0.20Le\n", r0_pos, r0_neg);
	printf("v0_pos = %0.20Le v0_neg = %0.20Le\n", v0_pos, v0_neg);
	printf("a0_pos = %0.20Le a0_neg = %0.20Le\n", a0_pos, a0_neg);
#ifdef USE_NORM
	// q = q_calc * k_q
	// m = m_calc * k_m
	// r = r_calc * k_r
	// t = t_calc * k_t
	// v = v_valc * k_v
	// c = c_calc * k_v -> LIGHT_VELONCITY = g_c * k_v
	// a = a_calc * k_a
	// E = E_calc * k_E

    // SI :
	// k_E = (LIGHT_VELONCITY * LIGHT_VELONCITY / 10000000.0) * k_q / (k_r * k_r)
	// k_E = SI_multiplier_E * k_q / (k_r * k_r)
	// sgs :
	// k_E = k_q / (k_r * k_r)

	// k_a = k_E * k_q  / k_m
	// k_t = k_v / k_a
	// k_v = k_t * k_a
	// k_a * k_r = k_v * k_v


	long double k_q = 1.0;
	long double k_m = 0.0001;
	long double k_r = 0.001;
	long double k_E = multiplier_E * k_q / (k_r * k_r);
	long double k_a = k_E * k_q  / k_m;
	long double k_v = sqrt(k_a * k_r);
	long double k_t = k_v  / k_a;
	g_c = LIGHT_VELONCITY / k_v;

	printf("k_q = %0.20Lf\n", k_q);
	printf("k_m = %0.20Lf\n", k_m);
	printf("k_r = %0.20Lf\n", k_r);
	printf("m_E = %0.20Lf\n", multiplier_E);
	printf("k_E = %0.20Lf\n", k_E);
	printf("k_a = %0.20Lf\n", k_a);
	printf("k_v = %0.20Lf\n", k_v);
	printf("g_c = %0.20Lf\n", g_c);
	printf("defc= %0.20f\n", LIGHT_VELONCITY);
	printf("k_t = %0.20Lf\n", k_t);

	charge q_calc = q / k_q;
	acceleration a0_pos_calc = a0_pos / k_a;
	acceleration a0_neg_calc = a0_neg / k_a;
	timespan t_a0_calc = t_a0 / k_t;
	mass m_pos_calc = m_pos / k_m, m_neg_calc = m_neg / k_m;
	coordinate r0_pos_calc = r0_pos / k_r, r0_neg_calc = r0_neg / k_r;
	velocity v0_pos_calc = v0_pos / k_v, v0_neg_calc = v0_neg / k_v;
#else
	charge q_calc = q;
	acceleration a0_pos_calc = a0_pos;
	acceleration a0_neg_calc = a0_neg;
	timespan t_a0_calc = t_a0;
	mass m_pos_calc = m_pos, m_neg_calc = m_neg;
	coordinate r0_pos_calc = r0_pos, r0_neg_calc = r0_neg;
	velocity v0_pos_calc = v0_pos, v0_neg_calc = v0_neg;
#endif

	coordinate r_finish_calc = r0_pos_calc * 10;
	distance dr_calc = r_finish_calc / 1000;
	coordinate r_min_pos_calc = dr_calc;
	coordinate r_min_neg_calc = dr_calc;

	set_dr(dr_calc);
	set_r_finish(r_finish_calc);
	return do_v1_calc_priv(q_calc, m_pos_calc, m_neg_calc, r0_pos_calc, r0_neg_calc, v0_pos_calc, v0_neg_calc, a0_pos_calc, a0_neg_calc, t_a0_calc , r_min_pos_calc, r_min_neg_calc);
}

int calc_E(charge q, timevalue t, coordinate R0,
			coordinate r0_pos, coordinate r0_neg,
			velocity v0_pos, velocity v0_neg,
			acceleration a0_pos, acceleration a0_neg,
			coordinate r_min_pos, coordinate r_min_neg,
			field * E_minus_grad_phi_R0_pos, field * E_minus_1_c_dA_dt_R0_pos,
			field * E_minus_grad_phi_R0_neg, field * E_minus_1_c_dA_dt_R0_neg,
			field * E1, field * E2, field * E)
{
	int error = 0, err;
	potential phi_lw_pos;
	potential phi_lw_neg;

	printf("\n");
	error = integral_phi_and_E(+q, t, R0, r0_pos, v0_pos, a0_pos, E_minus_grad_phi_R0_pos, E_minus_1_c_dA_dt_R0_pos, r_min_pos, &phi_lw_pos);
	printf("pos err = %d phi_lw = %Lf E1=%Lf E2 = %0.20Lf\n", error, phi_lw_pos, *E_minus_grad_phi_R0_pos, *E_minus_1_c_dA_dt_R0_pos);

	printf("\n");
	error = integral_phi_and_E(-q, t, R0, r0_neg, v0_neg, a0_neg, E_minus_grad_phi_R0_neg, E_minus_1_c_dA_dt_R0_neg, r_min_neg, &phi_lw_neg);
	printf("neg err = %d phi_lw = %Lf E1=%Lf E2 = %0.20Lf\n", error, phi_lw_neg, *E_minus_grad_phi_R0_neg, *E_minus_1_c_dA_dt_R0_neg);

	*E1 = *E_minus_grad_phi_R0_pos + *E_minus_grad_phi_R0_neg;
	*E2 = *E_minus_1_c_dA_dt_R0_pos + *E_minus_1_c_dA_dt_R0_neg;

	*E = *E1 + *E2;
}


int do_v1_calc_priv(charge q, mass m_pos, mass m_neg, coordinate r0_pos, coordinate r0_neg, velocity v0_pos, velocity v0_neg, acceleration a0_pos, acceleration a0_neg, timespan t_a0, coordinate r_min_pos, coordinate r_min_neg)
{
	printf("sizeof(double) %ld\n", sizeof(double));
	printf("sizeof(long double) %ld\n", sizeof(long double));

	printf ("nr = %d ", get_nr());
	printf ("nt = %d\n", get_nt());
	printf ("dr = %0.20Lf ", get_dr());

	printf("do_v1_calc_priv\n");
	printf("q = %0.20Le\n", q);
	printf("t_a0 = %0.20Le\n", t_a0);
	printf("m_pos = %0.20Le m_neg = %0.20Le\n", m_pos, m_neg);
	printf("r0_pos = %0.20Le r0_neg = %0.20Le\n", r0_pos, r0_neg);
	printf("v0_pos = %0.20Le v0_neg = %0.20Le\n", v0_pos, v0_neg);
	printf("a0_pos = %0.20Le a0_neg = %0.20Le\n", a0_pos, a0_neg);

#ifdef ALGORITHM_VERSION_1
	init_array_1(a0_pos, v0_pos, r0_pos, a0_neg, v0_neg, r0_neg);
#endif
	v_n_t = 0;
	timespan dt = 1e-7;
	printf("dt = %Le\n", dt);
	while (v_n_t < get_nt())
	{
		int error = 0, err;
		potential phi_lw_pos;
		potential phi_lw_neg;
		field E_minus_grad_phi_R0_pos, E_minus_1_c_dA_dt_R0_pos;
		field E_minus_grad_phi_R0_neg, E_minus_1_c_dA_dt_R0_neg;
		coordinate r_pos;
		coordinate r_neg;
		field E1, E2, E;

		timevalue t = v_t[v_n_t];
		//printf("t = %Le v_n_t = %d\n", t, v_n_t);
		#if 0
		for (int v_n_r = 0; v_n_r < get_nr(); ++v_n_r)
		{
			coordinate R0 = v_n_r * get_dr();
			//printf("R0 = %Lf v_n_r = %d dr = %Lf\n", R0, v_n_r, dr);

			calc_E(q, t, R0,
				r0_pos, r0_neg,
				v0_pos, v0_neg,
				a0_pos, a0_neg,
				r_min_pos, r_min_neg,
				&E_minus_grad_phi_R0_pos, &E_minus_1_c_dA_dt_R0_pos,
				&E_minus_grad_phi_R0_neg, &E_minus_1_c_dA_dt_R0_neg,
				&E1, &E2, &E);
			if (fabs(E) > 1e-20){
				#if 0
				printf(
					"R0 = %Lf t = %Lf "
					"E = % 0.20f 		field E_pos, E_neg; // электрическое поле в облаасти нахождения положительной и отрицательной обкладки
E1 = % 0.20f E2 = % 0.20f "
					"phi_lw_pos = %Lf phi_lw_neg = %Lf\n"
					, R0, t
					, E
					, E1, E2
					, phi_lw_pos, phi_lw_neg);
				#endif
#if 0
				printf(
					"R0 = %0.10Lf t = %0.10Lf "
					"E1_pos % 0.20f "
					"E1_neg % 0.20f "
					"phi_lw_pos = %Lf phi_lw_neg = %Lf "
					//"E2_pos % 0.20f "
					//"E2_neg % 0.20f "
					"\n"
					, R0, t
					, E_minus_grad_phi_R0_pos
					, E_minus_grad_phi_R0_neg
					//, E_minus_1_c_dA_dt_R0_pos
					//, E_minus_1_c_dA_dt_R0_neg
					, phi_lw_pos, phi_lw_neg
					);
#endif
			}
#ifdef ALGORITHM_VERSION_1
			set_E_ex_1(v_n_t, v_n_r, E);
#endif
		}
		printf("\n");
		#endif

#ifdef ALGORITHM_VERSION_1
		acceleration a_pos = get_a_ex1(t, +q);
		acceleration a_neg = get_a_ex1(t, -q);

		velocity v_pos = get_v_ex1(t, v0_pos, +q);
		velocity v_neg = get_v_ex1(t, v0_neg, -q);

		distance s_pos = get_s_ex1(t, v0_pos, +q);
		distance s_neg = get_s_ex1(t, v0_neg, -q);
		distance ds = s_neg - s_pos;

		err = get_r_ex1(+q, t, r0_pos, v0_pos, r_min_pos, &r_pos, 1);
		err = get_r_ex1(-q, t, r0_neg, v0_neg, r_min_neg, &r_neg, 1);

		printf("t = %Lf\n"
			"a_pos = % 0.20Le, a_neg = % 0.20Le\n"
			"v_pos = % 0.20Le, v_neg = % 0.20Le\n"
			"vcpos = % 0.20Le, vcneg = % 0.20Le\n"
			"s_pos = % 0.20Le, s_neg = % 0.20Le s_neg - s_pos = % 0.20Le  \n"
			"r_pos = % 0.20Le, r_neg = % 0.20Le r_neg - r_pos = % 0.20Le \n"
			, t
			, a_pos, a_neg
			, v_pos, v_neg
			, v_pos / g_c, v_neg / g_c
			, s_pos, s_neg, (s_neg - s_pos)
			, r_pos, r_neg, (r_neg - r_pos)
			);


		//
		field E_pos, E_neg; // электрическое поле в облаасти нахождения положительной и отрицательной обкладки
		field E1_pos, E1_neg; // электрическое поле в облаасти нахождения положительной и отрицательной обкладки
		field E2_pos, E2_neg; // электрическое поле в облаасти нахождения положительной и отрицательной обкладки

		printf("\nCalc E on positive side");

		calc_E(q, t,
			r_pos, // координата наблюдения совпадает с координатой обкладки
			r0_pos, r0_neg,
			v0_pos, v0_neg,
			a0_pos, a0_neg,
			r_min_pos, r_min_neg,
			&E_minus_grad_phi_R0_pos, &E_minus_1_c_dA_dt_R0_pos,
			&E_minus_grad_phi_R0_neg, &E_minus_1_c_dA_dt_R0_neg,
			&E1_pos, &E2_pos, &E_pos);

		printf(
			"R0_pos = %0.10Lf t = %0.10Lf "
			"E1_pos %Le "
			"E1_neg %Le "
			"E1 %Le "
			"E2_pos %Le "
			"E2_neg %Le "
			"E2 %0.20Le "
			"\n"
			, r_pos, t
			, E_minus_grad_phi_R0_pos
			, E_minus_grad_phi_R0_neg
			, E_minus_grad_phi_R0_pos + E_minus_grad_phi_R0_neg

			, E_minus_1_c_dA_dt_R0_pos
			, E_minus_1_c_dA_dt_R0_neg
			, E_minus_1_c_dA_dt_R0_pos + E_minus_1_c_dA_dt_R0_neg
			);

		printf("\nCalc E on negative side");

		calc_E(q, t,
			r_neg, // координата наблюдения совпадает с координатой обкладки
			r0_pos, r0_neg,
			v0_pos, v0_neg,
			a0_pos, a0_neg,
			r_min_pos, r_min_neg,
			&E_minus_grad_phi_R0_pos, &E_minus_1_c_dA_dt_R0_pos, 
			&E_minus_grad_phi_R0_neg, &E_minus_1_c_dA_dt_R0_neg,
			&E1_neg, &E2_neg, &E_neg);

		printf(
			"R0_neg = %0.10Lf t = %0.10Lf "
			"E1_pos %Le "
			"E1_neg %Le "
			"E1 %Le "
			"E2_pos %Le "
			"E2_neg %Le "
			"E2 %0.20Le "
			"\n"
			, r_pos, t
			, E_minus_grad_phi_R0_pos
			, E_minus_grad_phi_R0_neg
			, E_minus_grad_phi_R0_pos + E_minus_grad_phi_R0_neg

			, E_minus_1_c_dA_dt_R0_pos
			, E_minus_1_c_dA_dt_R0_neg
			, E_minus_1_c_dA_dt_R0_pos + E_minus_1_c_dA_dt_R0_neg
			);

		//printf(
		//	"E_pos = % 0.20Le, E_neg = % 0.20Le (E_pos - E_neg) = % 0.20Le\n"
		//	, E_pos, E_neg, (E_pos - E_neg)
		//	);
		//

		if (v_n_t < get_nt() - 1)
		{
			v_t[v_n_t + 1] = v_t[v_n_t] + dt;
			timevalue t1 = v_t[v_n_t + 1];

			printf("t = %Le t1 = %Le v_n_t = %d dt = %Le\n", t, t1, v_n_t, dt);
			if (0 != set_a_ex1(t1, r_neg, a0_neg, t_a0, -q, m_neg, E_neg, E1_neg, E2_neg))
			{
				error += 1;
				int * p = 0;
				*p += 1;
			}
			if (0 != set_a_ex1(t1, r_pos, a0_pos, t_a0, +q, m_pos, E_pos, E1_pos, E2_pos))
			{
				error += 1;
			}

			velocity v_pos, v_neg;
			if (0 != set_v_ex1(t1, v0_neg, a0_neg, t_a0, -q, m_neg, &v_neg))
			{
				error += 1;
			}
			if (0 != set_v_ex1(t1, v0_pos, a0_pos, t_a0, +q, m_pos, &v_pos))
			{
				error += 1;
			}

			printf(
				"v_pos = % 0.20Le, v_neg = % 0.20Le\n"
				, v_pos, v_neg
				);
			

			if (0 == error) {
				set_s_ex1(t1, r0_pos, v0_pos, r_min_pos, +q);
				set_s_ex1(t1, r0_neg, v0_neg, r_min_neg, -q);
			}

		}
#endif
		if (0 == error)
		{
			++v_n_t;
		}
		else
		{
			v_n_t-=2;
			if (v_n_t < 0)
				v_n_t = 0;
			dt /= 2.0;
			printf("dt = %Le\n", dt);
		}
		
	}
}



