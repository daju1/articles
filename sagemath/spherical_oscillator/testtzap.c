#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define Pi M_PI
#include "tzap.h"
#include <assert.h>
//#define USE_DEBUG
#include "dbg_info.h"

double get_sigma(double q, double r)
{
	double sigma = q / (4*Pi*r*r);
	return sigma;
}

double get_dS_dtheta(double r, double theta)
{
	double dS_dtheta = 2*Pi*r*r*sin(theta);
	return dS_dtheta;
}

/* Радиус Лиенара Вихерта */
int calc_R_lw(double q, double t, double R0, double r0, double v0, double a0, double theta, double * pt_zap, double * pr_zap, double * pR_zap, double r_min, double * R_lw_zap)
{
//#define DBG_INFO printf
	int err, error = 0;
	double v;

	/* численный расчёта запаздывающего момента */
	err = calc_tzap(q, t, R0, r0, v0, a0, theta, r_min, pt_zap);
	if (0 != err)
	{
		error += 1;
	}
	DBG_INFO("theta = %f t_zap = %f ", theta, *pt_zap);
	/* Запаздывающий радиус в зависимости от текущего момента */
#ifdef ALGORITHM_VERSION_0
	err = get_r(q, *pt_zap, r0, v0, a0, r_min, pr_zap); /* расстояние от заряда до центра сферы в запаздывающий момент времени */
#endif
#ifdef ALGORITHM_VERSION_1
	err = get_r_ex1(q, *pt_zap, r0, v0, r_min, pr_zap);
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
	DBG_INFO("r_zap = %f ", *pr_zap);
	DBG_INFO("R_zap = %f ", *pR_zap);
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
	*R_lw_zap = get_R(R0, *pr_zap, theta) - (v / c) * (R0 * cos(theta) - *pr_zap);
	DBG_INFO("R_lw_zap = %f ", *R_lw_zap);

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

double get_cos_alpha_zap(double R0, double theta, double r_zap, double R_zap)
{
	double cos_alpha_zap = (R0 - r_zap*cos(theta)) / R_zap;
	return cos_alpha_zap;
}

/* Скалярное произведение ускорения частицы в запаздывающий момент времени на запаздывающий радиус-вектор (вектор из запаздывающего положения заряда в точку наблюдения) */
/*
aR__zap := proc (t__zap, r__0, v__0, a__0) options operator, arrow;
a__r(t__zap, r__0, v__0, a__0)*(R__0*cos(theta)-r(t__zap, r__0, v__0, a__0)) end proc;
*/

double get_aR_zap(double R0, double theta, double r_zap, double a_zap)
{
	double aR_zap = a_zap*(R0*cos(theta) - r_zap);
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

double get_E_minus_grad_phi_R0(double theta, double v_zap, double R_zap, double aR_zap, double R_lw_zap, double cos_alpha_zap)
{
	double E_minus_grad_phi_R0 =
		(
			(cos_alpha_zap * R_zap / R_lw_zap) * (1.0 + (aR_zap - v_zap * v_zap) / (c * 2) )
			- v_zap*cos(theta) / c
		)
		/ (R_lw_zap * R_lw_zap);
#ifdef SI
	E_minus_grad_phi_R0 *= (c*c)/(10000000.0);
#endif
	DBG_INFO("E_minus_grad_phi_R0 = %f ", E_minus_grad_phi_R0);
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

double get_E_minus_1_c_dA_dt_R0(double theta, double v_zap, double a_zap, double R_zap, double aR_zap, double R_lw_zap)
{
	double E_minus_1_c_dA_dt_R0 =
		cos(theta) *
		(
			(v_zap / (c * c)) * ( (R_zap / R_lw_zap) * ( (v_zap * v_zap - aR_zap) / c - c)  + c)
			- a_zap * R_zap / (c * c)
		)
		/
		(R_lw_zap * R_lw_zap);
	DBG_INFO("E_minus_1_c_dA_dt_R0 = %f ", E_minus_1_c_dA_dt_R0);
#ifdef SI
	E_minus_1_c_dA_dt_R0 *= (c*c)/(10000000.0);
#endif
	return E_minus_1_c_dA_dt_R0;
}

//#define OLD_DS_THETA_ALG

/* скалярный потенциал Лиенара Вихерта зарядов равномерно распределённых по сферической поверхности радиуса r0 и движущихся из центра. */
/*
varphi := proc (q, t, r__0, v__0, a__0, R__0) options operator, arrow;
int(2*Pi*r(t, r__0, v__0, a__0)^2*sin(theta)*sigma(q, r__0)/K__zap(tzap(t, r__0, v__0, a__0, R__0, theta), r__0, v__0, a__0, R__0, theta), theta = 0 .. Pi) end proc;
*/
int integral_phi(double q, double t, double R0, double r0, double v0, double a0, double r_min, double * result, double * r)
{
	int err, error = 0;
	int i;
	double theta,
		t_zap, r_zap, R_zap,
		R_lw_zap, dS_dtheta;

	int N = 1000;
	double dtheta = Pi / N;
	*result = 0.0;
	double sigma0 = get_sigma(q, r0);
	double ommited_S = 0.0;
	double S0 = 4*Pi*r0*r0;
	double S = 0.0;
	DBG_INFO("integral_phi(q=%f t=%f, R0=%f, r0=%f, v0=%f, a0=%f)\n", q, t, R0, r0, v0, a0);


#ifdef ALGORITHM_VERSION_0
	err = get_r(q, t, r0, v0, a0, r_min, r);
#endif
#ifdef ALGORITHM_VERSION_1
	err = get_r_ex1(q, t, r0, v0, r_min, r);
#endif
#ifdef ALGORITHM_VERSION_2
	err = get_r_ex2(q, t, r0, v0, r_min, r);
#endif
	if (0 != err)
	{
		error += 1;
	}

	//printf("r = %f err = %d ", *r, err);

	for (i = 0; i <= N; ++i)
	{
		theta = (i * dtheta);

		err = calc_R_lw(q, t, R0, r0, v0, a0, theta, &t_zap, &r_zap, &R_zap, r_min, &R_lw_zap);
		if (0 != err)
		{
			error += 1;
		}
		if (i % 100 == 0)
			printf("%d %f R_lw_zap = %f t_zap = %e\n", i, theta, R_lw_zap, t_zap);

#ifdef OLD_DS_THETA_ALG
		dS_dtheta = get_dS_dtheta(*r, theta);
		//printf("dS_dtheta = %f ", dS_dtheta);
#else
		dS_dtheta = get_dS_dtheta(r0, theta);
		//printf("r0 = %f dS_dtheta = %f \n", r0, dS_dtheta);
#endif

		S += dS_dtheta * dtheta;
		//printf("S = %f S0 = %f\n", S, S0);
		DBG_INFO("dS_dtheta = %f ", dS_dtheta);
		if (0.0 != R_lw_zap){
			*result += dS_dtheta / R_lw_zap * dtheta;
			DBG_INFO("result = %f ", *result);
		}
		else
		{
			ommited_S += dS_dtheta * dtheta;
			DBG_INFO("ommited_S = %e ", ommited_S);
		}

		DBG_INFO("\n");
	}
	DBG_INFO("result = %f\n", *result);
	//printf("S = %f ommited_S = %f S0 = %f\n", S, ommited_S, S0);
	if (0.0 != ommited_S)
	{
		S -= ommited_S;
		DBG_INFO("corrected S = %f\n", S);
	}
	double sigma = q / S;
	DBG_INFO("sigma0 = %f\n", sigma0);
	DBG_INFO("sigma = %f\n", sigma);
	*result *= sigma;
	DBG_INFO("result = %f\n", *result);
	return error;
}

int integral_phi_and_E(double q, double t, double R0, double r0, double v0, double a0, double * pE_minus_grad_phi_R0, double *pE_minus_1_c_dA_dt_R0, double r_min, double *phi, double *r)
{
	int err, error = 0;
	int i;
	double theta,
		t_zap, r_zap, R_zap,
		R_lw_zap, dS_dtheta;

	double v_zap;
	double a_zap;
	double aR_zap;
	double cos_alpha_zap;
	double E_minus_grad_varphi_R0;
	double E_minus_1_c_dA_dt_R0;

	int N = 1000;
	double dtheta = Pi / N;
	*phi = 0.0;
	double sigma0 = get_sigma(q, r0);
	double ommited_S = 0.0;
	double S0 = 4*Pi*r0*r0;
	double S = 0.0;
	printf("integral_phi_and_E(q=%f t=%f, R0=%f, r0=%f, v0=%f, a0=%f)\n", q, t, R0, r0, v0, a0);

	*pE_minus_grad_phi_R0 = 0.0;
	*pE_minus_1_c_dA_dt_R0 = 0.0;

#ifdef ALGORITHM_VERSION_0
	err = get_r(q, t, r0, v0, a0, r_min, r);
#endif
#ifdef ALGORITHM_VERSION_1
	err = get_r_ex1(q, t, r0, v0, r_min, r);
#endif
#ifdef ALGORITHM_VERSION_2
	err = get_r_ex2(q, t, r0, v0, r_min, r);
#endif
	if (0 != err)
	{
		error += 1;
	}
	printf("r = %f err = %d ", *r, err);

	for (i = 0; i <= N; ++i)
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
		DBG_INFO("v_zap = %f ", v_zap);
		DBG_INFO("a_zap = %f ", a_zap);

		aR_zap = get_aR_zap(R0, theta, r_zap, a_zap);
		DBG_INFO("aR_zap = %f ", aR_zap);
		cos_alpha_zap = get_cos_alpha_zap(R0, theta, r_zap, R_zap);
		DBG_INFO("cos_alpha_zap = %f ", cos_alpha_zap);
		E_minus_grad_varphi_R0 = get_E_minus_grad_phi_R0 (theta, v_zap, R_zap, aR_zap, R_lw_zap, cos_alpha_zap);
		E_minus_1_c_dA_dt_R0   = get_E_minus_1_c_dA_dt_R0(theta, v_zap, a_zap, R_zap, aR_zap, R_lw_zap);

#ifdef OLD_DS_THETA_ALG
		dS_dtheta = get_dS_dtheta(*r, theta);
		//printf("dS_dtheta = %f ", dS_dtheta);
#else
		dS_dtheta = get_dS_dtheta(r0, theta);
		//printf("r0 = %f dS_dtheta = %f ", r0, dS_dtheta);
#endif
		S += dS_dtheta * dtheta;
		//printf("S = %f S0 = %f\n", S, S0);
		DBG_INFO("dS_dtheta = %f ", dS_dtheta);
		if (0.0 != R_lw_zap){
			*phi                   += dS_dtheta * dtheta / R_lw_zap ;
			*pE_minus_grad_phi_R0  += dS_dtheta * dtheta * E_minus_grad_varphi_R0;
			*pE_minus_1_c_dA_dt_R0 += dS_dtheta * dtheta * E_minus_1_c_dA_dt_R0;
			DBG_INFO("phi = %f ", *phi);
			DBG_INFO("E1 = %f ", *pE_minus_grad_phi_R0);
			DBG_INFO("E2 = %f ", *pE_minus_1_c_dA_dt_R0);
		}
		else
		{
			ommited_S += dS_dtheta * dtheta;
			DBG_INFO("ommited_S = %e ", ommited_S);
		}

		DBG_INFO("\n");
	}
	DBG_INFO("phi = %f E1 = %f E2 = %f\n", *phi, *pE_minus_grad_phi_R0, *pE_minus_1_c_dA_dt_R0);
	DBG_INFO("sigma0 = %f\n", sigma0);
	//printf("S = %f ommited_S = %f S0 = %f\n", S, ommited_S, S0);
	if (0.0 != ommited_S)
	{
		S -= ommited_S;
		DBG_INFO("corrected S = %f\n", S);
	}
	double sigma = q / S;
	DBG_INFO("sigma = %f\n", sigma);
	*phi                   *= sigma;
	*pE_minus_grad_phi_R0  *= sigma;
	*pE_minus_1_c_dA_dt_R0 *= sigma;

	//printf("phi = %f E1 = %f E2 = %f\n", *phi, *pE_minus_grad_phi_R0, *pE_minus_1_c_dA_dt_R0);
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
	double r_zap, R_zap, R_lw_zap, phi_lw;
	double E_minus_grad_phi_R0, E_minus_1_c_dA_dt_R0;

	/* Текущий момент */
	double t = 5;
	/* расстояние от центра сферы к точке наблюдения. Точка наблюдения расположена на оси z в сферической системе координат */
	double R0 = 2;
	/* начальный радиус заряженной сферы в момент t=t_start */
	double r0 = 1;
	/* скорость в момент t_start*/
	double v0 = 0.0;
	/* ускорение */
	double a0 = 0.0;
	/* минимально возможный радиус заряженной сферы (из соображений упругости) в момент более ранний чем t=t_start */
	double r_min = 0.1;
	/* угловая координата заряда на заряженной сфере в сферической системе координат */
	double theta = Pi/2;
	/* Заряд сферы */
	double q = 1.0;
	double t_zap;
	double r;

	error = integral_phi(/*q*/-1, /*t*/0, /*R0*/2.0, /*r0*/2.0, /*v0*/c/3.0, /*a0*/0.0, /*r_min*/0.1, &phi_lw, &r);
	printf("\nphi_lw = %0.10f error = %d r = %f\n", phi_lw, error, r);
	#ifdef CALC_LW_WITHOUT_LAGGING
	printf ("should be -0.4315231087\n\n");
	#else
	printf ("should be -0.5198603854\n\n");
	#endif

	error = integral_phi(/*q*/-1, /*t*/0, /*R0*/2.0, /*r0*/2.0, /*v0*/0.0, /*a0*/0.0, /*r_min*/0.1, &phi_lw, &r);
	printf("\nphi_lw = %0.10f error = %d r = %f\n", phi_lw, error, r);
	#ifdef CALC_LW_WITHOUT_LAGGING
	printf ("should be -0.5\n\n");
	#else
	printf ("should be -0.5\n\n");
	#endif

	error = calc_tzap(q, /*t*/0.0, /*R0*/0.0, /*r0*/2.0, /*v0*/1.0, /*a0*/0.0, /*theta*/0.0, r_min, &t_zap);
	printf("t_zap = %f\n", t_zap);

	error = integral_phi(/*q*/-1, /*t*/0, /*R0*/0.0, /*r0*/2.0, /*v0*/c/3.0, /*a0*/0.0, /*r_min*/0.1, &phi_lw, &r);
	printf("\nphi_lw = %0.10f error = %d r = %f\n", phi_lw, error, r);
	#ifdef CALC_LW_WITHOUT_LAGGING
	printf ("should be -0.375\n\n");
	#else
	printf ("should be -0.5\n\n");
	#endif

	error = integral_phi(/*q*/-1, /*t*/0, /*R0*/1.0, /*r0*/2.0, /*v0*/c/3.0, /*a0*/0.0, /*r_min*/0.1, &phi_lw, &r);
	printf("\nphi_lw = %0.10f error = %d r = %f\n", phi_lw, error, r);
	#ifdef CALC_LW_WITHOUT_LAGGING
	printf ("should be -0.38348\n\n");
	#else
	printf ("should be -0.5047083549\n\n");
	#endif

	error = integral_phi(/*q*/-1, /*t*/0, /*R0*/3.0, /*r0*/2.0, /*v0*/c/3.0, /*a0*/0.0, /*r_min*/0.1, &phi_lw, &r);
	printf("\nphi_lw = %0.10f error = %d r = %f\n", phi_lw, error, r);
	#ifdef CALC_LW_WITHOUT_LAGGING
	printf ("should be -0.31720\n\n");
	#else
	printf ("should be -0.3465735903\n\n");
	#endif


	error = integral_phi(/*q*/-1, /*t*/0, /*R0*/0.0, /*r0*/2.0, /*v0*/1.0*c/3.0, /*a0*/0.0, /*r_min*/0.1, &phi_lw, &r);
	printf("\nphi_lw = %0.10f error = %d r = %f\n", phi_lw, error, r);
#ifndef CALC_LW_WITHOUT_LAGGING
	printf("should be -0.5\n\n");
#endif

	error = integral_phi(/*q*/-1, /*t*/0, /*R0*/0.0, /*r0*/2.0*c/3.0, /*v0*/0.0, /*a0*/0.1*c/3.0, /*r_min*/0.1, &phi_lw, &r);
	printf("\nphi_lw = %0.10f error = %d r = %f\n", phi_lw, error, r);
	#ifndef CALC_LW_WITHOUT_LAGGING
	printf ("should be -0.50576 tzap should be -0.67424\n\n");
	#endif

	error = integral_phi(/*q*/-1, /*t*/0, /*R0*/0.0, /*r0*/1.0*c / 3.0, /*v0*/0.0, /*a0*/0.1*c / 3.0, /*r_min*/0.1, &phi_lw, &r);
	printf("\nphi_lw = %0.10f error = %d r = %f\n", phi_lw, error, r);
#ifndef CALC_LW_WITHOUT_LAGGING
	printf("should be -1.0057 tzap should be -0.33521\n\n");
#endif

	error = integral_phi(/*q*/-1, /*t*/0, /*R0*/1.0*c/3.0, /*r0*/1.0*c / 3.0, /*v0*/0.0, /*a0*/0.1*c / 3.0, /*r_min*/0.1, &phi_lw, &r);
	printf("\nphi_lw = %0.10f error = %d r = %f\n", phi_lw, error, r);
#ifndef CALC_LW_WITHOUT_LAGGING
	printf("should be -0.98892 tzap should be 2.*sqrt(445.+5.*cos(theta)-5.*sqrt(cos(theta)^2+180.*cos(theta)+7919.))\n\n");
#endif

	error = integral_phi(/*q*/-1, /*t*/0, /*R0*/2.0*c / 3.0, /*r0*/1.0*c / 3.0, /*v0*/0.0, /*a0*/0.1*c / 3.0, /*r_min*/0.1, &phi_lw, &r);
	printf("\nphi_lw = %0.10f error = %d r = %f\n", phi_lw, error, r);
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
	error = integral_phi(-1, 5.5, 2.8421709430404007e-13, 2, 0, -0.150000000000000, r_min, &phi_lw, &r);
	printf("phi_lw = %f error = %d\n", phi_lw, error);

	for (double t = 0; t < 8.0; t += 0.1)
	{
		error = integral_phi(-1, t, 2.8421709430404007e-13, 2, 0, -0.150000000000000, r_min, &phi_lw, &r);
		printf("t = %f phi_lw = %f error = %d\n", t, phi_lw, error);
	}

	error = integral_phi(-1, 6.0, 2.8421709430404007e-13, 2, 0, -0.150000000000000, r_min, &phi_lw, &r);
	printf("phi_lw = %f\n", phi_lw);
	error = integral_phi(-1, 6.5, 2.8421709430404007e-13, 2, 0, -0.150000000000000, r_min, &phi_lw, &r);
	printf("phi_lw = %f\n", phi_lw);
	error = integral_phi(-1, 7.0, 2.8421709430404007e-13, 2, 0, -0.150000000000000, r_min, &phi_lw, &r);
	printf("phi_lw = %f\n", phi_lw);
	error = integral_phi(-1, 7.5, 2.8421709430404007e-13, 2, 0, -0.150000000000000, r_min, &phi_lw, &r);
	printf("phi_lw = %f\n", phi_lw);

	get_R(1.5, 1.499999999999999, 0.0);
	get_R(1.5, 1.4999999979313956, 0.0);

	error = calc_R_lw(q, t, R0, r0, v0, a0, theta, &t_zap, &r_zap, &R_zap, r_min, &R_lw_zap);

	/* скалярный потенциал Лиенара Вихерта зарядов равномерно распределённых по сферической поверхности радиуса r0 и движущихся из центра.  */
	/* varphi := proc (q, t, r__0, v__0, a__0, R__0) options operator, arrow;
	int(2*Pi*r(t, r__0, v__0, a__0)^2*sin(theta)*sigma(q, r__0)/K__zap(tzap(t, r__0, v__0, a__0, R__0, theta), r__0, v__0, a__0, R__0, theta), theta = 0 .. Pi) end proc;*/
	error = integral_phi(q, t, R0, r0, v0, a0, r_min, &phi_lw, &r);
	printf("phi_lw = %f\n", phi_lw);
	error = integral_phi_and_E(q, t, R0, r0, v0, a0, &E_minus_grad_phi_R0, &E_minus_1_c_dA_dt_R0, r_min, &phi_lw, &r);
	printf("phi_lw = %f E1=%f E2 = %f\n", phi_lw, E_minus_grad_phi_R0, E_minus_1_c_dA_dt_R0);

	/*infinity error result*/
	error = integral_phi(q, t, -1.0000000000000142, 1, 0, 0, r_min, &phi_lw, &r);
	printf("phi_lw = %f\n", phi_lw);
	phi_lw = integral_phi_and_E(q, t, -1.0000000000000142, 1, 0, 0, &E_minus_grad_phi_R0, &E_minus_1_c_dA_dt_R0, r_min, &phi_lw, &r);
	printf("phi_lw = %f E1=%f E2 = %f\n", phi_lw, E_minus_grad_phi_R0, E_minus_1_c_dA_dt_R0);

	/*hung
	calc_tzap(t=5.000000, R0=-1.500000, r0=2.000000, v0=0.0, a0=1.000000, theta=0.000000
	*/
	error = calc_tzap(q, /*t*/5.0, /*R0*/-1.5, /*r0*/2.0, /*v0*/0.0, /*a0*/1.0, /*theta*/0.0, r_min, &t_zap);
	printf("t_zap = %f\n", t_zap);

	/*hung

	integral_phi(q=-1.000000 t=5.000000, R0=-1.500000, r0=2.000000, v0=0.0, a0=1.000000)

	t2=1.000667 t1=0.999333 t=5.000000 R=3.999333 dR=1.333710e-03 dR_pre=-1.333710e-03
	t2=0.999333 t1=1.000667 t=5.000000 R=4.000667 dR=-1.333710e-03 dR_pre=1.333710e-03
	t2=1.000667 t1=0.999333 t=5.000000 R=3.999333 dR=1.333710e-03 dR_pre=-1.333710e-03
	t2=0.999333 t1=1.000667 t=5.000000 R=4.000667 dR=-1.333710e-03 dR_pre=1.333710e-03
	*/

	error = integral_phi(-q, 5.0, -1.5, 2.0, 0.0, 1.0, r_min, &phi_lw, &r);
	printf("phi_lw = %f\n", phi_lw);
	error = integral_phi_and_E(-q, 5.0, -1.5, 2.0, 0.0, 1.0, &E_minus_grad_phi_R0, &E_minus_1_c_dA_dt_R0, r_min, &phi_lw, &r);
	printf("phi_lw = %f E1=%f E2 = %f\n", phi_lw, E_minus_grad_phi_R0, E_minus_1_c_dA_dt_R0);

	{
		double E1_p, E1_n, E2_p, E2_n, phi_p, phi_n;
		double v0_p = 0.0;
		double v0_n = 0.0;
		double a0_p = 0.0001;
		double a0_n = 0.01;
		double E_p, E_n, E1, E2, E;
		double q_add;

		for (double ti = 0.0; ti <= 15.0; ti += 0.1)
		//for (double R0_i = r0+0.5; R0_i < 20.0; R0_i += 0.5)
		{
			for (double R0_i = r0+0.5; R0_i < 20.0; R0_i += 0.5)
			//for (double ti = 0.0; ti <= 15.0; ti += 0.1)
			{
				printf("ti=%03.1f R0_i=%03.1f ", ti, R0_i);

				error = integral_phi_and_E(+q, ti, R0_i, r0, v0_p, a0_p, &E1_p, &E2_p, r_min, &phi_p, &r);
				error = integral_phi_and_E(-q, ti, R0_i, r0, v0_n, a0_n, &E1_n, &E2_n, r_min, &phi_n, &r);
				//printf("phi_p = %f ", phi_p);
				//printf("phi_n = %f ", phi_n);

				printf("E1_p=%e E2_p=%e ", E1_p, E2_p);
				printf("E1_n=%e E2_n=%e ", E1_n, E2_n);

				E_p = E1_p + E2_p;
				E_n = E1_n + E2_n;

				printf("E_p=%e ", E_p);
				printf("E_n=%e ", E_n);

				E1 = E1_p + E1_n;
				E2 = E2_p + E2_n;

				printf("E1=%e ", E1);
				printf("E2=%e ", E2);

				E = E_p + E_n;

				printf("E=%f\n", E);

				q_add = E * 4 * Pi*R0_i*R0_i;
				printf("q_add=%f\n", q_add);
			}
			printf("\n");
		}
		error = integral_phi_and_E(q, t, R0, r0, 0.0, 0.0, &E_minus_grad_phi_R0, &E_minus_1_c_dA_dt_R0, r_min, &phi_lw, &r);
		printf("phi_lw = %f E1=%f E2 = %f\n", phi_lw, E_minus_grad_phi_R0, E_minus_1_c_dA_dt_R0);

		error = integral_phi_and_E(-q, t, R0, r0, 0.0, 0.0, &E_minus_grad_phi_R0, &E_minus_1_c_dA_dt_R0, r_min, &phi_lw, &r);
		printf("phi_lw = %f E1=%f E2 = %f\n", phi_lw, E_minus_grad_phi_R0, E_minus_1_c_dA_dt_R0);

	}

	return 0;
}

extern int v_n_t; // итератор полноты заполения двумерных массивов по координате времени

int test_v1()
{
	/* начальный радиус заряженной сферы в момент t=t_start */
	double r0_pos = 1.0;
	double r0_neg = 1.0;

	/* минимально возможный радиус заряженной сферы (из соображений упругости) в момент более ранний чем t=t_start */
	double r_min = 0.1;

	/* скорость в момент t_start*/
	velocity v0_pos = 0.0;
	velocity v0_neg = 0.0;

	/* ускорение вызванное причинами неэлектрического характера, например вледствие подвода энергии извне */
	double a0_neg = 1.0*c;
	double a0_pos = 0.1*c;

	double m_pos = 10.0;
	double m_neg = 1.0;

	/* Заряд сферы */
	double q = 1.0;

	int error = 0;
	double phi_lw_pos;
	double phi_lw_neg;
	double E_minus_grad_phi_R0_pos, E_minus_1_c_dA_dt_R0_pos;
	double E_minus_grad_phi_R0_neg, E_minus_1_c_dA_dt_R0_neg;
	double r_pos;
	double r_neg;
	double E1, E2, E;

	double dr = get_dr();
	double dt = get_dt();
	printf ("dr = %f ", dr);
	printf ("dt = %f ", dt);

	double t_a0 = dt * get_dt() / 2;

	init_array_1(a0_pos, v0_pos, r0_pos, a0_neg, v0_neg, r0_neg);

	for (v_n_t = 0; v_n_t < get_nt(); ++v_n_t)
	{
		double t = v_n_t * get_dt();
		for (int v_n_r = 0; v_n_r < get_nr(); ++v_n_r)
		{
			double R0 = v_n_r * dr;
			printf ("R0 = %f v_n_r = %d dr = %f\n", R0, v_n_r, dr);

			error = integral_phi_and_E(+q, t, R0, r0_pos, v0_pos, a0_pos, &E_minus_grad_phi_R0_pos, &E_minus_1_c_dA_dt_R0_pos, r_min, &phi_lw_pos, &r_pos);
			printf("pos err = %d phi_lw = %f E1=%f E2 = %f\n", error, phi_lw_pos, E_minus_grad_phi_R0_pos, E_minus_1_c_dA_dt_R0_pos);

			error = integral_phi_and_E(-q, t, R0, r0_neg, v0_neg, a0_neg, &E_minus_grad_phi_R0_neg, &E_minus_1_c_dA_dt_R0_neg, r_min, &phi_lw_neg, &r_neg);
			printf("neg err = %d phi_lw = %f E1=%f E2 = %f\n", error, phi_lw_neg, E_minus_grad_phi_R0_neg, E_minus_1_c_dA_dt_R0_neg);

			E1 = E_minus_grad_phi_R0_pos + E_minus_grad_phi_R0_neg;
			E2 = E_minus_1_c_dA_dt_R0_pos + E_minus_1_c_dA_dt_R0_neg;

			E = E1 + E2;

			set_E_ex1(t, R0, E);
			set_E_ex_1(v_n_t, v_n_r, E);
		}

		for (int v_n_r = 0; v_n_r < get_nr(); ++v_n_r)
		{
			double R0 = v_n_r * dr;
			error = get_r_ex1(+q, t, r0_pos, v0_pos, r_min, &r_pos);
			error = get_r_ex1(-q, t, r0_neg, v0_neg, r_min, &r_neg);

			set_a_ex1(t, r_pos, a0_pos, t_a0, +q, m_pos);
			set_a_ex1(t, r_neg, a0_neg, t_a0, -q, m_neg);

			set_v_ex1(t, v0_pos, a0_pos, t_a0, +q, m_pos);
			set_v_ex1(t, v0_neg, a0_neg, t_a0, -q, m_neg);

			set_s_ex1(t, v0_pos, +q);
			set_s_ex1(t, v0_neg, -q);
		}
	}
}

int main()
{
#ifdef ALGORITHM_VERSION_0
	return test_v0();
#endif
#ifdef ALGORITHM_VERSION_1
	return test_v1();
#endif
}


