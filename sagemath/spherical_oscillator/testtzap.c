#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define Pi M_PI
#include "tzap.h"

//#define USE_DEBUG
#include "dbg_info.h"

extern double c;

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
double calc_R_lw(double t, double R0, double r0, double a0, double theta, double * pt_zap, double * pr_zap, double * pR_zap)
{
	double R_lw_zap;
	/* численный расчёта запаздывающего момента */
	*pt_zap = calc_tzap(t, R0, r0, a0, theta);
	DBG_INFO("theta = %f t_zap = %f ", theta, t_zap);
	/* Запаздывающий радиус в зависимости от текущего момента */
	*pr_zap = get_r(*pt_zap, r0, a0); /* расстояние от заряда до центра сферы в запаздывающий момент времени */
	*pR_zap = get_R(R0, *pr_zap, theta); /* расстояние от заряда до точки наблюдения в запаздывающий момент времени */
	DBG_INFO("r_zap = %f ", r_zap);
	DBG_INFO("R_zap = %f ", R_zap);
	/* Радиус Лиенара Вихерта */
	R_lw_zap = get_R(R0, *pr_zap, theta) - (get_v(*pt_zap, a0) / c) * (R0 * cos(theta) - *pr_zap);
	DBG_INFO("R_lw_zap = %f ", R_lw_zap);

	return R_lw_zap;
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

/* Скалярное произведение ускорения частицы в запаздывающий момент времени на запаздывающий радиус-вектором (вектор из запаздывающего положения заряда в точку наблюдения) */
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

double get_E_minus_grad_varphi_R0(double q, double r0, double theta, double v_zap, double R_zap, double aR_zap, double R_lw_zap, double cos_alpha_zap)
{
	double E_minus_grad_varphi_R0 = 
		get_sigma(q, r0) * 
		(
			(cos_alpha_zap * R_zap / R_lw_zap) * (1.0 + (aR_zap - v_zap * v_zap) / (c * 2) )
			- v_zap*cos(theta) / c 
		)
		/ (R_lw_zap * R_lw_zap);
	return E_minus_grad_varphi_R0;
}

double calc_E_minus_grad_varphi_R0(double t, double R0, double r0, double a0, double theta, double * pt_zap, double * pr_zap, double * pR_zap)
{

}


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

double get_E_minus_1_c_dA_dt_R0(double q, double r0, double theta, double v_zap, double a_zap, double R_zap, double aR_zap, double R_lw_zap, double cos_alpha_zap)
{
	double E_minus_1_c_dA_dt_R0 =
		cos(theta)*get_sigma(q, r0) *
		( 
			(v_zap / (c * c)) * ( (R_zap / R_lw_zap) * ( (v_zap * v_zap - aR_zap) / c - c)  + c)
			- a_zap * R_zap / (c * c)
		)
		/
		(R_lw_zap * R_lw_zap);
	return E_minus_1_c_dA_dt_R0;
}

/* скалярный потенциал Лиенара Вихерта зарядов равномерно распределённых по сферической поверхности радиуса r0 и движущихся из центра. */
/* 
varphi := proc (q, t, r__0, v__0, a__0, R__0) options operator, arrow; 
int(2*Pi*r(t, r__0, v__0, a__0)^2*sin(theta)*sigma(q, r__0)/K__zap(tzap(t, r__0, v__0, a__0, R__0, theta), r__0, v__0, a__0, R__0, theta), theta = 0 .. Pi) end proc;
*/
double integral_phi(double q, double t, double R0, double r0, double a0)
{
	int i;
	double theta, 
		t_zap, r_zap, R_zap, 
		R_lw_zap, r, dS_dtheta;
	int N = 1000;
	double dtheta = Pi / N;
	double result = 0.0;
	double sigma = get_sigma(q, r0);
	double ommited_S = 0.0;
	double S = 4*Pi*r0*r0;
	DBG_INFO("integral_phi(q=%f t=%f, R0=%f, r0=%f, a0=%f)\n", q, t, R0, r0, a0);
	
	for (i = 0; i <= N; ++i)
	{
		theta = (i * dtheta);
#if 0
		/* численный расчёта запаздывающего момента */
		t_zap = calc_tzap(t, R0, r0, a0, theta);
		DBG_INFO("theta = %f t_zap = %f ", theta, t_zap);
		/* Запаздывающий радиус в зависимости от текущего момента */
		r_zap = get_r(t_zap, r0, a0); /* расстояние от заряда до центра сферы в запаздывающий момент времени */
		R_zap = get_R(R0, r_zap, theta); /* расстояние от заряда до точки наблюдения в запаздывающий момент времени */
		DBG_INFO("r_zap = %f ", r_zap);
		DBG_INFO("R_zap = %f ", R_zap);
		/* Радиус Лиенара Вихерта */
		R_lw_zap = get_R(R0, r_zap, theta) - (get_v(t_zap, a0) / c) * (R0 * cos(theta) - r_zap);
		DBG_INFO("R_lw_zap = %f ", R_lw_zap);
#else
		R_lw_zap = calc_R_lw(t, R0, r0, a0, theta, &t_zap, &r_zap, &R_zap);
#endif
		r = get_r(t, r0, a0);
		DBG_INFO("r = %f ", r);
		dS_dtheta = get_dS_dtheta(r, theta);
		DBG_INFO("dS_dtheta = %f ", dS_dtheta);
		if (0.0 != R_lw_zap){
			result += dS_dtheta / R_lw_zap * dtheta;
			DBG_INFO("result = %f ", result);
		}
		else
		{
			ommited_S += dS_dtheta * dtheta;
			DBG_INFO("ommited_S = %e ", ommited_S);
		}

		DBG_INFO("\n", 0);
	}
	DBG_INFO("result = %f\n", result);
	DBG_INFO("sigma = %f\n", sigma);
	if (0.0 != ommited_S)
	{
		S -= ommited_S;
		sigma = q / S;
		DBG_INFO("corrected sigma = %f\n", sigma);
	}
	result *= sigma;
	DBG_INFO("result = %f\n", result);
	return result;
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

int main()
{
	double r_zap, R_zap, R_lw_zap, phi_lw;

	/* Текущий момент */
	double t = 5;
	/* расстояние от центра сферы к точке наблюдения. Точка наблюдения расположена на оси z в сферической системе координат */
	double R0 = 2;
	/* начальный радиус заряженной сферы в момент t=0 */
	double r0 = 1;
	/* ускорение */
	double a0 = 0.0;
	/* угловая координата заряда на заряженной сфере в сферической системе координат */
	double theta = Pi/2;
	/* Заряд сферы */
	double q = 1.0;
	double t_zap;
#if 0
	/* численный расчёта запаздывающего момента */
	t_zap = calc_tzap(t, R0, r0, a0, theta);
	printf("t_zap = %f\n", t_zap );
	/* Запаздывающий радиус в зависимости от текущего момента */
	r_zap = get_r(t_zap, r0, a0); /* расстояние от заряда до центра сферы в запаздывающий момент времени */
	R_zap = get_R(R0, r_zap, theta); /* расстояние от заряда до точки наблюдения в запаздывающий момент времени */

	printf("R_zap = %f c*(t-t_zap) = %f delta = %e\n", R_zap, c*(t-t_zap), R_zap - c*(t-t_zap));
	
	/* Радиус Лиенара Вихерта */
	R_lw_zap = get_R(R0, r_zap, theta) - (get_v(t_zap, a0) / c) * (R0 * cos(theta) - r_zap);
	printf("R_lw_zap = %f\n", R_lw_zap);
#else
	R_lw_zap = calc_R_lw(t, R0, r0, a0, theta, &t_zap, &r_zap, &R_zap);
#endif
	/* скалярный потенциал Лиенара Вихерта зарядов равномерно распределённых по сферической поверхности радиуса r0 и движущихся из центра.  */
	/* varphi := proc (q, t, r__0, v__0, a__0, R__0) options operator, arrow; 
	int(2*Pi*r(t, r__0, v__0, a__0)^2*sin(theta)*sigma(q, r__0)/K__zap(tzap(t, r__0, v__0, a__0, R__0, theta), r__0, v__0, a__0, R__0, theta), theta = 0 .. Pi) end proc;*/
	phi_lw = integral_phi(q, t, R0, r0, a0);
	printf("phi_lw = %f\n", phi_lw);

	/*infinity error result*/
	phi_lw = integral_phi(q, t, -1.0000000000000142, 1, 0);
	printf("phi_lw = %f\n", phi_lw);

	/*hung
	calc_tzap(t=5.000000, R0=-1.500000, r0=2.000000, a0=1.000000, theta=0.000000
	*/

	t_zap = calc_tzap(5.0, -1.5, 2.0, 1.0, 0.0);
	printf("t_zap = %f\n", t_zap);

	/*hung

	integral_phi(q=-1.000000 t=5.000000, R0=-1.500000, r0=2.000000, a0=1.000000)

	t2=1.000667 t1=0.999333 t=5.000000 R=3.999333 dR=1.333710e-03 dR_pre=-1.333710e-03 
	t2=0.999333 t1=1.000667 t=5.000000 R=4.000667 dR=-1.333710e-03 dR_pre=1.333710e-03 
	t2=1.000667 t1=0.999333 t=5.000000 R=3.999333 dR=1.333710e-03 dR_pre=-1.333710e-03 
	t2=0.999333 t1=1.000667 t=5.000000 R=4.000667 dR=-1.333710e-03 dR_pre=1.333710e-03 
	*/

	phi_lw = integral_phi(-q, 5.0, -1.5, 2.0, 1.0);
	printf("phi_lw = %f\n", phi_lw);

	return 0;
}

