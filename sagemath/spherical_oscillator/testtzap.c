#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define Pi M_PI
#include "tzap.h"
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

/* скалярный потенциал Лиенара Вихерта зарядов равномерно распределённых по сферической поверхности радиуса r0 и движущихся из центра. */
/* 
varphi := proc (q, t, r__0, v__0, a__0, R__0) options operator, arrow; 
int(2*Pi*r(t, r__0, v__0, a__0)^2*sin(theta)*sigma(q, r__0)/K__zap(tzap(t, r__0, v__0, a__0, R__0, theta), r__0, v__0, a__0, R__0, theta), theta = 0 .. Pi) end proc;
*/
double integral_phi(double q, double t, double R0, double r0, double a0)
{
	int i;
	double theta, t_zap, r_zap, R_zap, R_lw_zap, r, dS_dtheta;
	int N = 1000;
	double dtheta = Pi / N;
	double result = 0.0;
	double sigma = get_sigma(q, r0);
	double ommited_S = 0.0;
	double S = 4*Pi*r0*r0;
	
	for (i = 0; i <= N; ++i)
	{
		theta = (i * dtheta);
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

		DBG_INFO("\n");
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

	/* численный расчёта запаздывающего момента */
	double t_zap = calc_tzap(t, R0, r0, a0, theta);
	printf("t_zap = %f\n", t_zap );
	/* Запаздывающий радиус в зависимости от текущего момента */
	r_zap = get_r(t_zap, r0, a0); /* расстояние от заряда до центра сферы в запаздывающий момент времени */
	R_zap = get_R(R0, r_zap, theta); /* расстояние от заряда до точки наблюдения в запаздывающий момент времени */

	printf("R_zap = %f c*(t-t_zap) = %f delta = %e\n", R_zap, c*(t-t_zap), R_zap - c*(t-t_zap));
	
	/* Радиус Лиенара Вихерта */
	R_lw_zap = get_R(R0, r_zap, theta) - (get_v(t_zap, a0) / c) * (R0 * cos(theta) - r_zap);
	printf("R_lw_zap = %f\n", R_lw_zap);
	
	/* скалярный потенциал Лиенара Вихерта зарядов равномерно распределённых по сферической поверхности радиуса r0 и движущихся из центра.  */
	/* varphi := proc (q, t, r__0, v__0, a__0, R__0) options operator, arrow; 
	int(2*Pi*r(t, r__0, v__0, a__0)^2*sin(theta)*sigma(q, r__0)/K__zap(tzap(t, r__0, v__0, a__0, R__0, theta), r__0, v__0, a__0, R__0, theta), theta = 0 .. Pi) end proc;*/
	phi_lw = integral_phi(q, t, R0, r0, a0);
	printf("phi_lw = %f\n", phi_lw);

	phi_lw = integral_phi(q, t, -1.0000000000000142, 1, 0);

	return 0;
}

