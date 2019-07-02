#include <stdio.h>
#define _USE_MATH_DEFINES
#include <math.h>
#define Pi M_PI
#include "tzap.h"

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

/* скалярный потенциал Лиенара Вихерта зарядов  */
/* varphi := proc (q, t, r__0, v__0, a__0, R__0) options operator, arrow; 
int(2*Pi*r(t, r__0, v__0, a__0)^2*sin(theta)*sigma(q, r__0)/K__zap(tzap(t, r__0, v__0, a__0, R__0, theta), r__0, v__0, a__0, R__0, theta), theta = 0 .. Pi) end proc;*/
double integral_phi(double q, double t, double R0, double r0, double a0)
{
	int i;
	double theta, t_zap, r_zap, R_zap, R_LW_zap, r, dS_dtheta;
	int N = 1000;
	double dtheta = Pi / N;
	double result = 0.0;
	double sigma = get_sigma(q, r0);
	
	for (i = 0; i <= N; ++i)
	{
		theta = (i * dtheta);
		/* численный расчёта запаздывающего момента */
		t_zap = calc_tzap(t, R0, r0, a0, theta);
		printf("t_zap = %f\n", t_zap );
		/* Запаздывающий радиус в зависимости от текущего момента */
		r_zap = get_r(t_zap, r0, a0); /* расстояние от заряда до центра сферы в запаздывающий момент времени */
		R_zap = get_R(R0, r_zap, theta); /* расстояние от заряда до точки наблюдения в запаздывающий момент времени */

		/* Радиус Лиенара Вихерта */
		R_LW_zap = get_R(R0, r_zap, theta) - (get_v(t_zap, a0) / c) * (R0 * cos(theta) - r_zap);
		
		r = get_r(t, r0, a0);
		dS_dtheta = get_dS_dtheta(r, theta);
		result += dS_dtheta / R_LW_zap;
	}
	result *= sigma;
	return result;
}

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

	return 0;
}