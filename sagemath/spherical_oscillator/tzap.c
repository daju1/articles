#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include "tzap.h"
#include "dbg_info.h"
 

double t_start = 0;
double c = 1.0;
	
/* ускорение заряда */
double get_a(double t_zap, double a0)
{
	if (t_zap < t_start)
		return 0;
	return a0;
}

/* радиальная скорость заряда */
double get_v(double t_zap, double a0)
{
	if (t_zap < t_start)
		return 0;
	return a0*t_zap;
}
/* перемещение заряда */
double get_s(double t_zap, double a0)
{
	if (t_zap < t_start)
	{
		/*DBG_INFO("get_s t_zap=%f a0=%f returns 0\n", t_zap, a0);*/
		return 0;
	}
	/*DBG_INFO("get_s t_zap=%f a0=%f returns %f\n", t_zap, a0, a0*t_zap*t_zap/2);*/
	return a0*t_zap*t_zap/2;
}

/* расстояние от заряда до центра сферы в запаздывающий момент времени */
double get_r(double t_zap, double r0, double a0)
{
	return r0+get_s(t_zap, a0);
}

/* расстояние от заряда до точки наблюдения в запаздывающий момент времени */
double get_R(double R0, double r, double theta)
{
	double R = sqrt(R0*R0 - 2*R0*r*cos(theta) + r*r);
	return R;
}

/* численный расчёта запаздывающего момента */
double calc_tzap(double t, double R0, double r0, double a0, double theta)
{
	double epsilon = 1.0e-15;
	double t1;
	double t2 = t;
	double r, R, R_pre = DBL_MAX;
	double dR, dR_pre = DBL_MAX;
	double R_tmp;

	int i = 0;
	double n = 0.9;

	DBG_INFO("calc_tzap(t=%f, R0=%f, r0=%f, a0=%f, theta=%f\n", t, R0, r0, a0, theta);
	/*
	DBG_INFO("epsilon=%e\n", epsilon);
	DBG_INFO("t1=%f\n", t1);
	DBG_INFO("t2=%f\n", t2);
	*/

	do
	{
		t1 = t2;                 /* итерационный момент времени - на первой итерации текущее время наблюдения */
		r = get_r(t1, r0, a0);   /* итерационная координата заряда                                                                                    */
		R = get_R(R0, r, theta); /* итерационный радиус - на первой итерации текущий радиус                                      */
		t2 = t -  R / c;         /* время прохождения сигнала от итерационной координаты в точку наблюдения      */
		dR = c*(t-t1) - R;       /**/

		DBG_INFO("t2=%f t1=%f t=%f R=%f dR=%e dR_pre=%e ", t2, t1, t, R, dR, dR_pre);

		if (i > 1 && fabs(dR) - fabs(dR_pre) < 1.0e-3)
		{
			int j = 0;

			do
			{
				R_tmp = R_pre + n * (R - R_pre);
				DBG_INFO("R_tmp = R_pre + n * (R - R_pre);= %f ", R_tmp);

				t2 = t -  R_tmp / c;
				dR = c*(t-t1) - R_tmp;
				DBG_INFO("t2 = %f ", t2);

				n *= 0.9;
				DBG_INFO("n = %f ", n);
				++j;
			}
			while (fabs(dR) > fabs(dR_pre));
	
			R = R_tmp;
		}

		dR_pre = dR;
		R_pre = R;
		
		DBG_INFO("\n");
		++i;
	}
	while (fabs(t1 - t2) > epsilon);
	 
	DBG_INFO("fabs(t1 - t2) = %e fabs(t - t2) = %e calc_tzap() result=%f\n", fabs(t1 - t2), fabs(t - t2), t2); 
	return t2;
}
 
float calc_tzap_float(float t, float R0, float r0, float a0, float theta)
{
	return calc_tzap(t, R0, r0, a0, theta);
}
