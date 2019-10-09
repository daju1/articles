#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include "tzap.h"
#include "dbg_info.h"
 

double t_start = 0; // момент включения ускорения
double c = 1.0;
double v_max = 0.999;
	
/* ускорение заряда */
double get_a(double t_zap, double a0)
{
	double t_max;
	if (t_zap < t_start)
		return 0;
	t_max = t_start + v_max / a0;
	if (t_zap >= t_max)
		return 0;
	return a0;
}

/* радиальная скорость заряда */
double get_v(double t_zap, double v0, double a0)
{
	assert(v0 < v_max);
	double t_max;
	if (t_zap < t_start)
		return v0;
	t_max = t_start + v_max / a0;
	if (t_zap >= t_max)
		return v_max;
	return v0 + a0*t_zap;
}
/* перемещение заряда */
double get_s(double t_zap, double v0, double a0)
{
	assert(v0 < v_max);
	double dt_start, dt_max;
	double t_max;
	if (t_zap < t_start)
	{
		/*DBG_INFO("get_s t_zap=%f a0=%f returns 0\n", t_zap, a0);*/
		return v0*(t_zap - t_start);
	}

	if (a0 > 0.0)
	{
		t_max = t_start + (v_max - v0) / a0;
		assert(t_max > t_start);
	}
	else if (a0 < 0.0)
	{
		t_max = t_start + (-v_max + v0) / a0;
		assert(t_max > t_start);
	}

	if (a0 != 0.0 && t_zap >= t_max)
	{
		dt_start = (t_max - t_start);
		dt_max = (t_zap - t_max);
		return v0 * dt_start + a0*dt_start*dt_start/2 + dt_max*v_max;
	}

	dt_start = (t_zap - t_start);
	/*DBG_INFO("get_s t_zap=%f a0=%f returns %f\n", t_zap, a0, a0*t_zap*t_zap/2);*/
	return v0 * dt_start + a0*dt_start*dt_start/2;
}

/* расстояние от заряда до центра сферы в запаздывающий момент времени */
double get_r(double t_zap, double r0, double v0, double a0)
{
	return r0 + get_s(t_zap, v0, a0);
}

/* расстояние от заряда до точки наблюдения в запаздывающий момент времени */
double get_R(double R0, double r, double theta)
{
	double RR = R0*R0 - 2 * R0*r*cos(theta) + r*r;
	double RR2 = R0*R0*cos(theta) - 2 * R0*r*cos(theta) + r*r*cos(theta) + R0*R0*(1 - cos(theta)) + r*r*(1 - cos(theta));
	double RR3 = (R0 - r)*(R0 - r)*cos(theta) + (R0*R0 + r*r)*(1 - cos(theta));
	DBG_INFO("RR=%f, RR2 = %f, RR3 = %f, R0=%f, r=%f, theta=%f)\n", RR, RR2, RR3, R0, r, theta);
	assert(RR3 >= 0.0);
	double R = sqrt(RR3);
	return R;
}

/* численный расчёта запаздывающего момента */
double calc_tzap(double t, double R0, double r0, double v0, double a0, double theta)
{
	double epsilon = 1.0e-15;
	double t1;
	double t2 = t;
	double dt;
	double v1,v2,v;
	double r, R, R_pre = DBL_MAX;
	double dR, dR_pre = DBL_MAX;
	double R_tmp;

	int i = 0;
	double n = 0.9;
	v = get_v(t, v0, a0);      /* скорость заряда в текущий момент времени t                          */

	DBG_INFO("calc_tzap(t=%f, v = %f, R0=%f, r0=%f, v0=%f, a0=%f, theta=%f)\n", t, v, R0, r0, v0, a0, theta);
	assert(v < c);
	/*
	DBG_INFO("epsilon=%e\n", epsilon);
	DBG_INFO("t1=%f\n", t1);
	DBG_INFO("t2=%f\n", t2);
	*/

	do
	{
		t1 = t2;                 /* итерационный "текущий" момент времени - на первой итерации текущее время наблюдения */
		r = get_r(t1, r0, v0, a0);   /* итерационная координата заряда                                                      */
		R = get_R(R0, r, theta); /* итерационный радиус - на первой итерации текущий радиус                             */
		t2 = t - R / c;          /* время прохождения сигнала от итерационной координаты в точку наблюдения             */
		                         /* итерационный "запаздывающий" момент времени t2                                      */
		dR = c*(t-t1) - R;       /**/
		v1 = get_v(t1, v0, a0);      /* скорость заряда в итерационный "текущий" момент времени t1                          */
		v2 = get_v(t2, v0, a0);      /* скорость заряда в итерационный "запаздывающий" момент времени t2                    */

		DBG_INFO("t2=%f t1=%f t=%f v1 = %f, v2 = %f, v = %f, R=%f dR=%e dR_pre=%e ", t2, t1, t, v1, v2, v, R, dR, dR_pre);
		assert(v1 < c);
		assert(v2 < c);
#if 1
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
#endif
		dt = t1 - t2;
		DBG_INFO("dt=%e \n", dt);
		dR_pre = dR;
		R_pre = R;
		
		DBG_INFO("\n");
		++i;
	}
	while (fabs(dt) > epsilon);
	 
	DBG_INFO("fabs(t1 - t2) = %e fabs(t - t2) = %e calc_tzap() result=%f\n", fabs(t1 - t2), fabs(t - t2), t2); 
	return t2;
}
 
/*float calc_tzap_float(float t, float R0, float r0, float v0, float a0, float theta)
{
	return calc_tzap(t, R0, r0, v0, a0, theta);
}*/
