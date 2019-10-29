#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include "tzap.h"
#include "dbg_info.h"
#include "stdlib.h"

#define T_START 0
#define T_FINISH 100
#define DT 0.0001



static const double t_start = T_START;    // момент включения ускорения
static const double dt = DT;              // шаг времени
static const double t_finish = T_FINISH; // момент времени окончания расчёта

// одномерные массивы для сохранения истории при интегрировании только лишь по времени
// без учёта эволючии объёмного распределения заряда - упрощённый случай сферического конденсатора
// сохраняется только лишь история положительной обкладки и отрицательной обкладки
static const int v_Nt = ((T_FINISH - T_START) / DT);
static int v_n = 0; // итератор полноты заполнения одномерных
static double epsilon_n = 1e-8;
static double * v_t;

static double * v_a_pos;
static double * v_v_pos;
static double * v_s_pos;
static double * v_r_pos;

static double * v_a_neg;
static double * v_v_neg;
static double * v_s_neg;
static double * v_r_neg;

// двумерные массивы для сохранения истории при интегрировани как по времени так и по r0 
// с учётом эволючии объёмного распределения заряда - случай объёмного взрыва плазмы
// сохраняется история положительного и отрицательного облака заряженных частиц

//static int v_N_t; // размер двумерного массива по координате времени
//static int v_N_r0; // размер двумерного массива по коодинате начального радиуса
static int v_n_t = 0; // итератор полноты заполения двумерных массивов по координате времени

static double ** vv_a_pos;
static double ** vv_v_pos;
static double ** vv_s_pos;
static double ** vv_r_pos;

static double ** vv_a_neg;
static double ** vv_v_neg;
static double ** vv_s_neg;
static double ** vv_r_neg;

static double v_max = 0.999 * c;

void init_array_1(double a0_pos, double v0_pos, double r0_pos, double a0_neg, double v0_neg, double r0_neg)
{
	v_t = malloc(v_Nt * sizeof(double *));

	v_a_pos = malloc(v_Nt * sizeof(double *));
	v_v_pos = malloc(v_Nt * sizeof(double *));
	v_s_pos = malloc(v_Nt * sizeof(double *));
	v_r_pos = malloc(v_Nt * sizeof(double *));

	v_a_neg = malloc(v_Nt * sizeof(double *));
	v_v_neg = malloc(v_Nt * sizeof(double *));
	v_s_neg = malloc(v_Nt * sizeof(double *));
	v_r_neg = malloc(v_Nt * sizeof(double *));

	v_a_pos[0] = a0_pos;
	v_v_pos[0] = v0_pos;
	v_s_pos[0] = 0.0;
	v_r_pos[0] = r0_pos;

	v_a_neg[0] = a0_neg;
	v_v_neg[0] = v0_neg;
	v_s_neg[0] = 0.0;
	v_r_neg[0] = r0_neg;
}


void init_array_2(int v_N_r0, int v_N_t, double a0_pos, double v0_pos, double r0_pos, double a0_neg, double v0_neg, double r0_neg)
{
	v_t = malloc(v_Nt * sizeof(double *));

	vv_a_pos = malloc(v_N_r0 * sizeof(double **));
	vv_v_pos = malloc(v_N_r0 * sizeof(double **));
	vv_s_pos = malloc(v_N_r0 * sizeof(double **));
	vv_r_pos = malloc(v_N_r0 * sizeof(double **));

	vv_a_neg = malloc(v_N_r0 * sizeof(double **));
	vv_v_neg = malloc(v_N_r0 * sizeof(double **));
	vv_s_neg = malloc(v_N_r0 * sizeof(double **));
	vv_r_neg = malloc(v_N_r0 * sizeof(double **));

	for (int i = 0; i < v_N_r0; ++i)
	{
		vv_a_pos[i] = malloc(v_N_t * sizeof(double *));
		vv_v_pos[i] = malloc(v_N_t * sizeof(double *));
		vv_s_pos[i] = malloc(v_N_t * sizeof(double *));
		vv_r_pos[i] = malloc(v_N_t * sizeof(double *));

		vv_a_pos[i] = malloc(v_N_t * sizeof(double *));
		vv_v_pos[i] = malloc(v_N_t * sizeof(double *));
		vv_s_pos[i] = malloc(v_N_t * sizeof(double *));
		vv_r_pos[i] = malloc(v_N_t * sizeof(double *));

		// TODO: rework initialization
		vv_a_pos[i][0] = a0_pos;
		vv_v_pos[i][0] = v0_pos;
		vv_s_pos[i][0] = 0.0;
		vv_r_pos[i][0] = r0_pos;

		vv_a_neg[i][0] = a0_neg;
		vv_v_neg[i][0] = v0_neg;
		vv_s_neg[i][0] = 0.0;
		vv_r_neg[i][0] = r0_neg;
	}
}
double get_c()
{
	return c;
}
	
/* ускорение заряда 
	a0   - ускорение сообщаемое заряду силами неэлектрического происхождения, например ускорение теплового разгона
	t_a0 - время действия сил неэлектрического происхождения
	E    - электрическое поле
	m    - масса заряда q

	для вычисления ускорения заряда вызванного электрическими силами
	F = m * a
	F = E * q
	m * a = E * q
	a = E * q / m
*/
double get_a_ex1(double t_zap, double a0, double t_a0, double E, double q, double m)
{
	double t_max;
	if (t_zap < t_start)
		return 0;
	// t_max = t_start + v_max / a0;
	// if (t_zap >= t_max)
	//     return 0;

	double n = (t_zap - t_start) / dt;
	double * v_a = q > 0 ? v_a_pos : v_a_neg;
	if (n <= (double)v_n)
	{
		// результат может быть взят с помощью линейной интерполяции ранее рассчитнных значений
		int n1 = floor(n);
		int n2 = ceil(n);
		double part = n - n1;

		double a = v_a[n1] + part * (v_a[n2] - v_a[n1]);
		return a;
	}

	if (n - v_n > 1.0 + epsilon_n)
	{
		assert(0);
	}

	double a = E * q / m;
	if (t_zap <= t_a0)
	{
		a += a0;
	}

	double part = n - v_n;
	double da = (a - v_a[v_n]) / part;
	v_a[v_n + 1] = v_a[v_n] + da;
	return a;
}

double get_v_ex1(double t_zap, double v0, double a0, double t_a0, double E, double q, double m)
{
	assert(v0 < v_max);
	double t_max;
	if (t_zap < t_start)
		return v0;

	double n = (t_zap - t_start) / dt;
	double * v_v = q > 0 ? v_v_pos : v_v_neg;
	if (n <= (double)v_n)
	{
		// результат может быть взят с помощью линейной интерполяции ранее рассчитнных значений
		int n1 = floor(n);
		int n2 = ceil(n);
		double part = n - n1;

		double a = v_v[n1] + part * (v_v[n2] - v_v[n1]);
		return a;
	}

	if (n - v_n > 1.0 + epsilon_n)
	{
		assert(0);
	}

	double a = get_a_ex1(t_zap, a0, t_a0, E, q, m);
	double v = v_v[v_n];
	// a = dv / dt
	// dv = a * dt
	double dv = a * dt;
	v += dv;

	double part = n - v_n;
	double dv_out = (v - v_v[v_n]) / part;
	v_v[v_n + 1] = v_v[v_n] + dv_out;
	return v;
}

double get_s_ex1(double t_zap, double v0, double a0, double t_a0, double E, double q, double m)
{
	assert(v0 < v_max);
	double dt_start, dt_max;
	double t_max;
	double s;
	if (t_zap < t_start)
	{
		s = v0*(t_zap - t_start);
		DBG_INFO("get_s1 t_zap=%f a0=%f returns %f\n", t_zap, a0, s);
		return s;
	}

	double n = (t_zap - t_start) / dt;
	double * v_s = q > 0 ? v_s_pos : v_s_neg;
	if (n <= (double)v_n)
	{
		// результат может быть взят с помощью линейной интерполяции ранее рассчитнных значений
		int n1 = floor(n);
		int n2 = ceil(n);
		double part = n - n1;

		double s = v_s[n1] + part * (v_s[n2] - v_s[n1]);
		return s;
	}

	if (n - v_n > 1.0 + epsilon_n)
	{
		assert(0);
	}

	double a = get_a_ex1(t_zap, a0, t_a0, E, q, m);
	double v = get_v_ex1(t_zap, v0, a0, t_a0, E, q, m);
	s = v_s[v_n];
	// v = ds / dt
	double ds = v * dt;// + a * dt * dt ?????
	s + ds;

	double part = n - v_n;
	double ds_out = (s - v_s[v_n]) / part;
	v_s[v_n + 1] = v_s[v_n] + ds_out;

	return s;
}


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
	return v0 + a0*(t_zap-t_start);
}
/* перемещение заряда */
double get_s(double t_zap, double v0, double a0)
{
	assert(v0 < v_max);
	double dt_start, dt_max;
	double t_max;
	double s;
	if (t_zap < t_start)
	{
		s = v0*(t_zap - t_start);
		DBG_INFO("get_s1 t_zap=%f a0=%f returns %f\n", t_zap, a0, s);
		return s;
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
		DBG_INFO("get_s t_max %f\n", t_max);
	
		dt_start = (t_max - t_start);
		dt_max = (t_zap - t_max);
		s = v0 * dt_start + a0*dt_start*dt_start/2 + dt_max*v_max;
		DBG_INFO("get_s2 t_zap=%f a0=%f returns %f\n", t_zap, a0, s);
		return s;
	}

	dt_start = (t_zap - t_start);
	s = v0 * dt_start + a0*dt_start*dt_start/2;
	DBG_INFO("get_s t_zap=%f a0=%f dt_start=%f returns %f\n", t_zap, a0, dt_start, s);
	return s;
}

/* расстояние от заряда до центра сферы в запаздывающий момент времени */
int get_r(double t_zap, double r0, double v0, double a0, double r_min, double * r)
{
	int error = 0;
	*r = r0 + get_s(t_zap, v0, a0);
	if (*r < r_min)
	{
		DBG_INFO("Warning: r %f < r_min %f\n", *r, r_min);
		*r = r_min;
		error = 1;
	}
	assert(*r > 0.0);
	return error;
}

/* расстояние от заряда до точки наблюдения в запаздывающий момент времени */
double get_R(double R0, double r, double theta)
{
	//double RR = R0*R0 - 2 * R0*r*cos(theta) + r*r;
	//double RR2 = R0*R0*cos(theta) - 2 * R0*r*cos(theta) + r*r*cos(theta) + R0*R0*(1 - cos(theta)) + r*r*(1 - cos(theta));
	double RR3 = (R0 - r)*(R0 - r)*cos(theta) + (R0*R0 + r*r)*(1 - cos(theta));
	//DBG_INFO("RR=%f, RR2 = %f, RR3 = %f, R0=%f, r=%f, theta=%f)\n", RR, RR2, RR3, R0, r, theta);
	assert(RR3 >= 0.0);
	double R = sqrt(RR3);
	return R;
}

/* численный расчёта запаздывающего момента */
int calc_tzap(double t, double R0, double r0, double v0, double a0, double theta, double r_min, double * t2)
{
	int err, error = 0;
#ifdef CALC_LW_WITHOUT_LAGGING
	*t2 = t_start;
#else
	double epsilon = 1.0e-15;
	double t1;
	* t2 = t;
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
		t1 = *t2;                 /* итерационный "текущий" момент времени - на первой итерации текущее время наблюдения */
		err = get_r(t1, r0, v0, a0, r_min, &r);   /* итерационная координата заряда                                                      */
		if (0 != err)
		{
			error += 1;
		}
		assert(r >= 0.0);
		R = get_R(R0, r, theta); /* итерационный радиус - на первой итерации текущий радиус                             */
		*t2 = t - R / c;          /* время прохождения сигнала от итерационной координаты в точку наблюдения             */
		                         /* итерационный "запаздывающий" момент времени t2                                      */
		dR = c*(t-t1) - R;       /**/
		v1 = get_v(t1, v0, a0);      /* скорость заряда в итерационный "текущий" момент времени t1                          */
		v2 = get_v(*t2, v0, a0);      /* скорость заряда в итерационный "запаздывающий" момент времени t2                    */

		DBG_INFO("t2=%f t1=%f t=%f v1 = %f, v2 = %f, v = %f, R=%f dR=%e dR_pre=%e ", *t2, t1, t, v1, v2, v, R, dR, dR_pre);
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

				*t2 = t -  R_tmp / c;
				dR = c*(t-t1) - R_tmp;
				DBG_INFO("t2 = %f ", *t2);

				n *= 0.9;
				DBG_INFO("n = %f ", n);
				++j;
			}
			while (fabs(dR) > fabs(dR_pre));
	
			R = R_tmp;
		}
#endif
		dt = t1 - *t2;
		DBG_INFO("dt=%e \n", dt);
		dR_pre = dR;
		R_pre = R;
		
		DBG_INFO("\n");
		++i;
	}
	while (fabs(dt) > epsilon);
	 
	DBG_INFO("fabs(t1 - t2) = %e fabs(t - t2) = %e calc_tzap() result=%f\n", fabs(t1 - *t2), fabs(t - *t2), *t2);
#endif
	return error;
}
 
/*float calc_tzap_float(float t, float R0, float r0, float v0, float a0, float theta)
{
	return calc_tzap(t, R0, r0, v0, a0, theta);
}*/
