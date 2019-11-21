#include <math.h>
#include <stdio.h>
#include <limits.h>
#include <float.h>
#include <assert.h>
#include "tzap.h"
#include "dbg_info.h"
#include "stdlib.h"

#define T_START 0.0
#define T_FINISH 10.0
//#define DT 0.01
#define R_START 0.0
#define R_FINISH 10.0
#define DR 0.01

#ifdef SI
#define LIGHT_VELONCITY 299792458.0
double multiplier_E = LIGHT_VELONCITY * LIGHT_VELONCITY / 10000000.0;
#ifdef MASS_IN_GRAMM
double multiplier_a = LIGHT_VELONCITY * LIGHT_VELONCITY / 10000.0;
#else
double multiplier_a = LIGHT_VELONCITY * LIGHT_VELONCITY / 10000000.0;
#endif
#else
#define LIGHT_VELONCITY 3.0
double multiplier_E = 1.0;
double multiplier_a = 1.0;
#endif

double g_c = LIGHT_VELONCITY;
static double v_max = 0.999 * LIGHT_VELONCITY;

static double g_t_start = T_START;    // момент включения ускорения
//static double g_dt = DT;              // шаг времени
static double g_t_finish = T_FINISH;  // момент времени окончания расчёта
static double g_dr = DR;              // шаг координаты

// одномерные массивы для сохранения истории при интегрировании только лишь по времени
// без учёта эволючии объёмного распределения заряда - упрощённый случай сферического конденсатора
// сохраняется только лишь история положительной обкладки и отрицательной обкладки
static double g_r_start = R_START;
static double g_r_finish = R_FINISH;
//static int v_Nt = (int)((T_FINISH - T_START) / DT);
static int v_Nt = 10000;
static int v_Nr = (int)((R_FINISH - R_START) / DR);

double get_dr()
{
	return g_dr;
}

void set_dr(double dr)
{
	g_dr = dr;
	v_Nr = (int)((g_r_finish - g_r_start) / g_dr);
}

void set_r_finish(double r_finish)
{
	g_r_finish = r_finish;
	v_Nr = (int)((g_r_finish - g_r_start) / g_dr);
}

int get_nt()
{
	return v_Nt;
}

int get_nr()
{
	return v_Nr;
}

static double epsilon_n = 1e-8;
static double epsilon_r = 1e-8;

double * v_t;
int v_n_t = 0; // итератор полноты заполения двумерных массивов по координате времени

#ifdef ALGORITHM_VERSION_1
static double * v_a_pos;
static double * v_v_pos;
static double * v_s_pos;
static double * v_r_pos;

static double * v_a_neg;
static double * v_v_neg;
static double * v_s_neg;
static double * v_r_neg;

static double ** v_E1;
static double ** v_E2;
static double ** v_E;
#endif /* ALGORITHM_VERSION_1 */


// двумерные массивы для сохранения истории при интегрировани как по времени так и по r0
// с учётом эволючии объёмного распределения заряда - случай объёмного взрыва плазмы
// сохраняется история положительного и отрицательного облака заряженных частиц

//static int v_N_t; // размер двумерного массива по координате времени
//static int v_N_r; // размер двумерного массива по коодинате радиуса точки наблюдения
//static int v_N_r0; // размер двумерного массива по коодинате начального радиуса

#ifdef ALGORITHM_VERSION_2
static double ** vv_a_pos;
static double ** vv_v_pos;
static double ** vv_s_pos;
static double ** vv_r_pos;

static double ** vv_a_neg;
static double ** vv_v_neg;
static double ** vv_s_neg;
static double ** vv_r_neg;

static double ** vv_E1;
static double ** vv_E2;
static double ** vv_E;
#endif

#ifdef ALGORITHM_VERSION_1
void init_array_1(double a0_pos, velocity v0_pos, double r0_pos, double a0_neg, velocity v0_neg, double r0_neg)
{
	v_t = malloc(v_Nt * sizeof(double));

	v_a_pos = malloc(v_Nt * sizeof(double));
	v_v_pos = malloc(v_Nt * sizeof(double));
	v_s_pos = malloc(v_Nt * sizeof(double));
	v_r_pos = malloc(v_Nt * sizeof(double));

	v_a_neg = malloc(v_Nt * sizeof(double));
	v_v_neg = malloc(v_Nt * sizeof(double));
	v_s_neg = malloc(v_Nt * sizeof(double));
	v_r_neg = malloc(v_Nt * sizeof(double));

	v_E1 = malloc(v_Nr * sizeof(double *));
	v_E2 = malloc(v_Nr * sizeof(double *));
	v_E  = malloc(v_Nr * sizeof(double *));

	v_t[0] = T_START;

	v_a_pos[0] = a0_pos;
	v_v_pos[0] = v0_pos;
	v_s_pos[0] = 0.0;
	v_r_pos[0] = r0_pos;

	v_a_neg[0] = a0_neg;
	v_v_neg[0] = v0_neg;
	v_s_neg[0] = 0.0;
	v_r_neg[0] = r0_neg;

	for (int i_r = 0; i_r < v_Nr; ++i_r)
	{
		v_E1[i_r] = malloc(v_Nt * sizeof(double));
		v_E2[i_r] = malloc(v_Nt * sizeof(double));
		v_E [i_r] = malloc(v_Nt * sizeof(double));

		// initialization
		v_E1[i_r][0] = 0.0;
		v_E2[i_r][0] = 0.0;
		v_E [i_r][0] = 0.0;
	}
}
#endif /* ALGORITHM_VERSION_1 */

#ifdef ALGORITHM_VERSION_2
void init_array_2(int v_N_r0, int v_N_t, double * a0_pos, velocity * v0_pos, double * r0_pos, double * a0_neg, velocity * v0_neg, double * r0_neg)
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

	vv_E1 = malloc(v_N_r0 * sizeof(double **));
	vv_E2 = malloc(v_N_r0 * sizeof(double **));
	vv_E  = malloc(v_N_r0 * sizeof(double **));

	v_t[0] = T_START;

	for (int i_r0 = 0; i_r0 < v_N_r0; ++i_r0)
	{
		vv_a_pos[i_r0] = malloc(v_N_t * sizeof(double *));
		vv_v_pos[i_r0] = malloc(v_N_t * sizeof(double *));
		vv_s_pos[i_r0] = malloc(v_N_t * sizeof(double *));
		vv_r_pos[i_r0] = malloc(v_N_t * sizeof(double *));

		vv_a_pos[i_r0] = malloc(v_N_t * sizeof(double *));
		vv_v_pos[i_r0] = malloc(v_N_t * sizeof(double *));
		vv_s_pos[i_r0] = malloc(v_N_t * sizeof(double *));
		vv_r_pos[i_r0] = malloc(v_N_t * sizeof(double *));

		vv_E1[i_r0] = malloc(v_N_t * sizeof(double *));
		vv_E2[i_r0] = malloc(v_N_t * sizeof(double *));
		vv_E [i_r0] = malloc(v_N_t * sizeof(double *));

		// initialization
		vv_a_pos[i_r0][0] = a0_pos[i_r0];
		vv_v_pos[i_r0][0] = v0_pos[i_r0];
		vv_s_pos[i_r0][0] = 0.0;
		vv_r_pos[i_r0][0] = r0_pos[i_r0];

		vv_a_neg[i_r0][0] = a0_neg[i_r0];
		vv_v_neg[i_r0][0] = v0_neg[i_r0];
		vv_s_neg[i_r0][0] = 0.0;
		vv_r_neg[i_r0][0] = r0_neg[i_r0];

		vv_E1[i_r0][0] = 0.0;
		vv_E2[i_r0][0] = 0.0;
		vv_E [i_r0][0] = 0.0;
	}
}
#endif /*ALGORITHM_VERSION_2*/
double get_c()
{
	return g_c;
}

#ifdef ALGORITHM_VERSION_1
double interpolate(timevalue t_zap, double * values)
{
	double n_t 
		= 0 == v_n_t 
		? t_zap - g_t_start
		: v_n_t * (t_zap - g_t_start) / (v_t[v_n_t] - v_t[0]);

	if (n_t < (double)v_n_t)
	{
		// результат может быть взят с помощью линейной интерполяции ранее рассчитнных значений
		int n1 = (int)floor(n_t);
		int n2 = (int)ceil(n_t);

		while (t_zap < v_t[n1] && n1 > 0)
		{
			//printf("t_zap %0.25e < v_t[%d] %0.25e\n", t_zap, n1, v_t[n1]);
			n1 -= 1;
			n2 -= 1;
		}

		if (t_zap < v_t[n1])
		{
			printf("t_zap\n%0.25e < v_t[%d]\n%0.25e\n", t_zap, n1, v_t[n1]);
			int * p = 0;
			*p += 1;
		}

		while (t_zap > v_t[n2] && n2 < v_n_t - 1)
		{
			//printf("t_zap %0.25e > v_t[%d] %0.25e\n", t_zap, n2, v_t[n2]);
			n1 += 1;
			n2 += 1;
		}
	
		if (t_zap > v_t[n2])
		{
			printf("t_zap %0.25e > v_t[%d] %0.25e\n", t_zap, n2, v_t[n2]);
			int * p = 0;
			*p += 1;
		}

		double part = (t_zap - v_t[n1]) / (v_t[n2] - v_t[n1]);

		double value = values[n1] + part * (values[n2] - values[n1]);

		if (isnan(value))
		{
			printf("n1 %d n2 %d v_n_t %d t_zap %0.25e n_t %f\n", n1, n2, v_n_t, t_zap, n_t);
			int * p = 0;
			*p += 1;			
		}
		assert(!isnan(value));		
		return value;
	}

	if (n_t <= (double)v_n_t + epsilon_n && (double)v_n_t <= n_t)
	{
		int n = v_n_t;
		double value = values[n];
		assert(!isnan(value));
		return value;
	}


	printf("n_t = %0.20f v_n_t = %d\n", n_t, v_n_t);

	int * p = 0;
	*p+=1;

	assert(0);	
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
double get_a_ex1(timevalue t_zap, double q)
{
	if (t_zap < g_t_start)
		return 0;

	double * v_a = q > 0 ? v_a_pos : v_a_neg;
	double a = interpolate(t_zap, v_a);
	return a;
}

/* установить значение поля в точке наблюдения R0 в момент t */

void set_E_ex_1(int v_n_t, int v_n_r, double E)
{
	v_E[v_n_r][v_n_t] = E;
}

/* установить значение ускорения слоя исходя из его текущего радиуса и значения поля в текущий момент на этом радиусе*/
double set_a_ex1(timevalue t, double r, acceleration a0, timevalue t_a0, double q, double m, double * E)
{
	double * v_a = q > 0 ? v_a_pos : v_a_neg;

	assert(!isnan(*E));

	double a = multiplier_a * (*E) * q / m;
	if (t <= t_a0)
	{
		a += a0;
	}

	assert(!isnan(a));
	v_a[v_n_t + 1] = a;
	return a;
}

double get_v_ex1(timevalue t_zap, velocity v0, double q)
{
	assert(v0 < v_max);

	if (t_zap < g_t_start)
		return v0;

	double * v_v = q > 0 ? v_v_pos : v_v_neg;
	double v = interpolate(t_zap, v_v);
	return v;
}

double set_v_ex1(timevalue t, double v0, acceleration a0, timevalue t_a0, double q, double m)
{
	double * v_v = q > 0 ? v_v_pos : v_v_neg;
	double * v_a = q > 0 ? v_a_pos : v_a_neg;

	double t1 = v_t[v_n_t];
	double t2 = v_t[v_n_t + 1];
	double dt = t2 - t1;
	double a1 = v_a[v_n_t];
	double a2 = v_a[v_n_t + 1];
	double a = (a1 + a2) / 2;
	double v1 = v_v[v_n_t];
	// a = dv / dt
	// dv = a * dt
	double dv = a * dt;
	assert(!isnan(dv));
	double v2 = v1 + dv;
	if (fabs(v2) >= g_c)
	{
		printf("v1/c = %f, v2/c = %f, dv/c = %f\n", v1/g_c, v2/g_c, dv/g_c);
		int * p = 0;
		*p += 1;
	}
	v_v[v_n_t + 1] = v2;
	return v2;
}

double get_s_ex1(timevalue t_zap, double v0, double q)
{
	assert(v0 < v_max);
	double s;
	if (t_zap < g_t_start)
	{
		s = v0*(t_zap - g_t_start);
		DBG_INFO("get_s1 t_zap=%f returns %f\n", t_zap, s);
		return s;
	}

	double * v_s = q > 0 ? v_s_pos : v_s_neg;
	s = interpolate(t_zap, v_s);
	return s;
}

double set_s_ex1(timevalue t, double r0, double v0, double r_min, double q)
{
	double * v_s = q > 0 ? v_s_pos : v_s_neg;
	double * v_v = q > 0 ? v_v_pos : v_v_neg;
	double * v_a = q > 0 ? v_a_pos : v_a_neg;

	double t1 = v_t[v_n_t];
	double t2 = v_t[v_n_t + 1];
	double dt = t2 - t1;
	double a1 = v_a[v_n_t];
	double a2 = v_a[v_n_t + 1];
	double a = (a1 + a2) / 2;
	double v1 = v_v[v_n_t];
	double v2 = v_v[v_n_t + 1];
	double v = (v1 + v2) / 2;

	double s1 = v_s[v_n_t];

	// ds = v * dt + a * dt * dt / 2
	double ds = v * dt + a * dt * dt / 2;
	double s2 = s1 + ds;

	if (r0 + s2 < r_min)
	{
		printf("Warning: r %0.20f < r_min %0.20f\n", r0 + s2, r_min);
		//r0 + s = r_min; // ???????
		int* p = 0;
		*p += 1;
		assert(0);

	}

	v_s[v_n_t + 1] = s2;

	return s2;
}


/* расстояние от заряда до центра сферы в запаздывающий момент времени */
int get_r_ex1(double q, timevalue t_zap, double r0, double v0, double r_min, double * r)
{
	int error = 0;
	double s = get_s_ex1(t_zap, v0, q);
	*r = r0 + s;
	if (*r < r_min)
	{
		printf("Warning: r %0.20f < r_min %0.20f\n", *r, r_min);
		*r = r_min;
		error = 1;
		int* p = 0;
		*p += 1;
		assert(0);
	}
	assert(*r > 0.0);
	return error;
}
#endif /*ALGORITHM_VERSION_1*/

#ifdef ALGORITHM_VERSION_0
double get_a(timevalue t_zap, acceleration a0)
{
	double t_max;
#ifdef WITHOUT_ACCELERATION_BEFORE_TSTART
	if (t_zap < g_t_start)
		return 0;
#endif
	t_max = g_t_start + v_max / a0;
	if (t_zap >= t_max)
		return 0;
	return a0;
}

/* радиальная скорость заряда */

double get_v(timevalue t_zap, double v0, acceleration a0)
{
	assert(v0 < v_max);
	double t_max;
#ifdef WITHOUT_ACCELERATION_BEFORE_TSTART
	if (t_zap < g_t_start)
		return v0;
#endif
	t_max = g_t_start + v_max / a0;
	if (t_zap >= t_max)
		return v_max;
	return v0 + a0*(t_zap-g_t_start);
}
/* перемещение заряда */
double get_s(timevalue t_zap, double v0, acceleration a0)
{
	assert(v0 < v_max);
	double dt_start, dt_max;
	double t_max;
	double s;
#ifdef WITHOUT_ACCELERATION_BEFORE_TSTART
	if (t_zap < g_t_start)
	{
		s = v0*(t_zap - g_t_start);
		DBG_INFO("get_s1 t_zap=%f a0=%f returns %f\n", t_zap, a0, s);
		return s;
	}
#endif
	if (a0 > 0.0)
	{
		t_max = g_t_start + (v_max - v0) / a0;
		assert(t_max > g_t_start);
	}
	else if (a0 < 0.0)
	{
		t_max = g_t_start + (-v_max + v0) / a0;
		assert(t_max > g_t_start);
	}

	if (a0 != 0.0 && t_zap >= t_max)
	{
		DBG_INFO("get_s t_max %f\n", t_max);

		dt_start = (t_max - g_t_start);
		dt_max = (t_zap - t_max);
		s = v0 * dt_start + a0*dt_start*dt_start/2 + dt_max*v_max;
		DBG_INFO("get_s2 t_zap=%f a0=%f returns %f\n", t_zap, a0, s);
		return s;
	}

	dt_start = (t_zap - g_t_start);
	s = v0 * dt_start + a0*dt_start*dt_start/2;
	DBG_INFO("get_s t_zap=%f a0=%f dt_start=%f returns %f\n", t_zap, a0, dt_start, s);
	return s;
}

/* расстояние от заряда до центра сферы в запаздывающий момент времени */
int get_r(double q, timevalue t_zap, double r0, double v0, acceleration a0, double r_min, double * r)
{
	int error = 0;
	*r = r0 + get_s(t_zap, v0, a0);
#ifdef USE_MINIMAL_RADIUS
	if (*r < r_min)
	{
		DBG_INFO("Warning: r %f < r_min %f\n", *r, r_min);
		*r = r_min;
		error = 1;
	}
#endif
	assert(*r > 0.0);
	return error;
}
#endif /*ALGORITHM_VERSION_0*/
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

/*
tzap2(t, r__0, v__0, a__0, R__0, theta) = solve(c^2*(t-t__zap)^2 = R__0^2-2*R__0*r(t__zap, r__0, v__0, a__0)*cos(theta)+r(t__zap, r__0, v__0, a__0)^2, t__zap);

RootOf(a__0^2*_Z^4+4*v__0*_Z^3*a__0+(-4*R__0*a__0*cos(theta)+4*a__0*r__0+4*v__0^2-36)*_Z^2+(-8*R__0*v__0*cos(theta)+8*r__0*v__0+72*t)*_Z-8*R__0*cos(theta)*r__0+4*R__0^2+4*r__0^2-36*t^2);

allvalues(subs(t = 0, v__0 = 0, tzap2(t, r__0, v__0, a__0, R__0, theta)));

sqrt(2)*sqrt(R__0*a__0*cos(theta)-a__0*r__0+sqrt(cos(theta)^2*R__0^2*a__0^2-R__0^2*a__0^2+18*R__0*a__0*cos(theta)-18*a__0*r__0+81)+9)/a__0,
-sqrt(2)*sqrt(R__0*a__0*cos(theta)-a__0*r__0+sqrt(cos(theta)^2*R__0^2*a__0^2-R__0^2*a__0^2+18*R__0*a__0*cos(theta)-18*a__0*r__0+81)+9)/a__0,
sqrt(2*R__0*a__0*cos(theta)-2*a__0*r__0-2*sqrt(cos(theta)^2*R__0^2*a__0^2-R__0^2*a__0^2+18*R__0*a__0*cos(theta)-18*a__0*r__0+81)+18)/a__0,
-sqrt(2*R__0*a__0*cos(theta)-2*a__0*r__0-2*sqrt(cos(theta)^2*R__0^2*a__0^2-R__0^2*a__0^2+18*R__0*a__0*cos(theta)-18*a__0*r__0+81)+18)/a__0;
*/

/* численный расчёта запаздывающего момента */
int calc_tzap(double q, timevalue t, double R0, double r0, double v0, acceleration a0, double theta, double r_min, double * t2)
{
	int err, error = 0;
#ifdef CALC_LW_WITHOUT_LAGGING
	*t2 = g_t_start;
#else
	double epsilon = 1.0e-15;
	double t1;
	*t2 = t;
	double dt;
	double v1,v2;
	double r, R, R_pre = DBL_MAX;
	double dR, dR_pre = DBL_MAX;
	double R_tmp;

	int i = 0;
	double n = 0.9;
#if 0
	double v;
	v = get_v(t, v0, a0);      /* скорость заряда в текущий момент времени t */

	DBG_INFO("calc_tzap(t=%f, v = %f, R0=%f, r0=%f, v0=%f, a0=%f, theta=%f)\n", t, v, R0, r0, v0, a0, theta);
	assert(v < c);
#endif

	/*
	DBG_INFO("epsilon=%e\n", epsilon);
	DBG_INFO("t1=%f\n", t1);
	DBG_INFO("t2=%f\n", *t2);
	*/

	do
	{
		t1 = *t2;                 /* итерационный "текущий" момент времени - на первой итерации текущее время наблюдения */


#ifdef ALGORITHM_VERSION_0
		err = get_r(q, t1, r0, v0, a0, r_min, &r);   /* итерационная координата заряда */
#endif
#ifdef ALGORITHM_VERSION_1
		err = get_r_ex1(q, t1, r0, v0, r_min, &r);
#endif
#ifdef ALGORITHM_VERSION_2
		err = get_r_ex2(q, t1, r0, v0, r_min, &r);
#endif
		if (0 != err)
		{
			error += 1;
		}

		assert(r >= 0.0);
		R = get_R(R0, r, theta); /* итерационный радиус - на первой итерации текущий радиус                             */
		*t2 = t - R / g_c;          /* время прохождения сигнала от итерационной координаты в точку наблюдения             */
		                         /* итерационный "запаздывающий" момент времени t2                                      */
		dR = g_c*(t-t1) - R;       /**/
#ifdef ALGORITHM_VERSION_0
		v1 = get_v(t1, v0, a0);      /* скорость заряда в итерационный "текущий" момент времени t1                          */
		v2 = get_v(*t2, v0, a0);      /* скорость заряда в итерационный "запаздывающий" момент времени t2                    */
#endif
#ifdef ALGORITHM_VERSION_1
		v1 = get_v_ex1(t1, v0, q);
		v2 = get_v_ex1(*t2, v0, q);
#endif
#ifdef ALGORITHM_VERSION_2
		v1 = get_v_ex2(t1, v0, q);
		v2 = get_v_ex2(*t2, v0, q);
#endif
		DBG_INFO("t2=%f t1=%f t=%f v1 = %f, v2 = %f, v = %f, R=%f dR=%e dR_pre=%e ", *t2, t1, t, v1, v2, v, R, dR, dR_pre);
		if (v1 >= g_c || v2 >= g_c)
		{
			double v;
#ifdef ALGORITHM_VERSION_0			
			v = get_v(t, v0, a0);      /* скорость заряда в текущий момент времени t */
#endif
#ifdef ALGORITHM_VERSION_1
			v = get_v_ex1(t, v0, q);
#endif
			printf("t2=%f t1=%f t=%f v1 = %f, v2 = %f, v = %f, R=%f dR=%e dR_pre=%e\n", *t2, t1, t, v1, v2, v, R, dR, dR_pre);
			printf("v1/c = %f, v2/c = %f, v/c = %f\n", v1/g_c, v2/g_c, v/g_c);
			int * p = 0;
			*p += 1;
		}
		assert(v1 < g_c);
		assert(v2 < g_c);
#if 1
		if (i > 1 && fabs(dR) - fabs(dR_pre) < 1.0e-3)
		{
			int j = 0;

			do
			{
				R_tmp = R_pre + n * (R - R_pre);
				DBG_INFO("R_tmp = R_pre + n * (R - R_pre);= %f ", R_tmp);

				*t2 = t -  R_tmp / g_c;
				dR = g_c*(t-t1) - R_tmp;
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

