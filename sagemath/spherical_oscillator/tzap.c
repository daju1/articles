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
long double multiplier_E = LIGHT_VELONCITY * LIGHT_VELONCITY / 10000000.0;
#else
long double multiplier_E = 1.0;
#endif

velocity g_c = LIGHT_VELONCITY;
static velocity v_max = 0.999 * LIGHT_VELONCITY;

static timevalue g_t_start = T_START;    // момент включения ускорения
//static double g_dt = DT;              // шаг времени
static timevalue g_t_finish = T_FINISH;  // момент времени окончания расчёта
static distance g_dr = DR;              // шаг координаты

// одномерные массивы для сохранения истории при интегрировании только лишь по времени
// без учёта эволючии объёмного распределения заряда - упрощённый случай сферического конденсатора
// сохраняется только лишь история положительной обкладки и отрицательной обкладки
static coordinate g_r_start = R_START;
static coordinate g_r_finish = R_FINISH;
//static int v_Nt = (int)((T_FINISH - T_START) / DT);
static int v_Nt = 10000;
static int v_Nr = (int)((R_FINISH - R_START) / DR);

distance get_dr()
{
	return g_dr;
}

void set_dr(distance dr)
{
	g_dr = dr;
	v_Nr = (int)((g_r_finish - g_r_start) / g_dr);
}

void set_r_finish(coordinate r_finish)
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

static timevalue epsilon_n = 1e-8;
static distance epsilon_r = 1e-8;

timevalue * v_t;
int v_n_t = 0; // итератор полноты заполения двумерных массивов по координате времени

#ifdef ALGORITHM_VERSION_1
static acceleration * v_a_pos;
static velocity * v_v_pos;
static distance * v_s_pos;
static coordinate * v_r_pos;

static acceleration * v_a_neg;
static velocity * v_v_neg;
static distance * v_s_neg;
static coordinate * v_r_neg;

static field ** v_E1;
static field ** v_E2;
static field ** v_E;
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
void init_array_1(acceleration a0_pos, velocity v0_pos, coordinate r0_pos, acceleration a0_neg, velocity v0_neg, coordinate r0_neg)
{
	v_t = malloc(v_Nt * sizeof(timevalue));

	v_a_pos = malloc(v_Nt * sizeof(acceleration));
	v_v_pos = malloc(v_Nt * sizeof(velocity));
	v_s_pos = malloc(v_Nt * sizeof(distance));
	v_r_pos = malloc(v_Nt * sizeof(coordinate));

	v_a_neg = malloc(v_Nt * sizeof(acceleration));
	v_v_neg = malloc(v_Nt * sizeof(velocity));
	v_s_neg = malloc(v_Nt * sizeof(distance));
	v_r_neg = malloc(v_Nt * sizeof(coordinate));

	v_E1 = malloc(v_Nr * sizeof(field *));
	v_E2 = malloc(v_Nr * sizeof(field *));
	v_E  = malloc(v_Nr * sizeof(field *));

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
		v_E1[i_r] = malloc(v_Nt * sizeof(field));
		v_E2[i_r] = malloc(v_Nt * sizeof(field));
		v_E [i_r] = malloc(v_Nt * sizeof(field));

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
	v_t = malloc(v_Nt * sizeof(timevalue *));

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
long double interpolate(timevalue t_zap, long double * values)
{
	long double n_t
		= 0 == v_n_t 
		? t_zap - g_t_start
		: v_n_t * (t_zap - g_t_start) / (v_t[v_n_t] - v_t[0]);

	if (n_t < (long double)v_n_t)
	{
		// результат может быть взят с помощью линейной интерполяции ранее рассчитнных значений
		long n1 = (long)floor(n_t);
		long n2 = (long)ceil(n_t);
		if (n1 == n2)
		{
			int n = n1;
			if (fabs(t_zap - v_t[n]) < epsilon_n )
			{
				long double value = values[n];
				assert(!isnan(value));
				return value;
			}
			else
			{
				assert(0);
			}
		}

		while (t_zap < v_t[n1] && n1 > 0)
		{
			//printf("t_zap %0.25e < v_t[%d] %0.25e\n", t_zap, n1, v_t[n1]);
			n1 -= 1;
			n2 -= 1;
		}

		if (t_zap < v_t[n1])
		{
			printf("t_zap\n%0.25Le < v_t[%ld]\n%0.25Le\n", t_zap, n1, v_t[n1]);
			int * p = 0;
			*p += 1;
		}

		while (t_zap > v_t[n2] && n2 < v_n_t)
		{
			//printf("t_zap %0.25e > v_t[%d] %0.25e\n", t_zap, n2, v_t[n2]);
			n1 += 1;
			n2 += 1;
		}
	
		if (t_zap > v_t[n2])
		{
			printf("t_zap %0.25Le > v_t[%ld] %0.25Le\n", t_zap, n2, v_t[n2]);
			int * p = 0;
			*p += 1;
		}

		timevalue t1 = v_t[n1];
		timevalue t2 = v_t[n2];
		timespan dt = t2 - t1;

		long double part = (t_zap - v_t[n1]) / dt;

		long double value = values[n1] + part * (values[n2] - values[n1]);

		if (isnan(value))
		{
			//printf("n1 %ld n2 %ld v_n_t %d t_zap %0.25Le n_t %Lf\n", n1, n2, v_n_t, t_zap, n_t);
			if (fabs(t_zap - t1) < fabs(t_zap - t2))
			{
				value = value = values[n1];
			}
			else
			{
				value = value = values[n2];
			}
		}
		assert(!isnan(value));
		return value;
	}

	int n = v_n_t;
	if (fabs(t_zap - v_t[n]) < epsilon_n )
	{
		long double value = values[n];
		assert(!isnan(value));
		return value;
	}

	printf("n_t = %0.20Lf v_n_t = %d\n", n_t, v_n_t);

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
acceleration get_a_ex1(timevalue t_zap, charge q)
{
	if (t_zap < g_t_start)
		return 0;

	acceleration * v_a = q > 0 ? v_a_pos : v_a_neg;
	acceleration a = interpolate(t_zap, v_a);
	return a;
}

/* установить значение поля в точке наблюдения R0 в момент t */

void set_E_ex_1(int v_n_t, int v_n_r, field E)
{
	v_E[v_n_r][v_n_t] = E;
}

/* установить значение ускорения слоя исходя из его текущего радиуса и значения поля в текущий момент на этом радиусе*/
int set_a_ex1(timevalue t, timespan dt, coordinate r, power pw, timevalue t_a0, charge q, mass m, field E, field E1, field E2)
{
	int error = 0;
	acceleration * v_a = q > 0 ? v_a_pos : v_a_neg;
	velocity * v_v = q > 0 ? v_v_pos : v_v_neg;
	velocity v = interpolate(t - dt, v_v);

	assert(!isnan(E));

	// power pw
	// work dA = pw * dt
	// on the other hand dA = F * ds
	// where
	// ds = v * dt + a * dt^2 / 2
	// and
	// F = m*a - q * E
	// pw * dt = (m*a - q * E) * (v * dt + a * dt^2 / 2)
	// pw = (m*a - q * E) * (v + a * dt / 2)
	// solve(pw = (-E*q+a*m)*(v+(1/2)*a*dt), a)

	// E*dt*q - 2*m*v0 +- sqrt(E^2*dt^2*q^2 + 4*E*dt*m*q*v + 4*m^2*v0^2 + 8*dt*m*pw)
	// -----------------------------------------------------------------------------
	//                   2*dt*m

	//  E*q     v            / E^2*q^2    E*q*v      v^2      2 * pw \
	// ----- - ---- +- sqrt | -------- + -------- + ------ + -------- |
	//  2*m     dt           \ 4*m^2       dt*m      dt^2      dt*m  /

	power p = 0.0;
	if (t <= t_a0)
	{
		p = pw;
	}


	long double radical
		= E * E * q * q / (4 * m * m)
		+ E * q * v / (dt * m)
		+ v * v / (dt * dt)
		+ 2 * p / (dt * m);

	//acceleration a = E * q / m;
	//acceleration a1 = E1 * q / m;
	//acceleration a2 = E2 * q / m;

	if (radical < 0.0)
	{
		printf("radical %Le < 0.0\n", radical);
		return 1;
	}

	acceleration a = E * q / (2 * m) - v / dt + sqrt(radical);


	assert(!isnan(a));
	v_a[v_n_t + 1] = a;
	return error;
}

velocity get_v_ex1(timevalue t_zap, velocity v0, charge q)
{
	assert(v0 < v_max);

	if (t_zap < g_t_start)
		return v0;

	velocity * v_v = q > 0 ? v_v_pos : v_v_neg;
	velocity v = interpolate(t_zap, v_v);
	return v;
}

int set_v_ex1(timevalue t, velocity v0, timevalue t_a0, charge q, mass m, velocity * v2)
{
	int error = 0;
	velocity * v_v = q > 0 ? v_v_pos : v_v_neg;
	acceleration * v_a = q > 0 ? v_a_pos : v_a_neg;

	timevalue t1 = v_t[v_n_t];
	timevalue t2 = v_t[v_n_t + 1];
	timespan dt = t2 - t1;
	acceleration a1 = v_a[v_n_t];
	acceleration a2 = v_a[v_n_t + 1];
	acceleration a = (a1 + a2) / 2;
	velocity v1 = v_v[v_n_t];
	// a = dv / dt
	// dv = a * dt
	velocity dv = a * dt;
	assert(!isnan(dv));
	*v2 = v1 + dv;
	if (fabs(*v2) >= g_c)
	{
		printf("v1/c = %Lf, v2/c = %Lf, dv/c = %Lf\n", v1/g_c, (*v2)/g_c, dv/g_c);
		error = 1;
	}
	v_v[v_n_t + 1] = *v2;
	return error;
}

distance get_s_ex1(timevalue t_zap, velocity v0, charge q)
{
	assert(v0 < v_max);
	distance s;
	if (t_zap < g_t_start)
	{
		s = v0*(t_zap - g_t_start);
		DBG_INFO("get_s1 t_zap=%Lf returns %Lf\n", t_zap, s);
		return s;
	}

	distance * v_s = q > 0 ? v_s_pos : v_s_neg;
	s = interpolate(t_zap, v_s);
	return s;
}

distance set_s_ex1(timevalue t, coordinate r0, velocity v0, coordinate r_min, charge q)
{
	distance * v_s = q > 0 ? v_s_pos : v_s_neg;
	velocity * v_v = q > 0 ? v_v_pos : v_v_neg;
	acceleration * v_a = q > 0 ? v_a_pos : v_a_neg;


	timevalue t1 = v_t[v_n_t];
	timevalue t2 = v_t[v_n_t + 1];
	timespan dt = t2 - t1;
	acceleration a1 = v_a[v_n_t];
	acceleration a2 = v_a[v_n_t + 1];
	acceleration a = (a1 + a2) / 2;
	velocity v1 = v_v[v_n_t];
	velocity v2 = v_v[v_n_t + 1];
	velocity v = (v1 + v2) / 2;

	distance s1 = v_s[v_n_t];

	// ds = v * dt + a * dt * dt / 2
	distance ds = v * dt + a * dt * dt / 2;
	distance s2 = s1 + ds;

	if (r0 + s2 < r_min)
	{
		printf("Warning: r %0.20Lf < r_min %0.20Lf\n", r0 + s2, r_min);
		//r0 + s = r_min; // ???????
		int* p = 0;
		*p += 1;
		assert(0);
	}

	v_s[v_n_t + 1] = s2;

	return s2;
}


/* расстояние от заряда до центра сферы в запаздывающий момент времени */
int get_r_ex1(charge q, timevalue t_zap, coordinate r0, velocity v0, coordinate r_min, coordinate * r, int log)
{
	int error = 0;
	distance s = get_s_ex1(t_zap, v0, q);
	*r = r0 + s;
	if (log)
		printf(
			"t_zap = %Lf "
			"q = % 0.20Le "
			"r0 = % 0.20Le "
			"s = % 0.20Le "
			"*r = % 0.20Le\n"
			, t_zap
			, q
			, r0
			, s
			, *r
			);
	if (*r < r_min)
	{
		printf("Warning: r %0.20Lf < r_min %0.20Lf\n", *r, r_min);
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
		DBG_INFO("get_s1 t_zap=%Lf a0=%Lf returns %Lf\n", t_zap, a0, s);
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
		DBG_INFO("get_s t_max %Lf\n", t_max);

		dt_start = (t_max - g_t_start);
		dt_max = (t_zap - t_max);
		s = v0 * dt_start + a0*dt_start*dt_start/2 + dt_max*v_max;
		DBG_INFO("get_s2 t_zap=%Lf a0=%Lf returns %Lf\n", t_zap, a0, s);
		return s;
	}

	dt_start = (t_zap - g_t_start);
	s = v0 * dt_start + a0*dt_start*dt_start/2;
	DBG_INFO("get_s t_zap=%Lf a0=%Lf dt_start=%Lf returns %Lf\n", t_zap, a0, dt_start, s);
	return s;
}

/* расстояние от заряда до центра сферы в запаздывающий момент времени */
int get_r(charge q, timevalue t_zap, double r0, double v0, acceleration a0, double r_min, double * r)
{
	int error = 0;
	*r = r0 + get_s(t_zap, v0, a0);
#ifdef USE_MINIMAL_RADIUS
	if (*r < r_min)
	{
		DBG_INFO("Warning: r %Lf < r_min %Lf\n", *r, r_min);
		*r = r_min;
		error = 1;
	}
#endif
	assert(*r > 0.0);
	return error;
}
#endif /*ALGORITHM_VERSION_0*/
/* расстояние от заряда до точки наблюдения в запаздывающий момент времени */
distance get_R(coordinate R0, coordinate r, angle theta)
{
	//distance2 RR = R0*R0 - 2 * R0*r*cos(theta) + r*r;
	//distance2 RR2 = R0*R0*cos(theta) - 2 * R0*r*cos(theta) + r*r*cos(theta) + R0*R0*(1 - cos(theta)) + r*r*(1 - cos(theta));
	distance2 RR3 = (R0 - r)*(R0 - r)*cos(theta) + (R0*R0 + r*r)*(1 - cos(theta));
	//DBG_INFO("RR=%Lf, RR2 = %Lf, RR3 = %Lf, R0=%Lf, r=%Lf, theta=%Lf)\n", RR, RR2, RR3, R0, r, theta);
	assert(RR3 >= 0.0);
	distance R = sqrt(RR3);
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
int calc_tzap(charge q, timevalue t, coordinate R0, coordinate r0, velocity v0, power pw, angle theta, coordinate r_min, timevalue * t2)
{
	int err, error = 0;
#ifdef CALC_LW_WITHOUT_LAGGING
	*t2 = g_t_start;
#else

#ifdef SI
	long double epsilon = 1.0e-16;
	long double epsilon_dr = 1.0e-16;
#else
	long double epsilon = 1.0e-15;
	long double epsilon_dr = 1.0e-3;
#endif
	timevalue t1;
	*t2 = t;
	timespan dt;
	velocity v1,v2;
	coordinate r, R, R_pre = LDBL_MAX;
	distance dR, dR_pre = LDBL_MAX;
	coordinate R_tmp;

	int i = 0;
#if 0
	double v;
	v = get_v(t, v0, a0);      /* скорость заряда в текущий момент времени t */

	DBG_INFO("calc_tzap(t=%Lf, v = %Lf, R0=%Lf, r0=%Lf, v0=%Lf, a0=%Lf, theta=%Lf)\n", t, v, R0, r0, v0, a0, theta);
	assert(v < c);
#endif

	/*
	DBG_INFO("epsilon=%Le\n", epsilon);
	DBG_INFO("t1=%Lf\n", t1);
	DBG_INFO("t2=%Lf\n", *t2);
	*/

	do
	{
		t1 = *t2;                 /* итерационный "текущий" момент времени - на первой итерации текущее время наблюдения */


#ifdef ALGORITHM_VERSION_0
		err = get_r(q, t1, r0, v0, a0, r_min, &r);   /* итерационная координата заряда */
#endif
#ifdef ALGORITHM_VERSION_1
		err = get_r_ex1(q, t1, r0, v0, r_min, &r, 0);
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
		DBG_INFO("t2=%Lf t1=%Lf t=%Lf v1 = %Lf, v2 = %Lf, v = %Lf, R=%Lf dR=%Le dR_pre=%Le ", *t2, t1, t, v1, v2, v, R, dR, dR_pre);
		if (v1 >= g_c || v2 >= g_c)
		{
			velocity v;
#ifdef ALGORITHM_VERSION_0			
			v = get_v(t, v0, a0);      /* скорость заряда в текущий момент времени t */
#endif
#ifdef ALGORITHM_VERSION_1
			v = get_v_ex1(t, v0, q);
#endif
			printf("t2=%Lf t1=%Lf t=%Lf v1 = %Lf, v2 = %Lf, v = %Lf, R=%Lf dR=%Le dR_pre=%Le\n", *t2, t1, t, v1, v2, v, R, dR, dR_pre);
			printf("v1/c = %Lf, v2/c = %Lf, v/c = %Lf\n", v1/g_c, v2/g_c, v/g_c);
			int * p = 0;
			*p += 1;
		}
		assert(v1 < g_c);
		assert(v2 < g_c);
#if 1
		long double n = 0.9;
		if (i > 1 && fabs(dR) - fabs(dR_pre) < epsilon_dr)
		{
			int j = 0;

			do
			{
				R_tmp = R_pre + n * (R - R_pre);
				DBG_INFO("R_tmp = R_pre + n * (R - R_pre);= %Lf ", R_tmp);

				*t2 = t -  R_tmp / g_c;
				dR = g_c*(t-t1) - R_tmp;
				DBG_INFO("t2 = %Lf ", *t2);

				n *= 0.9;
				DBG_INFO("n = %Lf ", n);
				++j;
			}
			while (fabs(dR) > fabs(dR_pre));

			R = R_tmp;
		}
#endif
		dt = t1 - *t2;
		DBG_INFO("dt=%Le \n", dt);
		dR_pre = dR;
		R_pre = R;

		DBG_INFO("\n");
		++i;
	}
	while (fabs(dt) > epsilon);

	DBG_INFO("fabs(t1 - t2) = %Le fabs(t - t2) = %Le calc_tzap() result=%Lf\n", fabs(t1 - *t2), fabs(t - *t2), *t2);
#endif
	return error;
}

