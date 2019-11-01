//#define CALC_LW_WITHOUT_LAGGING
//#define WITHOUT_ACCELERATION_BEFORE_TSTART
//#define SI
#ifdef SI
#define c 299792458.0
#else
#define c 3.0
#endif

typedef double time;
typedef double timespan;
typedef double acceleration;
typedef double velocity;

double get_c();

/* ускорение заряда */
double get_a(time t_zap, acceleration a0);

/* радиальная скорость заряда */
double get_v(time t_zap, velocity v0, acceleration a0);

/* перемещение заряда */
double get_s(time t_zap, velocity v0, acceleration a0);

/* расстояние от заряда до центра сферы в запаздывающий момент времени */
int get_r(double q, time t_zap, double r0, velocity v0, acceleration a0, double r_min, double * r);

/* расстояние от заряда до точки наблюдения в запаздывающий момент времени */
double get_R(double R0, double r, double theta);

/* численный расчёта запаздывающего момента */
int calc_tzap(double q, time t, double R0, double r0, velocity v0, acceleration a0, double theta, double r_min, double * t2);

