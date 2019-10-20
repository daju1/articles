//#define CALC_LW_WITHOUT_LAGGING

#define c 3.0

/* ускорение заряда */
double get_a(double t_zap, double a0);

/* радиальная скорость заряда */
double get_v(double t_zap, double v0, double a0);

/* перемещение заряда */
double get_s(double t_zap, double v0, double a0);

/* расстояние от заряда до центра сферы в запаздывающий момент времени */
int get_r(double t_zap, double r0, double v0, double a0, double r_min, double * r);

/* расстояние от заряда до точки наблюдения в запаздывающий момент времени */
double get_R(double R0, double r, double theta);

/* численный расчёта запаздывающего момента */
int calc_tzap(double t, double R0, double r0, double v0, double a0, double theta, double r_min, double * t2);

