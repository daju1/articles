//#define CALC_LW_WITHOUT_LAGGING
#define WITHOUT_ACCELERATION_BEFORE_TSTART
#define USE_MINIMAL_RADIUS

//#define ALGORITHM_VERSION_0
#define ALGORITHM_VERSION_1
//#define ALGORITHM_VERSION_2

#define SI

typedef double timevalue;
typedef double timespan;
typedef double acceleration;
typedef double velocity;

double get_c();
#ifdef ALGORITHM_VERSION_0
/* ускорение заряда */
double get_a(timevalue t_zap, acceleration a0);

/* радиальная скорость заряда */
double get_v(timevalue t_zap, velocity v0, acceleration a0);

/* перемещение заряда */
double get_s(timevalue t_zap, velocity v0, acceleration a0);

/* расстояние от заряда до центра сферы в запаздывающий момент времени */
int get_r(double q, timevalue t_zap, double r0, velocity v0, acceleration a0, double r_min, double * r);
#endif /* ALGORITHM_VERSION_0 */

double get_dt();
double get_dr();
int get_nt();
int get_nr();

#ifdef ALGORITHM_VERSION_1
void init_array_1(double a0_pos, velocity v0_pos, double r0_pos, double a0_neg, velocity v0_neg, double r0_neg);

double get_a_ex1(timevalue t_zap, double q);
double get_v_ex1(timevalue t_zap, velocity v0, double q);
double get_s_ex1(timevalue t_zap, double v0, double q);
/* расстояние от заряда до центра сферы в запаздывающий момент времени */
int get_r_ex1(double q, timevalue t_zap, double r0, double v0, double r_min, double * r);

void set_E_ex1(timevalue t, double R0, double E);
void set_E_ex_1(int v_n_t, int v_n_r, double E);

double set_a_ex1(timevalue t, double r, acceleration a0, timevalue t_a0, double q, double m, double * E);
double set_v_ex1(timevalue t, double v0, acceleration a0, timevalue t_a0, double q, double m);
double set_s_ex1(timevalue t, double r0, double v0, double r_min, double q);
#endif /*ALGORITHM_VERSION_1*/

/* расстояние от заряда до точки наблюдения в запаздывающий момент времени */
double get_R(double R0, double r, double theta);

/* численный расчёта запаздывающего момента */
int calc_tzap(double q, timevalue t, double R0, double r0, velocity v0, acceleration a0, double theta, double r_min, double * t2);

