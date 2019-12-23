//#define CALC_LW_WITHOUT_LAGGING
#define WITHOUT_ACCELERATION_BEFORE_TSTART
#define USE_MINIMAL_RADIUS

#define ALGORITHM_VERSION_0
//#define ALGORITHM_VERSION_1
//#define ALGORITHM_VERSION_2

#define SI
#ifdef SI
//#define USE_NORM
#define LIGHT_VELONCITY 299792458.0
#else
#define LIGHT_VELONCITY 3.0
#endif

typedef long double timevalue;
typedef long double timespan;
typedef long double power;
typedef long double acceleration;
typedef long double velocity;
typedef long double charge;
typedef long double coordinate;
typedef long double distance;
typedef long double distance2;
typedef long double mass;
typedef long double potential;
typedef long double field;
typedef long double angle;


velocity get_c();
#ifdef ALGORITHM_VERSION_0

void init_array_0();

/* ускорение заряда */
acceleration get_a(timevalue t_zap, acceleration a0);

/* радиальная скорость заряда */
velocity get_v(timevalue t_zap, velocity v0, acceleration a0);

/* перемещение заряда */
distance get_s(timevalue t_zap, velocity v0, acceleration a0);

/* расстояние от заряда до центра сферы в запаздывающий момент времени */
int get_r(charge q, timevalue t_zap, coordinate r0, velocity v0, acceleration a0, coordinate r_min, coordinate * r);
#endif /* ALGORITHM_VERSION_0 */

int get_r_common(charge q, timevalue t, coordinate r0, velocity v0, acceleration a0, coordinate r_min, coordinate * r);

distance get_dr();
int get_nt();
int get_nr();

#ifdef ALGORITHM_VERSION_1
void init_array_1(acceleration a0_pos, velocity v0_pos, coordinate r0_pos, acceleration a0_neg, velocity v0_neg, coordinate r0_neg);

acceleration get_a_ex1(timevalue t_zap, charge q);
velocity get_v_ex1(timevalue t_zap, velocity v0, charge q);
distance get_s_ex1(timevalue t_zap, velocity v0, charge q);
/* расстояние от заряда до центра сферы в запаздывающий момент времени */
int get_r_ex1(charge q, timevalue t_zap, coordinate r0, velocity v0, coordinate r_min, coordinate * r, int log);

void set_E_ex_1(int v_n_t, int v_n_r, field E);

int set_a_ex1(timevalue t, timespan dt, coordinate r, power pw, timevalue t_a0, charge q, mass m, field E, field E1, field E2);
int set_v_ex1(timevalue t, velocity v0, timevalue t_a0, charge q, mass m, velocity * v2);

distance set_s_ex1(timevalue t, coordinate r0, velocity v0, coordinate r_min, charge q);
#endif /*ALGORITHM_VERSION_1*/

/* расстояние от заряда до точки наблюдения в запаздывающий момент времени */
distance get_R(coordinate R0, coordinate r, angle theta);

/* численный расчёта запаздывающего момента */
int calc_tzap(charge q, timevalue t, coordinate R0, coordinate r0, velocity v0, acceleration a0, angle theta, coordinate r_min, timevalue * t2);


void set_dr(distance dr);
void set_r_finish(coordinate r_finish);