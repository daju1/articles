
void init (double m, double g, double R, double p0, double q0, double t0, double xc, double yc);

void alloc();

int apply(double ti, double * pt, double *pMomenta, double *pq,
          double * pdot_q, double * pddot_q, double * pdddot_q,
          double * psx, double * psy,
          double * pvx, double * pvy,
          double * pwx, double * pwy,
          double * pdot_wx, double * pdot_wy
         );

typedef double timevalue;
typedef double timespan;
typedef double power;
typedef double acceleration;
typedef double velocity;
typedef double charge;
typedef double coordinate;
typedef double distance;
typedef double distance2;
typedef double mass;
typedef double potential;
typedef double field;
typedef double angle;
typedef double anglevelocity;


int find_period_by_newton_root(timevalue * pT, double *pf);
void calc_pendulum_period(double * p_init_T, double max_ti, double dti);


void cset_c(long double _c);
velocity cget_c();

// расчет итерациями запаздывающего момента

void cset_timespan_Epsilon(long double _eps);
void cset_distance_Epsilon(long double _eps);
void cset_max_steps(int _max_steps);
timespan cget_timespan_Epsilon();
void cset_min_newton_step(long double min_step);
void cset_newton_step_multiplier(long double multiplier);

// отношение радиуса Лиенара Вихерта к длине радиус-вектора
int klw(coordinate x, coordinate y, coordinate z, timevalue t, timevalue * pt2,
        double *pk, coordinate * rlagerror);
// Радиус Лиенара Вихерта
int Rlw(coordinate x, coordinate y, coordinate z, timevalue t,
        double *pRlw, coordinate * rlagerror);

// phi_lw - скалярный потенциал Лиенара Вихерта
int philw(coordinate x, coordinate y, coordinate z, timevalue t,
          charge q,
          double *pphi, coordinate * rlagerror);

// A_lw - векторный потенциал Лиенара Вихерта
int Alw(coordinate x, coordinate y, coordinate z, timevalue t,
        charge q,
        field * A_x, field * A_y, field * A_z, coordinate * rlagerror);


int electr_magnet(coordinate x, coordinate y, coordinate z, timevalue t,
                   charge q,
                   field * E_x, field * E_y, field * E_z, field * B_x, field * B_y, field * B_z, coordinate * rlagerror);


int electr_magnet_ex(coordinate x, coordinate y, coordinate z, timevalue t,
                   charge q,
                   field * E_x, field * E_y, field * E_z,
                   field * B_x, field * B_y, field * B_z,
                   field * A_x, field * A_y, field * A_z,
                   coordinate * rlagerror);