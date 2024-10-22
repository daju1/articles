#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_math.h>

struct pendulum
{
    double m;
    double g;
    double R;
    double y[2];
    double t;
    double xc, yc;
};

void init (struct pendulum * p, double m, double g, double R, double p0, double q0, double t0, double xc, double yc);
void init_left (double m, double g, double R, double p0, double q0, double t0, double xc, double yc);
void init_right (double m, double g, double R, double p0, double q0, double t0, double xc, double yc);

void alloc(struct pendulum * pendulum, gsl_odeiv2_system * sys, gsl_odeiv2_driver ** driver);
void alloc_left();
void alloc_right();

int apply(struct pendulum * pendulum, gsl_odeiv2_driver * driver,
          double ti, double * pt, double *pMomenta, double *pq,
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
typedef double force;


int find_period_by_newton_root(struct pendulum * pendulum, gsl_odeiv2_driver * driver,
                               timevalue * pT, double *pf);
void calc_pendulum_period(struct pendulum * pendulum, gsl_odeiv2_driver * driver,
                          double * p_init_T, double max_ti, double dti);


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
int klw(struct pendulum * pendulum, gsl_odeiv2_driver * driver,
        coordinate x, coordinate y, coordinate z, timevalue t, timevalue * pt2,
        double *pk, coordinate * rlagerror);
// Радиус Лиенара Вихерта
int Rlw(struct pendulum * pendulum, gsl_odeiv2_driver * driver,
        coordinate x, coordinate y, coordinate z, timevalue t,
        double *pRlw, coordinate * rlagerror);

// phi_lw - скалярный потенциал Лиенара Вихерта
int philw(struct pendulum * pendulum, gsl_odeiv2_driver * driver,
          coordinate x, coordinate y, coordinate z, timevalue t,
          charge q,
          double *pphi, coordinate * rlagerror);

// A_lw - векторный потенциал Лиенара Вихерта
int Alw(struct pendulum * pendulum, gsl_odeiv2_driver * driver,
        coordinate x, coordinate y, coordinate z, timevalue t,
        charge q,
        field * A_x, field * A_y, field * A_z, coordinate * rlagerror);


int electr_magnet(struct pendulum * pendulum, gsl_odeiv2_driver * driver,
                  coordinate x, coordinate y, coordinate z, timevalue t,
                  charge q,
                  field * E_x, field * E_y, field * E_z, field * B_x, field * B_y, field * B_z, coordinate * rlagerror);

int calc_fields(double k, distance r,
                distance nx,
                distance ny,
                distance nz,
                charge q,
                coordinate sx, coordinate sy, coordinate sz,
                velocity vx, velocity vy, velocity vz,
                acceleration wx, acceleration wy, acceleration wz,
                double dot_wx, double dot_wy, double dot_wz,
                field * E_x, field * E_y, field * E_z,
                field * B_x, field * B_y, field * B_z,
                field * A_x, field * A_y, field * A_z);

int electr_magnet_ex(struct pendulum * pendulum, gsl_odeiv2_driver * driver,
                   coordinate x, coordinate y, coordinate z, timevalue t,
                   charge q,
                   field * E_x, field * E_y, field * E_z,
                   field * B_x, field * B_y, field * B_z,
                   field * A_x, field * A_y, field * A_z,
                   coordinate * rlagerror);