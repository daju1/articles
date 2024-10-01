
void cset_c(long double _c);

velocity cget_c();


// вращение заряда в плоскости xy вокруг центра в точке
// coordinate xc, coordinate yc, coordinate zc
// по радиусу distance R
// с угловой частотой anglevelocity omega
// начальный угол в момент времени ноль angle alpha

typedef coordinate (*Coordinate)(timevalue t_zap,
                                 coordinate xc, coordinate yc, coordinate zc,
                                 distance R, anglevelocity omega, angle alpha);
typedef velocity (*Velocity)(timevalue t_zap,
                                 coordinate xc, coordinate yc, coordinate zc,
                                 distance R, anglevelocity omega, angle alpha);
typedef acceleration (*Acceleration)(timevalue t_zap,
                                 coordinate xc, coordinate yc, coordinate zc,
                                 distance R, anglevelocity omega, angle alpha);

// расчет итерациями запаздывающего момента

void cset_timespan_Epsilon(long double _eps);
void cset_distance_Epsilon(long double _eps);
void cset_max_steps(int _max_steps);

void cset_no_retardation_test(int test);

timespan cget_timespan_Epsilon();

void cset_min_newton_step(long double min_step);
void cset_newton_step_multiplier(long double multiplier);


typedef timevalue (*Tlag)(coordinate x, coordinate y, coordinate z, timevalue t,
                          Coordinate sx, Coordinate sy, Coordinate sz,
                          Velocity vx, Velocity vy, Velocity vz,
                          coordinate xc, coordinate yc, coordinate zc,
                          distance R, anglevelocity omega, angle alpha);


/*
timevalue newton_root_func(coordinate x, coordinate y, coordinate z,
                           timevalue t, timevalue t2,
                           Coordinate sx, Coordinate sy, Coordinate sz,
                           coordinate xc, coordinate yc, coordinate zc,
                           distance R, anglevelocity omega, angle alpha);

timevalue newton_root_derivative(coordinate x, coordinate y, coordinate z,
                                 timevalue t, timevalue t2,
                                 Coordinate sx, Coordinate sy, Coordinate sz,
                                 Velocity vx, Velocity vy, Velocity vz,
                                 coordinate xc, coordinate yc, coordinate zc,
                                 distance R, anglevelocity omega, angle alpha);
*/

void newton_root_derivative(coordinate x, coordinate y, coordinate z,
                            timevalue t, timevalue t2,
                            Coordinate sx, Coordinate sy, Coordinate sz,
                            Velocity vx, Velocity vy, Velocity vz,
                            coordinate xc, coordinate yc, coordinate zc,
                            distance R, anglevelocity omega, angle alpha,
                            long double * cdt, long double * d, long double * f1,
                            long double * dfdt);

int NewtonIt(long double step,
             coordinate x, coordinate y, coordinate z,
             timevalue t, timevalue t2,
             Coordinate sx, Coordinate sy, Coordinate sz,
             Velocity vx, Velocity vy, Velocity vz,
             coordinate xc, coordinate yc, coordinate zc,
             distance R, anglevelocity omega, angle alpha,
             timevalue * res, long double * f);

int find_newton_root(coordinate x, coordinate y, coordinate z, timevalue t, timevalue * pt2,
                     Coordinate sx, Coordinate sy, Coordinate sz,
                     Velocity vx, Velocity vy, Velocity vz,
                     coordinate xc, coordinate yc, coordinate zc,
                     distance R, anglevelocity omega, angle alpha);

int tlag(coordinate x, coordinate y, coordinate z, timevalue t,
         Coordinate sx, Coordinate sy, Coordinate sz,
         Velocity vx, Velocity vy, Velocity vz,
         coordinate xc, coordinate yc, coordinate zc,
         distance R, anglevelocity omega, angle alpha, timevalue * pt2, coordinate * rlagerror);

timevalue tlag_test(coordinate x, coordinate y, coordinate z, timevalue t1, timevalue t2,
                    Coordinate sx, Coordinate sy, Coordinate sz,
                    coordinate xc, coordinate yc, coordinate zc,
                    distance R, anglevelocity omega, angle alpha, long double * d, long double * cdt);

/*
Здесь нужно отметить во избежание путаницы, что радиус Лиенара Вихерта $k$ в формуле для электрического и магнитного полей и $k$ в программе отличаются друг от друга тем, что в программе $k$ нормирован на единицу, а в формуле нет. При переходе от формул к програмным кодам по сути произведено преобразование $k \rightarrow k\cdot r$
*/

void calc_k(coordinate x, coordinate y, coordinate z, timevalue t,
            Coordinate sx, Coordinate sy, Coordinate sz,
            Velocity vx, Velocity vy, Velocity vz,
            timevalue t2,
            long double * k, distance * r, distance * nx, distance * ny, distance * nz,
            coordinate xc, coordinate yc, coordinate zc,
            distance R, anglevelocity omega, angle alpha);

// отношение радиуса Лиенара Вихерта к радиусу
int klw(coordinate x, coordinate y, coordinate z, timevalue t,
                Coordinate sx, Coordinate sy, Coordinate sz,
                Velocity vx, Velocity vy, Velocity vz,
                coordinate xc, coordinate yc, coordinate zc,
                distance R, anglevelocity omega, angle alpha, long double *, coordinate * rlagerror);

// Радиус Лиенара Вихерта
int Rlw(coordinate x, coordinate y, coordinate z, timevalue t,
                Coordinate sx, Coordinate sy, Coordinate sz,
                Velocity vx, Velocity vy, Velocity vz,
                coordinate xc, coordinate yc, coordinate zc,
                distance R, anglevelocity omega, angle alpha, long double *, coordinate * rlagerror);

// phi_lw - скалярный потенциал Лиенара Вихерта
int philw(coordinate x, coordinate y, coordinate z, timevalue t,
                  Coordinate sx, Coordinate sy, Coordinate sz,
                  Velocity vx, Velocity vy, Velocity vz,
                  charge q,
                  coordinate xc, coordinate yc, coordinate zc,
                  distance R, anglevelocity omega, angle alpha, long double *, coordinate * rlagerror);

// A_lw - векторный потенциал Лиенара Вихерта
int Alw(coordinate x, coordinate y, coordinate z, timevalue t,
         Coordinate sx, Coordinate sy, Coordinate sz,
         Velocity vx, Velocity vy, Velocity vz,
         charge q,
         field * A_x, field * A_y, field * A_z, coordinate * rlagerror,
         coordinate xc, coordinate yc, coordinate zc,
         distance R, anglevelocity omega, angle alpha
       );


int electr_magnet(coordinate x, coordinate y, coordinate z, timevalue t,
                   Coordinate sx, Coordinate sy, Coordinate sz,
                   Velocity vx, Velocity vy, Velocity vz,
                   Acceleration wx, Acceleration wy, Acceleration wz,
                   charge q,
                   field * E_x, field * E_y, field * E_z,
                   field * B_x, field * B_y, field * B_z,
                   coordinate * rlagerror,
                   coordinate xc, coordinate yc, coordinate zc,
                   distance R, anglevelocity omega, angle alpha);

