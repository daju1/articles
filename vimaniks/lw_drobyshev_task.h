#include "lw.h"

typedef coordinate (*Coordinate)(timevalue t_zap);
typedef velocity (*Velocity)(timevalue t_zap);
typedef acceleration (*Acceleration)(timevalue t_zap);
typedef dotacceleration (*DotAcceleration)(timevalue t_zap);

int tlag(coordinate x, coordinate y, coordinate z, timevalue t,
         Coordinate sx, Coordinate sy, Coordinate sz,
         Velocity vx, Velocity vy, Velocity vz,
         timevalue * pt2, coordinate * rlagerror);
/*
Здесь нужно отметить во избежание путаницы, что радиус Лиенара Вихерта $k$ в формуле для электрического и магнитного полей и $k$ в программе отличаются друг от друга тем, что в программе $k$ нормирован на единицу, а в формуле нет. При переходе от формул к програмным кодам по сути произведено преобразование $k \rightarrow k\cdot r$
*/

void calc_k(coordinate x, coordinate y, coordinate z, timevalue t,
            Coordinate sx, Coordinate sy, Coordinate sz,
            Velocity vx, Velocity vy, Velocity vz,
            timevalue t2,
            long double * k, distance * r, distance * nx, distance * ny, distance * nz);

// отношение радиуса Лиенара Вихерта к радиусу
int klw(coordinate x, coordinate y, coordinate z, timevalue t,
        Coordinate sx, Coordinate sy, Coordinate sz,
        Velocity vx, Velocity vy, Velocity vz,
        long double *pk, coordinate * rlagerror);

// Радиус Лиенара Вихерта
int Rlw(coordinate x, coordinate y, coordinate z, timevalue t,
        Coordinate sx, Coordinate sy, Coordinate sz,
        Velocity vx, Velocity vy, Velocity vz,
        long double *pRlw, coordinate * rlagerror);

// phi_lw - скалярный потенциал Лиенара Вихерта
int philw(coordinate x, coordinate y, coordinate z, timevalue t,
          Coordinate sx, Coordinate sy, Coordinate sz,
          Velocity vx, Velocity vy, Velocity vz,
          charge q,
          long double *pphi, coordinate * rlagerror);

// A_lw - векторный потенциал Лиенара Вихерта
int Alw(coordinate x, coordinate y, coordinate z, timevalue t,
        Coordinate sx, Coordinate sy, Coordinate sz,
        Velocity vx, Velocity vy, Velocity vz,
        charge q,
        field * A_x, field * A_y, field * A_z,
        coordinate * rlagerror
       );

int electr_magnet(coordinate x, coordinate y, coordinate z, timevalue t,
                  Coordinate sx, Coordinate sy, Coordinate sz,
                  Velocity vx, Velocity vy, Velocity vz,
                  Acceleration wx, Acceleration wy, Acceleration wz,
                  charge q,
                  field * E_x, field * E_y, field * E_z,
                  field * B_x, field * B_y, field * B_z,
                  coordinate * rlagerror);

int electr_magnet_ex(coordinate x, coordinate y, coordinate z, timevalue t,
                     Coordinate sx, Coordinate sy, Coordinate sz,
                     Velocity vx, Velocity vy, Velocity vz,
                     Acceleration wx, Acceleration wy, Acceleration wz,
                     DotAcceleration dot_wx, DotAcceleration dot_wy, DotAcceleration dot_wz,
                     charge q,
                     field * E1_x, field * E1_y, field * E1_z,
                     field * E2_x, field * E2_y, field * E2_z,
                     field * E_x, field * E_y, field * E_z,
                     field * B_x, field * B_y, field * B_z,
                     field * A_x, field * A_y, field * A_z,
                     field * j_x, field * j_y, field * j_z,
                     long double * ra_c2,
                     long double * four_a_four_R_c2,
                     coordinate * rlagerror);