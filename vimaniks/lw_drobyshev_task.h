#include "lw.h"

typedef coordinate (*Coordinate)(timevalue t_zap);
typedef velocity (*Velocity)(timevalue t_zap);
typedef acceleration (*Acceleration)(timevalue t_zap);

timevalue tlag(coordinate x, coordinate y, coordinate z, timevalue t,
               Coordinate sx, Coordinate sy, Coordinate sz);
/*
Здесь нужно отметить во избежание путаницы, что радиус Лиенара Вихерта $k$ в формуле для электрического и магнитного полей и $k$ в программе отличаются друг от друга тем, что в программе $k$ нормирован на единицу, а в формуле нет. При переходе от формул к програмным кодам по сути произведено преобразование $k \rightarrow k\cdot r$
*/

void calc_k(coordinate x, coordinate y, coordinate z, timevalue t,
            Coordinate sx, Coordinate sy, Coordinate sz,
            Velocity vx, Velocity vy, Velocity vz,
            timevalue t2,
            long double * k, distance * r, distance * nx, distance * ny, distance * nz);

// отношение радиуса Лиенара Вихерта к радиусу
long double klw(coordinate x, coordinate y, coordinate z, timevalue t,
                Coordinate sx, Coordinate sy, Coordinate sz,
                Velocity vx, Velocity vy, Velocity vz);

// Радиус Лиенара Вихерта
long double Rlw(coordinate x, coordinate y, coordinate z, timevalue t,
                Coordinate sx, Coordinate sy, Coordinate sz,
                Velocity vx, Velocity vy, Velocity vz);
// phi_lw - скалярный потенциал Лиенара Вихерта
long double philw(coordinate x, coordinate y, coordinate z, timevalue t,
           Coordinate sx, Coordinate sy, Coordinate sz,
           Velocity vx, Velocity vy, Velocity vz,
           charge q);


// A_lw - векторный потенциал Лиенара Вихерта
void Alw(coordinate x, coordinate y, coordinate z, timevalue t,
       Coordinate sx, Coordinate sy, Coordinate sz,
       Velocity vx, Velocity vy, Velocity vz,
       charge q,
       field * A_x, field * A_y, field * A_z
       );

void electr_magnet(coordinate x, coordinate y, coordinate z, timevalue t,
                Coordinate sx, Coordinate sy, Coordinate sz,
                Velocity vx, Velocity vy, Velocity vz,
                Acceleration wx, Acceleration wy, Acceleration wz,
                charge q,
                field * E_x, field * E_y, field * E_z,
                field * B_x, field * B_y, field * B_z);

int electr_magnet_ex(coordinate x, coordinate y, coordinate z, timevalue t,
                   Coordinate sx, Coordinate sy, Coordinate sz,
                   Velocity vx, Velocity vy, Velocity vz,
                   Acceleration wx, Acceleration wy, Acceleration wz,
                   charge q,
                   field * E_x, field * E_y, field * E_z,
                   field * B_x, field * B_y, field * B_z,
                   field * A_x, field * A_y, field * A_z);