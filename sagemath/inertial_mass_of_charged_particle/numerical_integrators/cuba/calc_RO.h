
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#if REALSIZE == 16
#include "cubaq.h"
#elif REALSIZE == 10
#include "cubal.h"
#else
#include "cuba.h"
#endif


static inline cubareal Sq(cubareal x) {
  return x*x;
}

static inline cubareal Cb(cubareal x) {
  return x*x*x;
}

//z - координаты заряда и точки наблюдения
static inline cubareal zq (cubareal rq, cubareal theta_q)
{
    return (rq)*cos(theta_q);
}
static inline cubareal za (cubareal ra, cubareal theta_a)
{
    return (ra)*cos(theta_a);
}

// введём вспомогательные переменные - цилиндрический радиус
// как координата r точки при переходе в цилиндрическую систему коорденат с тем же направлением оси z
static inline cubareal rcq (cubareal rq, cubareal theta_q )
{
    return (rq)*sin(theta_q);
}
static inline cubareal rca (cubareal ra, cubareal theta_a )
{
    return (ra)*sin(theta_a);
}

// выражение для расстояния между точкой заряда и точкой наблюдения примет вид
static inline cubareal R0 (cubareal ra, cubareal theta_a, cubareal rq, cubareal theta_q, cubareal phi_q)
{
    return sqrt(Sq(rca(ra, theta_a)) + Sq(rcq(rq, theta_q)) + Sq(za(ra, theta_a)-zq(rq, theta_q)) - 2*rca(ra, theta_a)*rcq(rq, theta_q)*cos(phi_q));
}
