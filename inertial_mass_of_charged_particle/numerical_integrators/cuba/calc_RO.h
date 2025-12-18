
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <float.h>

#ifdef INCLUDE_CUBA_H
#if REALSIZE == 16
#include "cubaq.h"
#elif REALSIZE == 10
#include "cubal.h"
#else
#include "cuba.h"
#endif
#else
typedef double cubareal;
#endif /*INCLUDE_CUBA_H*/

/*
Sq = function('Sq') # x^2
Cb = function('Cb') # x^3
Qu = function('Qu') # x^4
Fi = function('Fi') # x^5
Si = function('Si') # x^6
Se = function('Si') # x^7
Ei = function('Ei') # x^8
Ni = function('Ni') # x^9
*/
static inline cubareal Sq(cubareal x) {
  return x*x;
}

static inline cubareal Cb(cubareal x) {
  return x*x*x;
}

static inline cubareal Qu(cubareal x) {
  return x*x*x*x;
}

static inline cubareal Fi(cubareal x) {
  return x*x*x*x*x;
}

static inline cubareal Si(cubareal x) {
  return x*x*x*x*x*x;
}

static inline cubareal Se(cubareal x) {
  return x*x*x*x*x*x*x;
}

static inline cubareal Ei(cubareal x) {
  return x*x*x*x*x*x*x*x;
}

static inline cubareal Ni(cubareal x) {
  return x*x*x*x*x*x*x*x*x;
}

static inline cubareal Te(cubareal x) {
  return x*x*x*x*x*x*x*x*x*x;
}

static inline cubareal El(cubareal x) {
  return x*x*x*x*x*x*x*x*x*x*x;
}

static inline cubareal Tw(cubareal x) {
  return x*x*x*x*x*x*x*x*x*x*x*x;
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
    return sqrt(Sq(rca(ra, theta_a)) +
                Sq(rcq(rq, theta_q)) +
                Sq(za(ra, theta_a)-zq(rq, theta_q))
                - 2*rca(ra, theta_a)*rcq(rq, theta_q)*cos(phi_q));
}

// \section{Приближение малых скоростей с учётом запаздывания}
// Решение этой системы имеет весьма сложный вид,
// но если мы исследуем вопрос какова будет инертная масса покоя,
// то при решении этой системы мы можем положить $v = 0$.
// В этом случае для нахождения запаздывания нужно будет решить уравнение

/*$-\frac{1}{4} \, a^{2} \mathit{dt}^{4} + c^{2} \mathit{dt}^{2} - a \mathit{dt}^{2} {\left(z_{a} - z_{q}\right)} - R_{0}^{2} = 0$

где
$(t-t') = dt$

Это уравнение имеет 4 решения, но физически приемлемый смысл при положительном ускорении имеет решение


% $\mathit{dt} = -\frac{\sqrt{2 \, c^{2} - 2 \, a \mathit{dz} + 2 \, \sqrt{-R_{0}^{2} a^{2} + c^{4} - 2 \, a c^{2} \mathit{dz} + a^{2} \mathit{dz}^{2}}}}{a}$

% $\mathit{dt} = \frac{\sqrt{2 \, c^{2} - 2 \, a \mathit{dz} + 2 \, \sqrt{-R_{0}^{2} a^{2} + c^{4} - 2 \, a c^{2} \mathit{dz} + a^{2} \mathit{dz}^{2}}}}{a}$

% $\mathit{dt} = -\frac{\sqrt{2 \, c^{2} - 2 \, a \mathit{dz} - 2 \, \sqrt{-R_{0}^{2} a^{2} + c^{4} - 2 \, a c^{2} \mathit{dz} + a^{2} \mathit{dz}^{2}}}}{a}$

$\mathit{dt} = \frac{\sqrt{2 \, c^{2} - 2 \, a \mathit{dz} - 2 \, \sqrt{-R_{0}^{2} a^{2} + c^{4} - 2 \, a c^{2} \mathit{dz} + a^{2} \mathit{dz}^{2}}}}{a}$

где
$\mathit{dz} = z_{a'} - z_{q'}$


В приближении малых скоростей ${}^{v}/{}_{c}\ll 1$ но при учете запаздывания
\[\overrightarrow{E}=\int\limits_{{{r}_{q}}}\int\limits_{{{\varphi}_{q}}}\int\limits_{{{\theta}_{q}}}\\
{\left\{ -\frac{\overrightarrow{a}R}{{{c}^{2}}} \right\}\frac{\rho \left( {{r}_{q}} \right){{r}_{q}}^{2}\sin \left( {{\theta }_{q}} \right)}{{{R}^{*}}^{2}}\ }d{{\theta }_{q}}d{{\varphi }_{q}}d{{r}_{q}}\]
 Откуда
\[{{F}_{z}}=-\frac{\overrightarrow{a}}{{{c}^{^{2}}}}\int\limits_{{{V}_{a}}}{\int\limits_{{{V}_{q}}}{\frac{\rho \left( {{r}_{q}} \right)\rho \left( {{r}_{a}} \right)}{R}}}\ d{{V}_{q}}d{{V}_{a}}\]
где
$\mathit{R} = c\frac{\sqrt{2 \, c^{2} - 2 \, a \mathit{dz} - 2 \, \sqrt{-R_{0}^{2} a^{2} + c^{4} - 2 \, a c^{2} \mathit{dz} + a^{2} \mathit{dz}^{2}}}}{a}$
*/
#if 0

static inline cubareal R_a (cubareal ra, cubareal theta_a, cubareal rq, cubareal theta_q, cubareal phi_q, cubareal c, cubareal a)
{
    cubareal R_0 = R0 (ra, theta_a, rq, theta_q, phi_q);
    cubareal dz = za(ra, theta_a) - zq(rq, theta_q);

    //cubareal tt1 = -2*a*dz + 2*sqrt(-2*a*dz*Sq(c) - Sq(R_0)*Sq(a) + Sq(a)*Sq(dz) + Qu(c)) + 2*Sq(c);
    //cubareal tt2 = -2*a*dz - 2*sqrt(-2*a*dz*Sq(c) - Sq(R_0)*Sq(a) + Sq(a)*Sq(dz) + Qu(c)) + 2*Sq(c);
    /*
    res3 = c*dt == -sqrt(-2*a*dz + 2*sqrt(-2*a*dz*Sq(c) - Sq(R_0)*Sq(a) + Sq(a)*Sq(dz) + Qu(c)) + 2*Sq(c))*c/a
    res3 = c*dt ==  sqrt(-2*a*dz + 2*sqrt(-2*a*dz*Sq(c) - Sq(R_0)*Sq(a) + Sq(a)*Sq(dz) + Qu(c)) + 2*Sq(c))*c/a
    res3 = c*dt == -sqrt(-2*a*dz - 2*sqrt(-2*a*dz*Sq(c) - Sq(R_0)*Sq(a) + Sq(a)*Sq(dz) + Qu(c)) + 2*Sq(c))*c/a
    res3 = c*dt ==  sqrt(-2*a*dz - 2*sqrt(-2*a*dz*Sq(c) - Sq(R_0)*Sq(a) + Sq(a)*Sq(dz) + Qu(c)) + 2*Sq(c))*c/a
    */
    cubareal R = c/a*sqrt(2*Sq(c) - 2*a*dz - 2*sqrt(-Sq(R_0*a) + Qu(c) - 2*a*Sq(c)*dz + Sq(a*dz) ) );

    //printf("tt1 = %f, tt2 = %f, R = %f R0 = %f , ra = %f, theta_a = %f, rq = %f, theta_q = %f, phi_q = %f\n", tt1, tt2, R, R_0, ra, theta_a, rq, theta_q, phi_q);
    return R;
}

#else

static inline cubareal R_a/*_qwen*/(cubareal ra, cubareal theta_a, cubareal rq, cubareal theta_q, cubareal phi_q, cubareal c, cubareal a)
{
    cubareal xa = ra * sin(theta_a) * cos(0.0);
    cubareal ya = ra * sin(theta_a) * sin(0.0);
    cubareal za = ra * cos(theta_a);

    cubareal xq = rq * sin(theta_q) * cos(phi_q);
    cubareal yq = rq * sin(theta_q) * sin(phi_q);
    cubareal zq = rq * cos(theta_q);

    cubareal dx = xa - xq;
    cubareal dy = ya - yq;
    cubareal dz = za - zq;  // мгновенное z_a - z_q при t=0
    cubareal R_perp2 = dx*dx + dy*dy;

    // D = c^4 + 2*a*c^2*dz - a^2 * R_perp2
    cubareal D = Sq(Sq(c)) + 2.0 * a * Sq(c) * dz - Sq(a) * R_perp2;

    if (D < 0.0) {
        // fallback: разложение до O(a)
        cubareal R0 = sqrt(dx*dx + dy*dy + dz*dz);
        return R0 + a * dz * R0 / (2.0 * Sq(c));  // dt ≈ R0/c + a dz R0/(2c^3)
    }

    cubareal sqrtD = sqrt(D);
    // Физический корень: малый tau
    cubareal tau2 = (2.0 * (Sq(c) + a * dz) - 2.0 * sqrtD) / Sq(a);
    if (tau2 <= 0) return DBL_MAX;

    cubareal tau = sqrt(tau2);
    return c * tau;  // R = c * tau
}
#endif

static inline cubareal R_v (cubareal ra, cubareal theta_a, cubareal rq, cubareal theta_q, cubareal phi_q, cubareal c, cubareal v)
{
    cubareal R_0 = R0 (ra, theta_a, rq, theta_q, phi_q);
    cubareal dz = za(ra, theta_a) - zq(rq, theta_q);
    /*
    res4 =  [
    dt == (dz*v - sqrt(R_0^2*c^2 - (R_0^2 - dz^2)*v^2))/(c^2 - v^2),
    dt == (dz*v + sqrt(R_0^2*c^2 - (R_0^2 - dz^2)*v^2))/(c^2 - v^2)
    ]
    c*dt == (dz*v - sqrt(Sq(R_0)*Sq(c) - (Sq(R_0) - Sq(dz))*Sq(v)))*c/(Sq(c) - Sq(v))
    */
    cubareal R = (dz*v - sqrt(Sq(R_0)*Sq(c) - (Sq(R_0) - Sq(dz))*Sq(v)))*c/(Sq(c) - Sq(v));

    printf("R = %f R0 = %f , ra = %f, theta_a = %f, rq = %f, theta_q = %f, phi_q = %f\n", R, R_0, ra, theta_a, rq, theta_q, phi_q);
    return R;
}
