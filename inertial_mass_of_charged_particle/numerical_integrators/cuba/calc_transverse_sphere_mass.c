/*
    calc_transverse_sphere_mass.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef INCLUDE_CUBA_H
#include "cuba.h"
#endif

//#define Z_Z
#define Z_RHO
//#define Z_PHI

static inline cubareal Sq(cubareal x) {
  return x*x;
}

static inline cubareal Cb(cubareal x) {
  return x*x*x;
}

/* Структура для передачи параметров задачи */
typedef struct {
    double R0;     /* Радиус сферы */
    double rho0;   /* Радиус орбиты (расстояние от центра до сферы) */
    double omega;  /* Угловая скорость вращения */
    double c;      /* Скорость света */
    int use_delay;
    int use_lorentz_factor;
    int use_lorentz_general_factor;
    int use_fermi_factor_O;
    int use_fermi_factor;
    int use_fermi_general_factor;
    int use_fast_integrand;
} ProblemParams;

/* Функция для вычисления поправки Ферми по отношению к опорной точке О*/
static double fermi_correction_O(
    double R_rho_O,
    double R_phi_O,
    double R_z_O,
    double a_rho,
    double a_phi,
    double a_z,
    double c
) {
    /* Скалярное произведение R · a */
    double R_dot_a_O = R_rho_O * a_rho + R_phi_O * a_phi + R_z_O * a_z;

    /* Поправка Ферми */
    return (/*1.0*/ + R_dot_a_O / (c * c));
}

/* Функция для вычисления поправки Ферми */
static double fermi_correction(
    double R_rho,
    double R_phi,
    double R_z,
    double a_rho,
    double a_phi,
    double a_z,
    double c
) {
    /* Скалярное произведение R · a */
    double R_dot_a = R_rho * a_rho + R_phi * a_phi + R_z * a_z;

    /* Поправка Ферми */
    return (/*1.0*/ + 0.5 * R_dot_a / (c * c));
}

/* Функция для вычисления поправки Ферми в общем релятивистском случае */
static double fermi_correction_general(
    double t,           /* текущее время */
    double t_prime,     /* запаздывающее время */
    double R_rho,
    double R_phi,
    double R_z,
    double a_rho,
    double a_phi,
    double a_z,
    double v_rho,
    double v_phi,
    double v_z,
    double c
) {
    /* В 4-пространстве необходимо преобразовать цилиндрические координаты в декартовы */
    double R_x = R_rho * cos(R_phi);
    double R_y = R_rho * sin(R_phi);
    double a_x = a_rho * cos(R_phi) - a_phi * sin(R_phi);
    double a_y = a_rho * sin(R_phi) + a_phi * cos(R_phi);
    double v_x = v_rho * cos(R_phi) - v_phi * sin(R_phi);
    double v_y = v_rho * sin(R_phi) + v_phi * cos(R_phi);

    /* Вычисляем скорость */
    double v = sqrt(v_x*v_x + v_y*v_y + v_z*v_z);

    /* Фактор Лоренца */
    double gamma = 1.0 / sqrt(1.0 - v*v/(c*c));
    double gamma2 = gamma * gamma;
    double gamma4 = gamma2 * gamma2;

    /* Скалярное произведение v·a */
    double v_dot_a = v_x*a_x + v_y*a_y + v_z*a_z;

    /* 4-радиус-вектор */
    double R0 = c * (t - t_prime);  /* временная компонента */
    double R1 = R_x;
    double R2 = R_y;
    double R3 = R_z;

    /* 4-ускорение */
    double a0 = gamma4 * v_dot_a / c;
    double a1 = gamma2 * a_x + gamma4 * v_dot_a * v_x / (c*c);
    double a2 = gamma2 * a_y + gamma4 * v_dot_a * v_y / (c*c);
    double a3 = gamma2 * a_z + gamma4 * v_dot_a * v_z / (c*c);

    /* 4-скалярное произведение */
    double R_dot_a = R0 * a0 - R1 * a1 - R2 * a2 - R3 * a3;

    /* Поправка Ферми в общем случае */
    return (/*1.0*/ - 0.5 * R_dot_a / (c * c));
}

static void computing_electric_field(
    double t_prime,     /* запаздывающее время */
    double R,
    double R_rho_O, double R_phi_O, double R_z_O,
    double R_rho, double R_phi, double R_z,
    double v_phi,
    double a_rho,
    double c,
    double *E1_rho, double *E1_phi, double *E1_z,
    double *E2_rho, double *E2_phi, double *E2_z,
    double *E_total_rho, double *E_total_phi, double *E_total_z,
    int use_fermi_factor_O,
    int use_fermi_factor,
    int use_fermi_general_factor,
    int use_fast_integrand
)
{
    /* Радиус Лиенара-Вихерта (учет запаздывания) */
    double R_star = R - (R_phi * v_phi) / c;

    /* Защита от деления на ноль */
    if (R_star < 1e-10 || R < 1e-10) {
        *E1_rho = *E1_phi = *E1_z = 0.0;
        *E2_rho = *E2_phi = *E2_z = 0.0;
        *E_total_rho = *E_total_phi = *E_total_z = 0.0;
        return;
    }

    double fermi_factor = 0.0;

    double a_phi = 0.0;
    double a_z = 0.0;

    double v_rho = 0.0;
    double v_z = 0.0;

    if (use_fermi_general_factor)
    {
        double t = 0;    /* текущее время */

        fermi_factor = fermi_correction_general(
            t,           /* текущее время */
            t_prime,     /* запаздывающее время */
            R_rho, R_phi, R_z,
            a_rho, a_phi, a_z,
            v_rho, v_phi, v_z,
            c);
    }
    else if (use_fermi_factor)
    {
        fermi_factor = fermi_correction(
            R_rho, R_phi, R_z,
            a_rho, a_phi, a_z,
            c);
    }
    else if (use_fermi_factor_O)
    {
        fermi_factor = fermi_correction_O(
            R_rho_O, R_phi_O, R_z_O,
            a_rho, a_phi, a_z,
            c);
    }
    double v2_c2 = (Sq(v_rho) + Sq(v_phi) + Sq(v_z) ) / Sq(c);
    double Ra_c2 = (R_rho * a_rho) / Sq(c);

    double f_one = 1.0;
    if (use_fast_integrand)
    {
        /*
            Режим быстрого интегранда выключает из процедуры численного интегрирования
            симметричную часть интегранда
            $$R_z / Cb(R - (R_z * v_z) / c)$$
            которая теоретически при интегрировании по сфере должна дать значение близкое к нулю.
        */

        /*
            двойной интеграл:
            $$I = \int \int \frac{v_z}{R^{*2}} \, dV_q \, dV_a$$

            Где $R^* = R - \frac{v_z}{c} R_z$,
             $R_z = z_a - z_q$,
             $R = \sqrt{(x_a-x_q)^2 + (y_a-y_q)^2 + (z_a-z_q)^2}$.

             равен нулю из-за симметрии сферы

       */
        f_one = 0.0;
    }

    /* Вычисление градиентного поля E1 */
    double common_factor1 = 1.0 / pow(R_star, 2);
#if 0
    double velocity_factor1 = /*1.0 +*/ Ra_c2 /*- v2_c2*/;

    *E1_rho = common_factor1 * (R_rho * velocity_factor1 / R_star);
    *E1_phi = common_factor1 * (R_phi * velocity_factor1 / R_star/* - v_phi / c*/);
    *E1_z   = common_factor1 * (R_z   * velocity_factor1 / R_star);
#else
    /* Вычисление градиентного поля E1 (только слагаемое с ускорением) */
    /* (1.0 + Ra_c2 - v2_c2) * (1.0+fermi_factor)  = 1.0 + Ra_v2_c2_fermi */
    double Ra_v2_c2_fermi = fermi_factor + Ra_c2 + fermi_factor * (Ra_c2 - v2_c2);

    *E1_rho = common_factor1 * (R_rho * (f_one + Ra_v2_c2_fermi) / R_star - f_one * v_rho / c - v_rho / c * fermi_factor);
    *E1_phi = common_factor1 * (R_phi * (f_one + Ra_v2_c2_fermi) / R_star - f_one * v_phi / c - v_phi / c * fermi_factor);
    *E1_z   = common_factor1 * (R_z   * (f_one + Ra_v2_c2_fermi) / R_star - f_one * v_z / c   - v_z   / c * fermi_factor);
#endif

    /* Вычисление поля самоиндукции E2 */
    double common_factor2   = 1.0 / pow(R_star, 2);
    double velocity_factor2 = (R / R_star) * (v2_c2 - Ra_c2 - 1.0) + 1.0;

    common_factor2 *= (1.0+fermi_factor);

    *E2_rho = common_factor2 * (- a_rho * R / pow(c, 2));
    *E2_phi = common_factor2 * (v_phi / c * velocity_factor2);
    *E2_z = 0;

    /* Вычисление суммарного поля */
    double common_factor = 1.0 / pow(R_star, 3);
#if 0
    double velocity_factor = /*1.0 +*/ Ra_c2 /*- v2_c2*/;

    *E_total_rho = common_factor * ((R_rho                  ) * velocity_factor - (a_rho * R_star * R) / pow(c, 2));
    *E_total_phi = common_factor * ((R_phi - (R * v_phi) / c) * velocity_factor);
    *E_total_z   = common_factor * ((R_z                    ) * velocity_factor);
#else
    *E_total_rho = common_factor * ( ((f_one + Ra_v2_c2_fermi) * R_rho - (f_one + Ra_v2_c2_fermi) * R*v_rho/c)  - (a_rho * R_star * R) / Sq(c) - (a_rho * R_star * R) * fermi_factor / Sq(c));
    *E_total_phi = common_factor * ( ((f_one + Ra_v2_c2_fermi) * R_phi - (f_one + Ra_v2_c2_fermi) * R*v_phi/c)  - (a_phi * R_star * R) / Sq(c) - (a_phi * R_star * R) * fermi_factor / Sq(c));
    *E_total_z   = common_factor * ( ((f_one + Ra_v2_c2_fermi) * R_z   - (f_one + Ra_v2_c2_fermi) * R*v_z/c)    - (a_z   * R_star * R) / Sq(c) - (a_z   * R_star * R) * fermi_factor / Sq(c));
#endif
}

/* Функция для вычисления фактора Лоренца */
static inline double lorentz_factor(double v, double c) {
    double beta = v / c;
    return 1.0 / sqrt(1.0 - beta * beta);
}

static void compute_electric_field(
    double ra, double theta_a, double psi_a, /*воображаемые координаты несплюснутой сферы*/
    double rq, double theta_q, double psi_q, /*по которым производится интегрирование*/
    const ProblemParams *params,
    double *E1_rho, double *E1_phi, double *E1_z,
    double *E2_rho, double *E2_phi, double *E2_z,
    double *E_total_rho, double *E_total_phi, double *E_total_z
) {
    /* Преобразование в декартовы координаты в локальной системе */
#ifdef Z_Z
    /* Преобразование в декартовы координаты в локальной системе */
    double x_a = params->rho0 + ra * sin(theta_a) * cos(psi_a);
    double y_a = ra * sin(theta_a) * sin(psi_a);
    double z_a = ra * cos(theta_a);

    double x_q = params->rho0 + rq * sin(theta_q) * cos(psi_q);
    double y_q = rq * sin(theta_q) * sin(psi_q);
    double z_q = rq * cos(theta_q);
#endif

#ifdef Z_RHO
    /* Ось z = r * cos(theta) сферической системы связанной с зарядом
        направлена вдоль радиального направления
        лабораторной цилиндрической системы координат
        Поэтому x_a и y_a - это лабораторные декартовы координаты в момент t = 0
     */
    double x_a = params->rho0 + ra * cos(theta_a);
    double y_a = ra * sin(theta_a) * sin(psi_a);
    double z_a = ra * sin(theta_a) * cos(psi_a);

    double x_q = params->rho0 + rq * cos(theta_q);
    double y_q = rq * sin(theta_q) * sin(psi_q);
    double z_q = rq * sin(theta_q) * cos(psi_q);

    double x_O = params->rho0;
    double y_O = 0.0;
    double z_O = 0.0;
#endif

#ifdef Z_PHI
    /* Ось z = r * cos(theta) сферической системы связанной с зарядом
        направлена азимутально (вдоль phi)
        сонаправлено скорости - для упрощения расчёта изменения
        формы сферы в результате лоренцева сокращения
     */
    double x_a = params->rho0 + ra * sin(theta_a) * cos(psi_a);
    double y_a = ra * cos(theta_a);
    double z_a = ra * sin(theta_a) * sin(psi_a);

    double x_q = params->rho0 + rq * sin(theta_q) * cos(psi_q);
    double y_q = rq * cos(theta_q);
    double z_q = rq * sin(theta_q) * sin(psi_q);
#endif

    /*
        Переменные x_a, y_a и z_a относятся к лабораторной
        цилиндрической системе следующим образом

        z_a - как раз соответствует оси z лабораторной
        цилиндрической системы координат,

        в рассматриваемый текущий момент времени t = 0
        центр сферы находится в точке
        X = params->rho0, Y = 0 - лабораторной декартовой системы.

        Поэтому x_a и y_a - это лабораторные декартовы координаты в момент t = 0
    */

    /* Преобразование в цилиндрическую систему (исходная) */
    double rho_a = sqrt(pow(x_a, 2) + pow(y_a, 2));
    double phi_a = atan2(y_a, x_a);
    if (phi_a < 0) phi_a += 2.0 * M_PI;

    double rho_q = sqrt(pow(x_q, 2) + pow(y_q, 2));
    double phi_q = atan2(y_q, x_q);
    if (phi_q < 0) phi_q += 2.0 * M_PI;

    // double rho_O = sqrt(pow(x_O, 2) + pow(y_O, 2));
    // double phi_O = atan2(y_O, x_O);
    // if (phi_O < 0) phi_O += 2.0 * M_PI;

    double rho_O = x_O;
    double phi_O = 0.0;

    /* Скорость и ускорение источника */
    double v_phi = rho_q * params->omega;
    double a_rho = -rho_q * pow(params->omega, 2);

    if (params->use_lorentz_factor || params->use_lorentz_general_factor)
    {
        /* Учет Лоренц-сокращения в направлении движения */
        double gamma = lorentz_factor(v_phi, params->c);
        /* Учет Лоренц-сокращения: в системе отсчета наблюдателя
        сфера сжата вдоль направления движения
        (ось z - локальной сферической системы)
        (она же ось y - после преобразований в декартовы координаты в локальной системе)
         */
        y_q = y_q / gamma;
        y_a = y_a / gamma;

        /* Преобразование в цилиндрическую систему (исходная) */
        rho_a = sqrt(pow(x_a, 2) + pow(y_a, 2));
        phi_a = atan2(y_a, x_a);
        if (phi_a < 0) phi_a += 2.0 * M_PI;

        rho_q = sqrt(pow(x_q, 2) + pow(y_q, 2));
        phi_q = atan2(y_q, x_q);
        if (phi_q < 0) phi_q += 2.0 * M_PI;

        v_phi = rho_q * params->omega;
        a_rho = -rho_q * pow(params->omega, 2);
    }

    /* Вектор от точки О к наблюдателю */
    double R_rho_O = rho_a * cos(phi_a - phi_O) - rho_O;
    double R_phi_O = rho_a * sin(phi_a - phi_O);
    double R_z_O = z_a - z_O;

    /* Вектор от источника к наблюдателю */
    double R_rho = rho_a * cos(phi_a - phi_q) - rho_q;
    double R_phi = rho_a * sin(phi_a - phi_q);
    double R_z = z_a - z_q;
    double R = sqrt(pow(R_rho, 2) + pow(R_phi, 2) + pow(R_z, 2));

    double t_prime = 0.0;

    computing_electric_field(
        t_prime,     /* запаздывающее время */
        R,
        R_rho_O, R_phi_O, R_z_O,
        R_rho, R_phi, R_z,
        v_phi,
        a_rho,
        params->c,
        E1_rho, E1_phi, E1_z,
        E2_rho, E2_phi, E2_z,
        E_total_rho, E_total_phi, E_total_z,
        params->use_fermi_factor_O,
        params->use_fermi_factor,
        params->use_fermi_general_factor,
        params->use_fast_integrand
    );
}

/* Функция для вычисления координат источника в запаздывающий момент времени */
static void get_source_position(
    double t_prime,
    double rq,
    double theta_q,
    double psi_q,
    double rho0,
    double omega,
    double c,
    double *x_q,
    double *y_q,
    double *z_q,
    int use_lorentz_factor,
    int use_lorentz_general_factor
) {
    /* Координаты точки на сфере в системе покоя сферы */
#ifdef Z_Z
    /* Преобразование в декартовы координаты в локальной системе */
    double x_q0 = rho0 + rq * sin(theta_q) * cos(psi_q);
    double y_q0 = rq * sin(theta_q) * sin(psi_q);
    double z_q0 = rq * cos(theta_q);
#endif

#ifdef Z_RHO
    /* Ось z = r * cos(theta) сферической системы связанной с зарядом
        направлена вдоль радиального направления
        лабораторной цилиндрической системы координат
        Поэтому x_a и y_a - это лабораторные декартовы координаты в момент t = 0
     */
    double x_q0 = rho0 + rq * cos(theta_q);
    double y_q0 = rq * sin(theta_q) * sin(psi_q);
    double z_q0 = rq * sin(theta_q) * cos(psi_q);
#endif

#ifdef Z_PHI
    /* Преобразование в декартовы координаты в локальной системе */
    /* Ось z сферической системы направлена азимутально (вдоль phi) */
    /* - сонаправлено скорости - для упрощения расчёта изменения */
    /* формы сферы в результате лоренцева сокращения */

    double x_q0 = rho0 + rq * sin(theta_q) * cos(psi_q);
    double z_q0 = rq * sin(theta_q) * sin(psi_q);
    double y_q0 = rq * cos(theta_q);
#endif

    double rho_q = sqrt(pow(x_q0, 2) + pow(y_q0, 2));
    if (use_lorentz_factor || use_lorentz_general_factor)
    {
    	/* Скорость источника */
    	double v_phi = rho_q * omega;
        /* Учет Лоренц-сокращения в направлении движения */
        double gamma = lorentz_factor(v_phi, c);
        /* Учет Лоренц-сокращения: в системе отсчета наблюдателя
        сфера сжата вдоль направления движения (ось z) */
        y_q0 /= gamma;
#if 0
        if (use_lorentz_general_factor) {
            /* Учет дополнительного искривления из-за ускорения */
            double proper_time = compute_proper_time(t_prime, t0, v0, a, c);
            double acceleration_correction = 0.5 * a * pow(proper_time, 2) / gamma;
            *y_q0 += acceleration_correction;
        }
#endif
        rho_q = sqrt(pow(x_q0, 2) + pow(y_q0, 2));
    }
    double phi_q = atan2(y_q0, x_q0);
    if (phi_q < 0) phi_q += 2.0 * M_PI;

    /* Угловая позиция источника в запаздывающий момент времени */
    double phi_prime = omega * t_prime;

    /* Поворот системы координат сферы на угол phi_prime */
    phi_q += phi_prime;

    *x_q = rho_q * sin(phi_q);
    *y_q = rho_q * cos(phi_q);
    *z_q = z_q0;
}

/* Функция для вычисления R(t') - расстояния от источника к наблюдателю в момент t' */
static double compute_R(
    double t_prime,
    double ra,
    double theta_a,
    double psi_a,
    double rq,
    double theta_q,
    double psi_q,
    double rho0,
    double omega,
    double c,
    int use_lorentz_factor,
    int use_lorentz_general_factor
) {
    /* Координаты наблюдателя в момент 0*/
    double x_a, y_a, z_a;
    get_source_position(0, ra, theta_a, psi_a, rho0, omega, c, &x_a, &y_a, &z_a,
        use_lorentz_factor,
        use_lorentz_general_factor);

    /* Координаты источника в момент t' */
    double x_q, y_q, z_q;
    get_source_position(t_prime, rq, theta_q, psi_q, rho0, omega, c, &x_q, &y_q, &z_q,
        use_lorentz_factor,
        use_lorentz_general_factor);

    /* Расстояние */
    double dx = x_a - x_q;
    double dy = y_a - y_q;
    double dz = z_a - z_q;
    return sqrt(dx*dx + dy*dy + dz*dz);
}

/* Аналитическое вычисление dR/dt' с использованием Лоренц-фактора */
static double compute_dR_dt_prime_analytical(
    double t_prime,
    double ra,
    double theta_a,
    double psi_a,
    double rq,
    double theta_q,
    double psi_q,
    double rho0,
    double omega,
    double c,
    int use_lorentz_factor,
    int use_lorentz_general_factor
) {
    /* Координаты наблюдателя в момент 0*/
    double x_a, y_a, z_a;
    get_source_position(0, ra, theta_a, psi_a, rho0, omega, c, &x_a, &y_a, &z_a,
        use_lorentz_factor,
        use_lorentz_general_factor);

    /* Координаты источника в запаздывающий момент времени */
    double x_q, y_q, z_q;
    get_source_position(t_prime, rq, theta_q, psi_q, rho0, omega, c, &x_q, &y_q, &z_q,
        use_lorentz_factor,
        use_lorentz_general_factor);
    double rho_q = sqrt(pow(x_q, 2) + pow(y_q, 2));

    /* Вектор от источника к наблюдателю */
    double R_x = x_a - x_q;
    double R_y = y_a - y_q;
    double R_z = z_a - z_q;
    double R = sqrt(R_x*R_x + R_y*R_y + R_z*R_z);

    /* Скорость источника в запаздывающий момент времени */
    double phi_prime = omega * t_prime;
    double v_x = rho_q * omega * cos(phi_prime);
    double v_y = -rho_q * omega * sin(phi_prime);
    double v_z = 0.0;
    double v = sqrt(v_x*v_x + v_y*v_y + v_z*v_z);

    double gamma = lorentz_factor(v, c);
    double v_r = v / gamma;  /* Радиальная скорость с учетом Лоренц-сокращения */

    /* Вычисляем проекцию скорости на направление к наблюдателю */
    double cos_theta = (R_x * v_x + R_y * v_y + R_z * v_z) / (R * v);
    double v_radial = v_r * cos_theta;

    /* Аналитическое выражение для dR/dt' */
    return -v_radial;
}

/* Улучшенное вычисление запаздывающего времени с методом Ньютона */
static double compute_retarded_time_newton(
    double ra,
    double theta_a,
    double psi_a,
    double rq,
    double theta_q,
    double psi_q,
    double rho0,
    double omega,
    double c,
    int max_iter,
    double tolerance,
    int use_lorentz_factor,
    int use_lorentz_general_factor
) {
    /* Координаты наблюдателя в момент 0*/
    double x_a, y_a, z_a;
    get_source_position(0, ra, theta_a, psi_a, rho0, omega, c, &x_a, &y_a, &z_a,
        use_lorentz_factor,
        use_lorentz_general_factor);

    /* Координаты источника в момент времени 0 */
    double x_q, y_q, z_q;
    get_source_position(0, rq, theta_q, psi_q, rho0, omega, c, &x_q, &y_q, &z_q,
        use_lorentz_factor,
        use_lorentz_general_factor);

    double R = sqrt(pow(x_a - x_q, 2) + pow(y_a - y_q, 2) + pow(z_a - z_q, 2));
    double t = 0.0;
    double t_prime = t - R / c;

    /* Итерационное решение */
    for (int i = 0; i < max_iter; i++) {
        get_source_position(t_prime, rq, theta_q, psi_q, rho0, omega, c, &x_q, &y_q, &z_q,
            use_lorentz_factor,
            use_lorentz_general_factor);

        R = sqrt(pow(x_a - x_q, 2) + pow(y_a - y_q, 2) + pow(z_a - z_q, 2));
        double t_prime_new = t - R / c;

        if (fabs(t_prime_new - t_prime) < tolerance) {
            t_prime = t_prime_new;
            break;
        }

        t_prime = t_prime_new;
    }

    /* Уточнение методом Ньютона с аналитической производной */
    for (int i = 0; i < 5; i++) {
        double R_val = compute_R(t_prime, ra, theta_a, psi_a, rq, theta_q, psi_q, rho0, omega, c,
            use_lorentz_factor,
            use_lorentz_general_factor);
        double f = t_prime - t + R_val / c;
        double df_dt = 1.0 + compute_dR_dt_prime_analytical(t_prime, ra, theta_a, psi_a, rq, theta_q, psi_q, rho0, omega, c,
            use_lorentz_factor,
            use_lorentz_general_factor) / c;

        double delta = -f / df_dt;
        t_prime += delta;

        if (fabs(delta) < tolerance) {
            break;
        }
    }

    return t_prime;
}

/* Функция для вычисления электрического поля с учетом запаздывания */
static void compute_electric_field_with_delay(
    double ra, double theta_a, double psi_a,
    double rq, double theta_q, double psi_q,
    const ProblemParams *params,
    double *E1_rho, double *E1_phi, double *E1_z,
    double *E2_rho, double *E2_phi, double *E2_z,
    double *E_total_rho, double *E_total_phi, double *E_total_z
) {
    /* Вычисление запаздывающего времени точки О */
    double tO_prime = compute_retarded_time_newton(
        ra, theta_a, psi_a, 0, 0, 0,
        params->rho0, params->omega, params->c, 100, 1e-10,
        params->use_lorentz_factor, params->use_lorentz_general_factor
    );

    /* Вычисление запаздывающего времени */
    double t_prime = compute_retarded_time_newton(
        ra, theta_a, psi_a, rq, theta_q, psi_q,
        params->rho0, params->omega, params->c, 100, 1e-10,
        params->use_lorentz_factor, params->use_lorentz_general_factor
    );

    /* Координаты наблюдателя в момент t = 0 */
    double x_a, y_a, z_a;
    get_source_position(0, ra, theta_a, psi_a, params->rho0, params->omega, params->c, &x_a, &y_a, &z_a,
        params->use_lorentz_factor,
        params->use_lorentz_general_factor);

    /* Получение координат точки О в запаздывающий момент времени */
    double x_O, y_O, z_O;
    get_source_position(tO_prime, 0, 0, 0,
                       params->rho0, params->omega, params->c, &x_O, &y_O, &z_O,
                       params->use_lorentz_factor, params->use_lorentz_general_factor);

    /* Получение координат источника в запаздывающий момент времени */
    double x_q, y_q, z_q;
    get_source_position(t_prime, rq, theta_q, psi_q,
                       params->rho0, params->omega, params->c, &x_q, &y_q, &z_q,
                       params->use_lorentz_factor, params->use_lorentz_general_factor);

    /* Преобразование в цилиндрическую систему (исходная) */
    double rho_a = sqrt(pow(x_a, 2) + pow(y_a, 2));
    double phi_a = atan2(y_a, x_a);
    if (phi_a < 0) phi_a += 2.0 * M_PI;

    double rho_q = sqrt(pow(x_q, 2) + pow(y_q, 2));
    double phi_q = atan2(y_q, x_q);
    if (phi_q < 0) phi_q += 2.0 * M_PI;

    double rho_O = sqrt(pow(x_O, 2) + pow(y_O, 2));
    double phi_O = atan2(y_O, x_O);
    if (phi_O < 0) phi_O += 2.0 * M_PI;

    /* Вектор от точки О к наблюдателю */
    double R_rho_O = rho_a * cos(phi_a - phi_O) - rho_O;
    double R_phi_O = rho_a * sin(phi_a - phi_O);
    double R_z_O = z_a - z_O;

    /* Вектор от источника к наблюдателю */
    double R_rho = rho_a * cos(phi_a - phi_q) - rho_q;
    double R_phi = rho_a * sin(phi_a - phi_q);
    double R_z = z_a - z_q;
    double R = sqrt(pow(R_rho, 2) + pow(R_phi, 2) + pow(R_z, 2));

    /* Скорость и ускорение источника */
    double v_phi = rho_q * params->omega;
    double a_rho = -rho_q * pow(params->omega, 2);

    computing_electric_field(
        t_prime,     /* запаздывающее время */
        R,
        R_rho_O, R_phi_O, R_z_O,
        R_rho, R_phi, R_z,
        v_phi,
        a_rho,
        params->c,
        E1_rho, E1_phi, E1_z,
        E2_rho, E2_phi, E2_z,
        E_total_rho, E_total_phi, E_total_z,
        params->use_fermi_factor_O,
        params->use_fermi_factor,
        params->use_fermi_general_factor,
        params->use_fast_integrand
    );
}

// the charge distribution
static inline cubareal rho_q (cubareal r0, cubareal q)
{
    return 3 * (q) / (4*M_PI*(r0)*(r0)*(r0));
}


// интегрирование по координатам заряда источника потенциала
static inline void int_q (const ProblemParams *params, cubareal q, cubareal psi_a, cubareal theta_a, cubareal ra, cubareal psi_q, cubareal theta_q, cubareal rq,
    double *int_q_E1_rho, double *int_q_E1_phi, double *int_q_E1_z,
    double *int_q_E2_rho, double *int_q_E2_phi, double *int_q_E2_z,
    double *int_q_E_rho,  double *int_q_E_phi,  double *int_q_E_z
 )
{
    cubareal r0 = params->R0;

    /* Вычисление электрического поля */
    double E1_rho, E1_phi, E1_z;
    double E2_rho, E2_phi, E2_z;
    double E_rho,  E_phi,  E_z;

    if (params->use_delay)
    {
        compute_electric_field_with_delay(ra, theta_a, psi_a,
                           rq, theta_q, psi_q, params,
                          &E1_rho, &E1_phi, &E1_z,
                          &E2_rho, &E2_phi, &E2_z,
                          &E_rho,  &E_phi,  &E_z);
    }
    else {
        compute_electric_field(ra, theta_a, psi_a,
                           rq, theta_q, psi_q, params,
                          &E1_rho, &E1_phi, &E1_z,
                          &E2_rho, &E2_phi, &E2_z,
                          &E_rho,  &E_phi,  &E_z);
    }

    double k_q = rho_q(r0, q) * Sq(rq) * sin(theta_q);

    *int_q_E1_rho = k_q * E1_rho;
    *int_q_E1_phi = k_q * E1_phi;
    *int_q_E1_z   = k_q * E1_z;
    *int_q_E2_rho = k_q * E2_rho;
    *int_q_E2_phi = k_q * E2_phi;
    *int_q_E2_z   = k_q * E2_z;
    *int_q_E_rho  = k_q * E_rho;
    *int_q_E_phi  = k_q * E_phi;
    *int_q_E_z    = k_q * E_z;
}

// интегрирование по координатам точек наблюдения
void int_a (
    const ProblemParams *params, cubareal q,
    cubareal phi_a, cubareal theta_a, cubareal ra,
    cubareal phi_q, cubareal theta_q, cubareal rq,
    double *int_a_E1_rho, double *int_a_E1_phi, double *int_a_E1_z,
    double *int_a_E2_rho, double *int_a_E2_phi, double *int_a_E2_z,
    double *int_a_E_total_rho, double *int_a_E_total_phi, double *int_a_E_total_z
)
{
    cubareal r0 = params->R0;
    double k_a = rho_q(r0, q) * Sq(ra) * sin(theta_a);

    double int_q_E1_rho, int_q_E1_phi, int_q_E1_z;
    double int_q_E2_rho, int_q_E2_phi, int_q_E2_z;
    double int_q_E_total_rho, int_q_E_total_phi, int_q_E_total_z;
    int_q(params, q, phi_a, theta_a, ra, phi_q, theta_q, rq,
        &int_q_E1_rho, &int_q_E1_phi, &int_q_E1_z,
        &int_q_E2_rho, &int_q_E2_phi, &int_q_E2_z,
        &int_q_E_total_rho, &int_q_E_total_phi, &int_q_E_total_z
        );

    *int_a_E1_rho      = k_a * int_q_E1_rho;
    *int_a_E1_phi      = k_a * int_q_E1_phi;
    *int_a_E1_z        = k_a * int_q_E1_z;

    *int_a_E2_rho      = k_a * int_q_E2_rho;
    *int_a_E2_phi      = k_a * int_q_E2_phi;
    *int_a_E2_z        = k_a * int_q_E2_z;

    *int_a_E_total_rho = k_a * int_q_E_total_rho;
    *int_a_E_total_phi = k_a * int_q_E_total_phi;
    *int_a_E_total_z   = k_a * int_q_E_total_z;
}

int Integrand(const int *ndim, const cubareal xx[],
  const int *ncomp, cubareal ff[], void *userdata) {
    ProblemParams *params = (ProblemParams *)userdata;

    // theta_a, ra, phi_q, theta_q, rq
#if 0
    #define psi_a   xx[0]
    #define theta_a xx[1]
    #define ra      xx[2]
    #define psi_q   xx[3]
    #define theta_q xx[4]
    #define rq      xx[5]
#else
    #define psi_a 0.0
    #define theta_a xx[0]
    #define ra      xx[1]
    #define psi_q   xx[2]
    #define theta_q xx[3]
    #define rq      xx[4]
#endif
    #define observed_ratio_1_rho ff[0]
    #define observed_ratio_1_phi ff[1]
    #define observed_ratio_1_z   ff[2]

    #define observed_ratio_2_rho ff[3]
    #define observed_ratio_2_phi ff[4]
    #define observed_ratio_2_z   ff[5]

    #define observed_ratio_rho ff[6]
    #define observed_ratio_phi ff[7]
    #define observed_ratio_z   ff[8]

    if ( fabs(theta_a - theta_q) < 1e-12 && fabs(ra - rq) < 1e-12 && fabs(psi_q - psi_a) < 1e-12)
    {
        printf ("theta_a = %e ", theta_a);
        printf ("theta_q = %e ", theta_q);
        printf ("ra = %e ", ra);
        printf ("rq = %e ", rq);
        printf ("psi_a = %e\n", psi_a);
        printf ("psi_q = %e\n", psi_q);

        return -999;
    }

    cubareal r0 = params->R0;
    cubareal q  = 1.0;

    /* Вычисление ускорения */
    double Gamma = -params->rho0 * pow(params->omega, 2);

    double k = r0 * r0 * (2 * M_PI) * M_PI * (2*M_PI) * M_PI;

    double int_a_E1_rho, int_a_E1_phi, int_a_E1_z;
    double int_a_E2_rho, int_a_E2_phi, int_a_E2_z;
    double int_a_E_rho,  int_a_E_phi,  int_a_E_z;
    int_a(params, q, psi_a * (2*M_PI), theta_a * M_PI, ra * r0, psi_q * (2*M_PI), theta_q * M_PI, rq * r0,
        &int_a_E1_rho, &int_a_E1_phi, &int_a_E1_z,
        &int_a_E2_rho, &int_a_E2_phi, &int_a_E2_z,
        &int_a_E_rho,  &int_a_E_phi,  &int_a_E_z
        );

    double f1_rho = k * int_a_E1_rho;
    double f1_phi = k * int_a_E1_phi;
    double f1_z   = k * int_a_E1_z;

    double f2_rho = k * int_a_E2_rho;
    double f2_phi = k * int_a_E2_phi;
    double f2_z   = k * int_a_E2_z;

    double f_rho = k * int_a_E_rho;
    double f_phi = k * int_a_E_phi;
    double f_z   = k * int_a_E_z;

    /*
        Твою энегрию электрического поля для сравнения с моим результатом самодействия я умножил на 2
        double U = 3.0 / (5.0 * params.R0);     // Энергия электрического поля 
        double m_perp_B = 2 * U / (params.c * params.c);  // Вариация типа В (теоретическое значение) 
        Потому что никто за сто с лишним лет на это не обращает внимание, но даже у Ферми есть ошибка в строке 
        ```
            или, заменив $O$ на $Р'$ (что ничего не изменяет) и взяв полусумму двух полученных таким образом значений,

            $$\frac{1}{2} \int \int\frac{\overrightarrow{P-P'}}{r^3} \, \frac{\vec {Г} \cdot \overrightarrow{\left(P - P'\right)}}{c^2} \, de\,de'.$$
        ```
        Ошибка в том, что при вычислении самодействия вычисляя например силу электромагнитной индукции с помощью которой кусочек одного и того же зарадя 1 действует на кусочек того же заряда 2 мы не должны брать полусумму, потому что сила взаимоиндукции ЭМ поля (вызванная ускорением кусочка заряда 2) которая действует на кусочек заряда 1 и точно так же сила взаимоиндукции ЭМ поля (вызванная ускорением кусочка заряда 1) которая действует на кусочек заряда 2 они обе должны суммироваться в общий интеграл самодействия.

        Что касается множителя 1/2 который вычислячется при вычислении электромагнитного поля например в параграфах 15, 16 у И.Е.Тамма то (мы можем позже отдельно проанализировать это) но мне кажется что для здесь этот множитель взят законно. А при вычислении силы самодействия он берётся незаконно. Если анализировать дальше следствия которые вытекают из этого факта, то нам придётся очень глубоко углубляться в теорию.

        Кстати косвенным указанием моей правоты насчёт этого коэффициента 1/2 является цитата из русской Википедии

        Русская википедия в статье "Классический радиус электрона"  на момент написания данной работы даёт следующее определение:

        ```
            Классический радиус электрона равен радиусу полой сферы, на которой равномерно распределён заряд, если этот заряд равен заряду электрона, а потенциальная энергия электростатического поля ${U}_{0}$ полностью эквивалентна половине массы электрона (без учета квантовых эффектов):

            ${\displaystyle U_{0}={\frac {1}{2}}{\frac {1}{4\pi \varepsilon _{0}}}\cdot {\frac {e^{2}}{r_{0}}}={\frac {1}{2}}m_{0}c^{2}}$.
        ```

        Возникает закономерный вопрос: а чему соответствует вторая половина массы электрона?

        Но об этом мы потом поговорим - это очень большой разговор, сейчас я просто обьясняю свои потивы почему я в код твоей программы добавляю множитель 2
    */

    double U = 3.0 / (5.0 * params->R0);     /* Энергия электрического поля */
    double m_perp_B = /*2 **/ U / (params->c * params->c);                 /* Вариация типа В (теоретическое значение) */

    /* Вычисление поперечной массы */
    //double m_perp_A = integral[0] / Gamma / (params.c * params.c);  /* Вариация типа А */
    //double m = f/Gamma/ (params->c * params->c);

    double m1_rho = f1_rho / Gamma;
    double m1_phi = f1_phi / Gamma;
    double m1_z   = f1_z   / Gamma;

    double m2_rho = f2_rho / Gamma;
    double m2_phi = f2_phi / Gamma;
    double m2_z   = f2_z   / Gamma;

    double m_rho = f_rho / Gamma;
    double m_phi = f_phi / Gamma;
    double m_z   = f_z   / Gamma;

    /* Проверка коэффициента 4/3 */
    double expected_ratio = 4.0 / 3.0;

    observed_ratio_1_rho = m1_rho / m_perp_B;
    observed_ratio_1_phi = m1_phi / m_perp_B;
    observed_ratio_1_z   = m1_z   / m_perp_B;

    observed_ratio_2_rho = m2_rho / m_perp_B;
    observed_ratio_2_phi = m2_phi / m_perp_B;
    observed_ratio_2_z   = m2_z   / m_perp_B;

    observed_ratio_rho = m_rho / m_perp_B;
    observed_ratio_phi = m_phi / m_perp_B;
    observed_ratio_z   = m_z   / m_perp_B;

    return 0;
}

/*********************************************************************/
#define NDIM 6
#define NCOMP 9
#define USERDATA NULL
#define NVEC 1
#if 1
#define EPSREL 1e-8
#define EPSABS 1e-16
#else
#define EPSREL 1e-3
#define EPSABS 1e-12
#endif
#define VERBOSE 2
#define LAST 4
#define SEED 0
#define MINEVAL 0
#define MAXEVAL 500000

#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

#define NNEW 1000
#define NMIN 2
#define FLATNESS 25.

#define KEY1 47
#define KEY2 1
#define KEY3 1
#define MAXPASS 5
#define BORDER 0.
#define MAXCHISQ 10.
#define MINDEVIATION .25
#define NGIVEN 0
#define LDXGIVEN NDIM
#define NEXTRA 0

#define KEY 0

int integrate(
    double R0,       /* Радиус сферы */
    double rho0,     /* Радиус орбиты */
    double omega,    /* Угловая скорость */
    double c,        /* Скорость света */
    int use_delay,
    int use_lorentz_factor,
    int use_lorentz_general_factor,
    int use_fermi_factor_O,
    int use_fermi_factor,
    int use_fermi_general_factor,
    int use_fast_integrand,
    cubareal* integral, cubareal* error, cubareal* prob)
{
    int comp, nregions, neval, fail;

    /* Параметры задачи */
    ProblemParams params;

    params.R0    = R0;    /* радиус сферы */
    params.rho0  = rho0;  /* радиус орбиты */
    params.omega = omega; /* угловая скорость */
    params.c    = c;      /* скорость света в м/с */

    params.use_delay                  = use_delay;
    params.use_lorentz_factor         = use_lorentz_factor;
    params.use_lorentz_general_factor = use_lorentz_general_factor;
    params.use_fermi_factor_O         = use_fermi_factor_O;
    params.use_fermi_factor           = use_fermi_factor;
    params.use_fermi_general_factor   = use_fermi_general_factor;
    params.use_fast_integrand         = use_fast_integrand;

    double v = params.rho0 * params.omega;
    double v_c =  v / params.c;  /* отношение скорости заряда к скорости света*/

#if 1
  printf("-------------------- Vegas test --------------------\n");

    Vegas(NDIM, NCOMP, Integrand, &params, NVEC,
        EPSREL, EPSABS, VERBOSE, SEED,
        MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
        GRIDNO, STATEFILE, SPIN,
        &neval, &fail, integral, error, prob);
#endif

#if 0
  printf("\n-------------------- Suave test --------------------\n");

  Suave(NDIM, NCOMP, Integrand, USERDATA, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST, SEED,
    MINEVAL, MAXEVAL, NNEW, NMIN, FLATNESS,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);
#endif

#if 0
  printf("\n------------------- Divonne test -------------------\n");

  Divonne(NDIM, NCOMP, Integrand, USERDATA, NVEC,
    EPSREL, EPSABS, VERBOSE, SEED,
    MINEVAL, MAXEVAL, KEY1, KEY2, KEY3, MAXPASS,
    BORDER, MAXCHISQ, MINDEVIATION,
    NGIVEN, LDXGIVEN, NULL, NEXTRA, NULL,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);
#endif

#if 0
  printf("\n-------------------- Cuhre test --------------------\n");

  Cuhre(NDIM, NCOMP, Integrand, USERDATA, NVEC,
    EPSREL, EPSABS, VERBOSE | LAST,
    MINEVAL, MAXEVAL, KEY,
    STATEFILE, SPIN,
    &nregions, &neval, &fail, integral, error, prob);
#endif

}
