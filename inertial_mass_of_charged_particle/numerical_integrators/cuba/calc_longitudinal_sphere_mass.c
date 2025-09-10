/*
	calc_longitudinal_sphere_mass.c
*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef INCLUDE_CUBA_H
#include "cuba.h"
#endif

static inline cubareal Sq(cubareal x) {
  return x*x;
}

static inline cubareal Cb(cubareal x) {
  return x*x*x;
}

/* Структура для передачи параметров задачи */
typedef struct {
    double R0;    /* Радиус сферы */
    double v0;    /* Начальная продольная скорость */
    double a;     /* Продольное ускорение */
    double c;     /* Скорость света */
    double t;     /* Текущее время */
    double t0;    /* Начальное время */
    int use_delay;
    int use_lorentz_factor;
    int use_lorentz_general_factor;
    int use_fermi_factor;
    int use_fermi_general_factor;
    int use_fast_integrand;
} ProblemParams;

/* Функция для вычисления поправки Ферми */
static double fermi_correction(
    double R_x,
    double R_y,
    double R_z,
    double a_x,
    double a_y,
    double a_z,
    double c
) {
    /* Скалярное произведение R · a */
    double R_dot_a = R_x * a_x + R_y * a_y + R_z * a_z;

    /* Поправка Ферми */
    return (/*1.0*/ + 0.5 * R_dot_a / (c * c));
}

/* Функция для вычисления поправки Ферми в общем релятивистском случае */
static double fermi_correction_general(
    double t,           /* текущее время */
    double t_prime,     /* запаздывающее время */
    double R_x,
    double R_y,
    double R_z,
    double a_x,
    double a_y,
    double a_z,
    double v_x,
    double v_y,
    double v_z,
    double c
) {
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
    double t,           /* текущее время */
    double t_prime,     /* запаздывающее время */
    double R,
    double R_x, double R_y, double R_z,
    double v_x, double v_y, double v_z,
    double a_x, double a_y, double a_z,
    double c,
    double *E1_x, double *E1_y, double *E1_z,
    double *E2_x, double *E2_y, double *E2_z,
    double *E_total_x, double *E_total_y, double *E_total_z,
    int use_fermi_factor,
    int use_fermi_general_factor,
    int use_fast_integrand
)
{
    /* Радиус Лиенара-Вихерта (учет запаздывания) */
    double R_star = R - (R_z * v_z) / c;

    /* Защита от деления на ноль */
    if (R_star < 1e-10 || R < 1e-10) {
        *E1_x = *E1_y = *E1_z = 0.0;
        *E2_x = *E2_y = *E2_z = 0.0;
        *E_total_x = *E_total_y = *E_total_z = 0.0;
        return;
    }

    double fermi_factor = 0.0;

    if (use_fermi_general_factor)
    {
        fermi_factor = fermi_correction_general(
            t,           /* текущее время */
            t_prime,     /* запаздывающее время */
            R_x,
            R_y,
            R_z,
            a_x,
            a_y,
            a_z,
            v_x,
            v_y,
            v_z,
            c);
    }
    else if (use_fermi_factor)
    {
        fermi_factor = fermi_correction(R_x, R_y, R_z, a_x, a_y, a_z, c);
    }
    double v2_c2 = (Sq(v_x) + Sq(v_y) + Sq(v_z) ) / Sq(c);
    double Ra_c2 = (R_z * a_z) / Sq(c);

    /* Вычисление градиентного поля E1 */
    double common_factor1 = 1.0 / pow(R_star, 2);
    /* Вычисление градиентного поля E1 (только слагаемое с ускорением) */
    /* (1.0 + Ra_c2 - v2_c2) * (1.0+fermi_factor)  = 1.0 + Ra_v2_c2_fermi */
    double Ra_v2_c2_fermi = fermi_factor + Ra_c2 + fermi_factor * (Ra_c2 - v2_c2);

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

    *E1_x = common_factor1 * (R_x * (f_one + Ra_v2_c2_fermi) / R_star - f_one * v_x / c - v_x / c * fermi_factor);
    *E1_y = common_factor1 * (R_y * (f_one + Ra_v2_c2_fermi) / R_star - f_one * v_y / c - v_y / c * fermi_factor);
    *E1_z = common_factor1 * (R_z * (f_one + Ra_v2_c2_fermi) / R_star - f_one * v_z / c - v_z / c * fermi_factor);

    /* Вычисление поля самоиндукции E2 */
    double common_factor2 = 1.0 / pow(R_star, 2);
    double velocity_factor2 = (R / R_star) * (v2_c2 - Ra_c2 - 1.0) + 1.0;

    common_factor2 *= (1.0+fermi_factor);

    *E2_x = common_factor2 * (v_x / c * velocity_factor2 - a_x * R / Sq(c));
    *E2_y = common_factor2 * (v_y / c * velocity_factor2 - a_y * R / Sq(c));
    *E2_z = common_factor2 * (v_z / c * velocity_factor2 - a_z * R / Sq(c));

    /* Вычисление суммарного поля по полной формуле Лиенара-Вихерта */
    double common_factor = 1.0 / Cb(R_star);

    //E_x = common_factor * (1.0 + fermi_factor) * ((1.0 + Ra_v2_c2) * (R_x - R*v_x/c)  - (a_x * R_star * R) / Sq(c));
    //E_x = common_factor * ((1.0 + Ra_v2_c2) * (1.0 + fermi_factor) * (R_x - R*v_x/c)  - (a_x * R_star * R) * (1.0+fermi_factor) / Sq(c));
    //E_x = common_factor * ((1.0 + Ra_v2_c2 + fermi_factor + fermi_factor * Ra_v2_c2) * (R_x - R*v_x/c)  - (a_x * R_star * R) / Sq(c) - (a_x * R_star * R) * fermi_factor / Sq(c));
    //E_x = common_factor * ((1.0 + Ra_v2_c2_fermi) * (R_x - R*v_x/c)  - (a_x * R_star * R) / Sq(c) - (a_x * R_star * R) * fermi_factor / Sq(c));
    *E_total_x = common_factor * ( ((f_one + Ra_v2_c2_fermi) * R_x - (f_one + Ra_v2_c2_fermi) * R*v_x/c)  - (a_x * R_star * R) / Sq(c) - (a_x * R_star * R) * fermi_factor / Sq(c));
    *E_total_y = common_factor * ( ((f_one + Ra_v2_c2_fermi) * R_y - (f_one + Ra_v2_c2_fermi) * R*v_y/c)  - (a_y * R_star * R) / Sq(c) - (a_y * R_star * R) * fermi_factor / Sq(c));
    *E_total_z = common_factor * ( ((f_one + Ra_v2_c2_fermi) * R_z - (f_one + Ra_v2_c2_fermi) * R*v_z/c)  - (a_z * R_star * R) / Sq(c) - (a_z * R_star * R) * fermi_factor / Sq(c));
}

/* Функция для вычисления фактора Лоренца */
static inline double lorentz_factor(double v, double c) {
    double beta = v / c;
    return 1.0 / sqrt(1.0 - beta * beta);
}

/* Функция для вычисления электрического поля по Лиенару-Вихерту для продольного случая */
static void compute_electric_field(
    double ra, double theta_a, double phi_a, /*воображаемые координаты несплюснутой сферы*/
    double rq, double theta_q, double phi_q, /*по которым производится интегрирование*/
    const ProblemParams *params,
    double *E1_x, double *E1_y, double *E1_z,
    double *E2_x, double *E2_y, double *E2_z,
    double *E_total_x, double *E_total_y, double *E_total_z
) {
    /* Преобразование в декартовы координаты */
    double x_a = ra * sin(theta_a) * cos(phi_a);
    double y_a = ra * sin(theta_a) * sin(phi_a);
    double z_a = ra * cos(theta_a);

    double x_q = rq * sin(theta_q) * cos(phi_q);
    double y_q = rq * sin(theta_q) * sin(phi_q);
    double z_q = rq * cos(theta_q);

    /* Скорость и ускорение источника */
    double v_x = 0;
    double v_y = 0;
    // double v_z = params->v;  /* Продольная скорость */
    double v_z = params->v0 + params->a * (params->t - params->t0);

    double a_x = 0;
    double a_y = 0;
    double a_z = params->a;  /* Продольное ускорение */

    if (params->use_lorentz_factor || params->use_lorentz_general_factor)
    {
        /* Учет Лоренц-сокращения в направлении движения */
        double gamma = lorentz_factor(v_z, params->c);
        /* Учет Лоренц-сокращения: в системе отсчета наблюдателя
        сфера сжата вдоль направления движения (ось z) */
        z_q = z_q / gamma;
        z_a = z_a / gamma;
    }

    /* Вектор от источника к наблюдателю */
    double R_x = x_a - x_q;
    double R_y = y_a - y_q;
    double R_z = z_a - z_q;
    double R = sqrt(R_x*R_x + R_y*R_y + R_z*R_z);

    computing_electric_field(
        0,           /* текущее время */
        0,           /* запаздывающее время */
        R,
        R_x, R_y, R_z,
        v_x, v_y, v_z,
        a_x, a_y, a_z,
        params->c,
        E1_x, E1_y, E1_z,
        E2_x, E2_y, E2_z,
        E_total_x, E_total_y, E_total_z,
        params->use_fermi_factor,
        params->use_fermi_general_factor,
        params->use_fast_integrand
    );
}


/* Функция для вычисления собственного времени */
static double compute_proper_time(
    double t,
    double t0,
    double v0,
    double a,
    double c
) {
    /* Для постоянного ускорения собственное время вычисляется как */
    double v = v0 + a * (t - t0);
    double beta = v / c;
    double gamma = lorentz_factor(v, c);

    /* Интеграл от 1/gamma(t) по координатному времени */
    double proper_time = 0.0;
    if (fabs(a) > 1e-10) {
        /* Для постоянного ускорения */
        proper_time = (c / a) * log((a * (t - t0) + v0 + sqrt(pow(c, 2) + pow(v0, 2))) / 
                          (v0 + sqrt(pow(c, 2) + pow(v0, 2))));
    } else {
        /* Для случая a = 0 (постоянная скорость) */
        proper_time = (t - t0) / gamma;
    }

    return proper_time;
}

/* Функция для вычисления координат источника в запаздывающий момент времени */
/* Функция для вычисления координат источника с учетом Лоренц-сокращения */
static void get_source_position(
    double t_prime,
    double t0,
    double rq, double theta_q, double phi_q, /*воображаемые координаты несплюснутой сферы*/
    double v0,
    double a,
    double c,
    double *x_q,
    double *y_q,
    double *z_q,
    int use_lorentz_factor,
    int use_lorentz_general_factor
) {
    /* Вычисляем мгновенную скорость в запаздывающий момент времени */
    double v = v0 + a * (t_prime - t0);

    /* Положение центра сферы в запаздывающий момент времени */
    double z_center = v0 * (t_prime - t0) + 0.5 * a * pow(t_prime - t0, 2);

    /* Координаты точки на сфере в системе покоя сферы */
    double z_q0 = rq * cos(theta_q);

    *x_q = rq * sin(theta_q) * cos(phi_q);
    *y_q = rq * sin(theta_q) * sin(phi_q);
    if (use_lorentz_factor || use_lorentz_general_factor)
    {
        /* Учет Лоренц-сокращения в направлении движения */
        double gamma = lorentz_factor(v, c);
        /* Учет Лоренц-сокращения: в системе отсчета наблюдателя
        сфера сжата вдоль направления движения (ось z) */
        double z_q_rel = z_q0 / gamma;
        *z_q = z_center + z_q_rel;

        if (use_lorentz_general_factor) {
            /* Учет дополнительного искривления из-за ускорения */
            double proper_time = compute_proper_time(t_prime, t0, v0, a, c);
            double acceleration_correction = 0.5 * a * pow(proper_time, 2) / gamma;
            *z_q += acceleration_correction;
        }
    }
    else
    {
        /* В продольном случае движение происходит только вдоль оси z */
        *z_q = z_center + z_q0;
    }
}

/* Функция для вычисления R(t') - расстояния от источника к наблюдателю в момент t' */
static double compute_R(
    double t_prime,
    double t,
    double t0,
    double ra,
    double theta_a,
    double phi_a,
    double rq,
    double theta_q,
    double phi_q,
    double v0,
    double a,
    double c,
    int use_lorentz_factor,
    int use_lorentz_general_factor
) {
    /* Координаты наблюдателя в момент t */
    double x_a, y_a, z_a;
    get_source_position(t, t0, ra, theta_a, phi_a, v0, a, c, &x_a, &y_a, &z_a,
        use_lorentz_factor,
        use_lorentz_general_factor);

    /* Координаты источника в момент t' */
    double x_q, y_q, z_q;
    get_source_position(t_prime, t0, rq, theta_q, phi_q, v0, a, c, &x_q, &y_q, &z_q,
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
    double t,
    double t0,
    double ra,
    double theta_a,
    double phi_a,
    double rq,
    double theta_q,
    double phi_q,
    double v0,
    double a,
    double c,
    int use_lorentz_factor,
    int use_lorentz_general_factor
) {
    /* Координаты наблюдателя в момент t */
    double x_a, y_a, z_a;
    get_source_position(t, t0, ra, theta_a, phi_a, v0, a, c, &x_a, &y_a, &z_a,
        use_lorentz_factor,
        use_lorentz_general_factor);

    /* Координаты источника в запаздывающий момент времени t_prime */
    double x_q, y_q, z_q;
    get_source_position(t_prime, t0, rq, theta_q, phi_q, v0, a, c, &x_q, &y_q, &z_q,
        use_lorentz_factor,
        use_lorentz_general_factor);

    /* Вектор от источника к наблюдателю */
    double R_x = x_a - x_q;
    double R_y = y_a - y_q;
    double R_z = z_a - z_q;
    double R = sqrt(R_x*R_x + R_y*R_y + R_z*R_z);

    /* Скорость источника в запаздывающий момент времени */
    double v = v0 + a * (t_prime - t0);
    double gamma = lorentz_factor(v, c);
    double v_r = v / gamma;  /* Радиальная скорость с учетом Лоренц-сокращения */

    /* Вычисляем проекцию скорости на направление к наблюдателю */
    double cos_theta = R_z / R;  /* Для продольного случая */
    double v_radial = v_r * cos_theta;

    /* Аналитическое выражение для dR/dt' */
    return -v_radial;
}

/* Улучшенное вычисление запаздывающего времени с методом Ньютона */
static double compute_retarded_time_newton(
    double t,
    double t0,
    double ra,
    double theta_a,
    double phi_a,
    double rq,
    double theta_q,
    double phi_q,
    double v0,
    double a,
    double c,
    int max_iter,
    double tolerance,
    int use_lorentz_factor,
    int use_lorentz_general_factor
) {
    /* Координаты наблюдателя в момент t */
    double x_a, y_a, z_a;
    get_source_position(t, t0, ra, theta_a, phi_a, v0, a, c, &x_a, &y_a, &z_a,
        use_lorentz_factor,
        use_lorentz_general_factor);

    /* Координаты источника в момент времени t */
    double x_q, y_q, z_q;
    get_source_position(t, t0, rq, theta_q, phi_q, v0, a, c, &x_q, &y_q, &z_q,
        use_lorentz_factor, use_lorentz_general_factor);

    double R = sqrt(pow(x_a - x_q, 2) + pow(y_a - y_q, 2) + pow(z_a - z_q, 2));
    double t_prime = t - R / c;

    /* Итерационное решение */
    for (int i = 0; i < max_iter; i++) {
        get_source_position(t_prime, t0, rq, theta_q, phi_q, v0, a, c, &x_q, &y_q, &z_q,
            use_lorentz_factor, use_lorentz_general_factor);

        R = sqrt(pow(x_a - x_q, 2) + pow(y_a - y_q, 2) + pow(z_a - z_q, 2));
        double t_prime_new = t - R / c;

        if (fabs(t_prime_new - t_prime) < tolerance) {
            t_prime = t_prime_new;
            break;
        }

        t_prime = t_prime_new;
    }

    /* Уточнение методом Ньютона */
    /* Уточнение методом Ньютона с аналитической производной */
    for (int i = 0; i < 5; i++) {
        double R_val = compute_R(t_prime, t, t0,
            ra, theta_a, phi_a,
            rq, theta_q, phi_q,
            v0, a, c,
            use_lorentz_factor, use_lorentz_general_factor);
        double f = t_prime - t + R_val / c;
        double df_dt = 1.0 + compute_dR_dt_prime_analytical(t_prime, t, t0,
            ra, theta_a, phi_a,
            rq, theta_q, phi_q,
            v0, a, c,
            use_lorentz_factor, use_lorentz_general_factor) / c;

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
    double ra, double theta_a, double phi_a, /*воображаемые координаты несплюснутой сферы*/
    double rq, double theta_q, double phi_q,
    const ProblemParams *params,
    double *E1_x, double *E1_y, double *E1_z,
    double *E2_x, double *E2_y, double *E2_z,
    double *E_total_x,
    double *E_total_y,
    double *E_total_z
) {
    /* Вычисление запаздывающего времени */
    double t_prime = compute_retarded_time_newton(
        params->t, params->t0, ra, theta_a, phi_a, rq, theta_q, phi_q,
        params->v0, params->a, params->c, 100, 1e-10,
        params->use_lorentz_factor, params->use_lorentz_general_factor
    );

    /* Координаты наблюдателя в момент t */
    double x_a, y_a, z_a;
    get_source_position(params->t, params->t0, ra, theta_a, phi_a, params->v0, params->a, params->c, &x_a, &y_a, &z_a,
        params->use_lorentz_factor,
        params->use_lorentz_general_factor);

    /* Получение координат источника в запаздывающий момент времени */
    double x_q, y_q, z_q;
    get_source_position(t_prime, params->t0, rq, theta_q, phi_q,
                       params->v0, params->a, params->c, &x_q, &y_q, &z_q,
                       params->use_lorentz_factor, params->use_lorentz_general_factor);

    /* Вектор от источника к наблюдателю */
    double R_x = x_a - x_q;
    double R_y = y_a - y_q;
    double R_z = z_a - z_q;
    double R = sqrt(R_x*R_x + R_y*R_y + R_z*R_z);

    /* Скорость источника в запаздывающий момент времени */
    double v_z = params->v0 + params->a * (t_prime - params->t0);
    double v_x = 0.0;
    double v_y = 0.0;

    /* Ускорение источника (предполагаем постоянное) */
    double a_x = 0.0;
    double a_y = 0.0;
    double a_z = params->a;

    computing_electric_field(
        params->t,           /* текущее время */
        t_prime,     /* запаздывающее время */
        R,
        R_x, R_y, R_z,
        v_x, v_y, v_z,
        a_x, a_y, a_z,
        params->c,
        E1_x, E1_y, E1_z,
        E2_x, E2_y, E2_z,
        E_total_x, E_total_y, E_total_z,
        params->use_fermi_factor,
        params->use_lorentz_general_factor,
        params->use_fast_integrand);
}

// the charge distribution
static inline cubareal rho_q (cubareal r0, cubareal q)
{
    return 3 * (q) / (4*M_PI*(r0)*(r0)*(r0));
}

/* Функция для вычисления плотности заряда с учетом Лоренц-сокращения */
/*static inline double rho_q_lorentz(double r0, double q, double v, double c) {
    double gamma = lorentz_factor(v, c);
    return gamma * 3.0 * q / (4.0 * M_PI * pow(r0, 3));
}*/

// интегрирование по координатам заряда источника потенциала
static inline void int_q (const ProblemParams *params, cubareal q, cubareal psi_a, cubareal theta_a, cubareal ra, cubareal psi_q, cubareal theta_q, cubareal rq,
    double *int_q_E1_x, double *int_q_E1_y, double *int_q_E1_z,
    double *int_q_E2_x, double *int_q_E2_y, double *int_q_E2_z,
    double *int_q_E_x,  double *int_q_E_y,  double *int_q_E_z

 )
{
    cubareal r0 = params->R0;

    /* Вычисление электрического поля */
    double E1_x, E1_y, E1_z;
    double E2_x, E2_y, E2_z;
    double E_x,  E_y,  E_z;

    if (params->use_delay)
    {
        compute_electric_field_with_delay(ra, theta_a, psi_a,
            rq, theta_q, psi_q, params,
            &E1_x, &E1_y, &E1_z,
            &E2_x, &E2_y, &E2_z,
            &E_x,  &E_y,  &E_z);
    }
    else {
        compute_electric_field(ra, theta_a, psi_a,
            rq, theta_q, psi_q, params,
            &E1_x, &E1_y, &E1_z,
            &E2_x, &E2_y, &E2_z,
            &E_x,  &E_y,  &E_z);
    }

    double rho;
    // Если вы интегрируете по "воображаемой" сфере (сферическим координатам без учета сжатия):
    // Для получения правильного результата вам нужна плотность rho​ , а не gamma * rho
    /*if (params->use_lorentz_factor || params->use_lorentz_general_factor)
    {
        double v_z = params->v0 + params->a * (params->t - params->t0);
        rho = rho_q_lorentz(r0, q, v_z, params->c);
    }
    else*/
    {
        rho = rho_q(r0, q);
    }

    double k_q = rho * Sq(rq) * sin(theta_q);

    *int_q_E1_x = k_q * E1_x;
    *int_q_E1_y = k_q * E1_y;
    *int_q_E1_z = k_q * E1_z;
    *int_q_E2_x = k_q * E2_x;
    *int_q_E2_y = k_q * E2_y;
    *int_q_E2_z = k_q * E2_z;
    *int_q_E_x  = k_q * E_x;
    *int_q_E_y  = k_q * E_y;
    *int_q_E_z  = k_q * E_z;
}

// интегрирование по координатам точек наблюдения
void int_a (
    const ProblemParams *params, cubareal q,
    cubareal phi_a, cubareal theta_a, cubareal ra,
    cubareal phi_q, cubareal theta_q, cubareal rq,
    double *int_a_E1_x, double *int_a_E1_y, double *int_a_E1_z,
    double *int_a_E2_x, double *int_a_E2_y, double *int_a_E2_z,
    double *int_a_E_total_x, double *int_a_E_total_y, double *int_a_E_total_z
)
{
    cubareal r0 = params->R0;

    double rho;
    /*if (params->use_lorentz_factor || params->use_lorentz_general_factor)
    {
        double v_z = params->v0 + params->a * (params->t - params->t0);
        rho = rho_q_lorentz(r0, q, v_z, params->c);
    }
    else*/
    {
        rho = rho_q(r0, q);
    }

    double k_a = rho * Sq(ra) * sin(theta_a);

    double int_q_E1_x, int_q_E1_y, int_q_E1_z;
    double int_q_E2_x, int_q_E2_y, int_q_E2_z;
    double int_q_E_total_x, int_q_E_total_y, int_q_E_total_z;
    int_q(params, q, phi_a, theta_a, ra, phi_q, theta_q, rq,
        &int_q_E1_x, &int_q_E1_y, &int_q_E1_z,
        &int_q_E2_x, &int_q_E2_y, &int_q_E2_z,
        &int_q_E_total_x, &int_q_E_total_y, &int_q_E_total_z
        );

    *int_a_E1_x      = k_a * int_q_E1_x;
    *int_a_E1_y      = k_a * int_q_E1_y;
    *int_a_E1_z      = k_a * int_q_E1_z;

    *int_a_E2_x      = k_a * int_q_E2_x;
    *int_a_E2_y      = k_a * int_q_E2_y;
    *int_a_E2_z      = k_a * int_q_E2_z;

    *int_a_E_total_x = k_a * int_q_E_total_x;
    *int_a_E_total_y = k_a * int_q_E_total_y;
    *int_a_E_total_z = k_a * int_q_E_total_z;
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
    #define observed_ratio_1_x ff[0]
    #define observed_ratio_1_y ff[1]
    #define observed_ratio_1_z   ff[2]

    #define observed_ratio_2_x ff[3]
    #define observed_ratio_2_y ff[4]
    #define observed_ratio_2_z   ff[5]

    #define observed_ratio_x ff[6]
    #define observed_ratio_y ff[7]
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
    double Gamma = params->a;

    double k = r0 * r0 * (2 * M_PI) * M_PI * (2*M_PI) * M_PI;

    double int_a_E1_x, int_a_E1_y, int_a_E1_z;
    double int_a_E2_x, int_a_E2_y, int_a_E2_z;
    double int_a_E_x,  int_a_E_y,  int_a_E_z;
    int_a(params, q, psi_a * (2*M_PI), theta_a * M_PI, ra * r0, psi_q * (2*M_PI), theta_q * M_PI, rq * r0,
        &int_a_E1_x, &int_a_E1_y, &int_a_E1_z,
        &int_a_E2_x, &int_a_E2_y, &int_a_E2_z,
        &int_a_E_x,  &int_a_E_y,  &int_a_E_z
        );

    double f1_x = k * int_a_E1_x;
    double f1_y = k * int_a_E1_y;
    double f1_z = k * int_a_E1_z;

    double f2_x = k * int_a_E2_x;
    double f2_y = k * int_a_E2_y;
    double f2_z = k * int_a_E2_z;

    double f_x = k * int_a_E_x;
    double f_y = k * int_a_E_y;
    double f_z = k * int_a_E_z;

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

    /* Вычисление продольной массы */
    //double m_perp_A = integral[0] / Gamma / (params.c * params.c);  /* Вариация типа А */
    //double m = f/Gamma/ (params->c * params->c);

    double m1_x = f1_x / Gamma ;
    double m1_y = f1_y / Gamma ;
    double m1_z = f1_z / Gamma ;

    double m2_x = f2_x / Gamma ;
    double m2_y = f2_y / Gamma ;
    double m2_z = f2_z / Gamma ;

    double m_x = f_x / Gamma ;
    double m_y = f_y / Gamma ;
    double m_z = f_z / Gamma ;

    /* Проверка коэффициента 4/3 */
    double expected_ratio = 4.0 / 3.0;

    observed_ratio_1_x = m1_x / m_perp_B;
    observed_ratio_1_y = m1_y / m_perp_B;
    observed_ratio_1_z = m1_z / m_perp_B;

    observed_ratio_2_x = m2_x / m_perp_B;
    observed_ratio_2_y = m2_y / m_perp_B;
    observed_ratio_2_z = m2_z / m_perp_B;

    observed_ratio_x = m_x / m_perp_B;
    observed_ratio_y = m_y / m_perp_B;
    observed_ratio_z = m_z / m_perp_B;

    return 0;
}


/*********************************************************************/
#if 0
#define NDIM 6
#else
#define NDIM 5
#endif
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
#if 0
#define MAXEVAL 500000000
#else
#define MAXEVAL 500000
#endif

#if 0
#define NSTART 5000
#define NINCREASE 1000
#else
#define NSTART 1000
#define NINCREASE 500
#endif
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

int integrate(
    double R0,    /* Радиус сферы */
    double v0,    /* Начальная продольная скорость */
    double a,     /* Продольное ускорение */
    double c,     /* Скорость света */
    double t,     /* Текущее время */
    double t0,    /* Начальное время */
    int use_delay,
    int use_lorentz_factor,
    int use_lorentz_general_factor,
    int use_fermi_factor,
    int use_fermi_general_factor,
    int use_fast_integrand,
    cubareal* integral, cubareal* error, cubareal* prob)
{
    int comp, nregions, neval, fail;

    /* Параметры задачи */
    ProblemParams params;

    params.R0 = R0; /* радиус сферы (половина R1) */
    params.v0 = v0;
    params.t0 = t0;
    params.t = t;

    params.c = c;   /* скорость света в м/с */
    params.a = a;   /* продольное ускорение*/

    params.use_delay                  = use_delay;
    params.use_lorentz_factor         = use_lorentz_factor;
    params.use_lorentz_general_factor = use_lorentz_general_factor;
    params.use_fermi_factor           = use_fermi_factor;
    params.use_fermi_general_factor   = use_fermi_general_factor;
    params.use_fast_integrand         = use_fast_integrand;

    double v_z = params.v0 + params.a * (params.t - params.t0);
    double v_c =  v_z / params.c;  /* отношение скорости заряда к скорости света*/

    Vegas(NDIM, NCOMP, Integrand, &params, NVEC,
        EPSREL, EPSABS, VERBOSE, SEED,
        MINEVAL, MAXEVAL, NSTART, NINCREASE, NBATCH,
        GRIDNO, STATEFILE, SPIN,
        &neval, &fail, integral, error, prob);
}
