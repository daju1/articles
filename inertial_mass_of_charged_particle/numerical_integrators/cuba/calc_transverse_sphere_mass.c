/*
	demo-c.c
		test program for the Cuba library
		last modified 13 Mar 15 th
*/

#include "calc_RO.h"


#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#ifdef INCLUDE_CUBA_H
#include "cuba.h"
#endif


/* Структура для передачи параметров задачи */
typedef struct {
    double R0;     /* Радиус сферы */
    double rho0;   /* Расстояние от оси вращения до центра сферы */
    double omega;  /* Угловая скорость вращения */
    double c;      /* Скорость света */
} ProblemParams;

/* Функция для вычисления электрического поля по Лиенару-Вихерту */
static void compute_electric_field(
    double ra, double theta_a, double psi_a,
    double rq, double theta_q, double psi_q,
    const ProblemParams *params,
    double *E1_rho, double *E1_phi, double *E1_z,
    double *E2_rho, double *E2_phi, double *E2_z,
    double *E_total_rho, double *E_total_phi, double *E_total_z
) {
    /* Преобразование в декартовы координаты в локальной системе */
    double x_a = ra * sin(theta_a) * cos(psi_a);
    double y_a = ra * sin(theta_a) * sin(psi_a);
    double z_a = ra * cos(theta_a);

    double x_q = rq * sin(theta_q) * cos(psi_q);
    double y_q = rq * sin(theta_q) * sin(psi_q);
    double z_q = rq * cos(theta_q);

    /* Преобразование в цилиндрическую систему (исходная) */
    double rho_a = sqrt(pow(params->rho0 + x_a, 2) + pow(y_a, 2));
    double phi_a = atan2(y_a, params->rho0 + x_a);
    if (phi_a < 0) phi_a += 2.0 * M_PI;

    double rho_q = sqrt(pow(params->rho0 + x_q, 2) + pow(y_q, 2));
    double phi_q = atan2(y_q, params->rho0 + x_q);
    if (phi_q < 0) phi_q += 2.0 * M_PI;

    /* Вектор от источника к наблюдателю */
    double R_rho = rho_a * cos(phi_a - phi_q) - rho_q;
    double R_phi = rho_a * sin(phi_a - phi_q);
    double R_z = z_a - z_q;
    double R = sqrt(pow(R_rho, 2) + pow(R_phi, 2) + pow(R_z, 2));

    /* Скорость и ускорение источника */
    double v_q = rho_q * params->omega;
    double a_q = -rho_q * pow(params->omega, 2);

    /* Радиус Лиенара-Вихерта */
    double R_star = R - (R_phi * v_q) / params->c;

    /* Защита от деления на ноль */
    if (R_star < 1e-10) {
        *E1_rho = *E1_phi = *E1_z = 0.0;
        *E2_rho = *E2_phi = *E2_z = 0.0;
        *E_total_rho = *E_total_phi = *E_total_z = 0.0;
        return;
    }

    /* Вычисление градиентного поля E1 */
    double common_factor1 = 1.0 / pow(R_star, 2);
    double velocity_factor1 = /*1.0 +*/ (R_rho * a_q) / pow(params->c, 2) /*- pow(v_q / params->c, 2)*/;

    *E1_rho = common_factor1 * (R_rho * velocity_factor1 / R_star);
    *E1_phi = common_factor1 * (R_phi * velocity_factor1 / R_star - v_q / params->c);
    *E1_z   = common_factor1 * (R_z   * velocity_factor1 / R_star);

    /* Вычисление поля самоиндукции E2 */
    double common_factor2 = 1.0 / pow(R_star, 2);
    double velocity_factor2 = (R / R_star) * (pow(v_q / params->c, 2) - (R_rho * a_q) / pow(params->c, 2) - 1.0) + 1.0;

    *E2_rho = common_factor2 * (- a_q * R / pow(params->c, 2));
    *E2_phi = common_factor2 * (v_q / params->c * velocity_factor2);
    *E2_z = 0;

    /* Вычисление суммарного поля */
    *E_total_rho = *E1_rho + *E2_rho;
    *E_total_phi = *E1_phi + *E2_phi;
    *E_total_z = *E1_z + *E2_z;

    /* Основная формула Лиенара-Вихерта */
    double common_factor = 1.0 / pow(R_star, 3);
    double velocity_factor = 1.0 + (R_rho * a_q) / pow(params->c, 2) - pow(v_q, 2) / pow(params->c, 2);

    *E_total_rho = common_factor * (R_rho * velocity_factor - (a_q * R_star * R) / pow(params->c, 2));
    *E_total_phi = common_factor * ((R_phi - (R * v_q) / params->c) * velocity_factor);
    *E_total_z = common_factor * (R_z * velocity_factor);
}

/* Функция для вычисления электрического поля по Лиенару-Вихерту */
static void compute_electric_field_z_rho(
    double ra, double theta_a, double psi_a,
    double rq, double theta_q, double psi_q,
    const ProblemParams *params,
    double *E1_rho, double *E1_phi, double *E1_z,
    double *E2_rho, double *E2_phi, double *E2_z,
    double *E_total_rho, double *E_total_phi, double *E_total_z
) {
    /* Преобразование в декартовы координаты в локальной системе */
    /* Ось z сферической системы направлена вдоль радиального направления */
    double x_a = params->rho0 + ra * cos(theta_a);
    double y_a = ra * sin(theta_a) * sin(psi_a);
    double z_a = ra * sin(theta_a) * cos(psi_a);

    double x_q = params->rho0 + rq * cos(theta_q);
    double y_q = rq * sin(theta_q) * sin(psi_q);
    double z_q = rq * sin(theta_q) * cos(psi_q);

    /* Преобразование в цилиндрическую систему (исходная) */
    double rho_a = sqrt(pow(x_a, 2) + pow(y_a, 2));
    double phi_a = atan2(y_a, x_a);
    if (phi_a < 0) phi_a += 2.0 * M_PI;

    double rho_q = sqrt(pow(x_q, 2) + pow(y_q, 2));
    double phi_q = atan2(y_q, x_q);
    if (phi_q < 0) phi_q += 2.0 * M_PI;

    /* Вектор от источника к наблюдателю */
    double R_rho = rho_a * cos(phi_a - phi_q) - rho_q;
    double R_phi = rho_a * sin(phi_a - phi_q);
    double R_z = z_a - z_q;
    double R = sqrt(pow(R_rho, 2) + pow(R_phi, 2) + pow(R_z, 2));

    /* Скорость и ускорение источника */
    double v_q = rho_q * params->omega;
    double a_q = -rho_q * pow(params->omega, 2);

    /* Радиус Лиенара-Вихерта */
    double R_star = R - (R_phi * v_q) / params->c;

    /* Защита от деления на ноль */
    if (R_star < 1e-10 || R < 1e-10) {
        *E1_rho = *E1_phi = *E1_z = 0.0;
        *E2_rho = *E2_phi = *E2_z = 0.0;
        *E_total_rho = *E_total_phi = *E_total_z = 0.0;
        return;
    }

    /* Вычисление градиентного поля E1 */
    double common_factor1 = 1.0 / pow(R_star, 2);
    double velocity_factor1 = /*1.0 +*/ (R_rho * a_q) / pow(params->c, 2) /*- pow(v_q / params->c, 2)*/;

    *E1_rho = common_factor1 * (R_rho * velocity_factor1 / R_star);
    *E1_phi = common_factor1 * (R_phi * velocity_factor1 / R_star - v_q / params->c);
    *E1_z   = common_factor1 * (R_z   * velocity_factor1 / R_star);

    /* Вычисление поля самоиндукции E2 */
    double common_factor2 = 1.0 / pow(R_star, 2);
    double velocity_factor2 = (R / R_star) * (pow(v_q / params->c, 2) - (R_rho * a_q) / pow(params->c, 2) - 1.0) + 1.0;

    *E2_rho = common_factor2 * (- a_q * R / pow(params->c, 2));
    *E2_phi = common_factor2 * (v_q / params->c * velocity_factor2);
    *E2_z = 0;

    /* Вычисление суммарного поля */
    double common_factor = 1.0 / pow(R_star, 3);
    double velocity_factor = /*1.0 +*/ (R_rho * a_q) / pow(params->c, 2) /*- pow(v_q, 2) / pow(params->c, 2)*/;

    *E_total_rho = common_factor * (R_rho * velocity_factor - (a_q * R_star * R) / pow(params->c, 2));
    *E_total_phi = common_factor * ((R_phi - (R * v_q) / params->c) * velocity_factor);
    *E_total_z = common_factor * (R_z * velocity_factor);
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
    /* Преобразование в декартовы координаты в локальной системе */
    double x_a = ra * sin(theta_a) * cos(psi_a);
    double y_a = ra * sin(theta_a) * sin(psi_a);
    double z_a = ra * cos(theta_a);

    double x_q = rq * sin(theta_q) * cos(psi_q);
    double y_q = rq * sin(theta_q) * sin(psi_q);
    double z_q = rq * cos(theta_q);

    /* Преобразование в цилиндрическую систему (исходная) */
    double rho_a = sqrt(pow(params->rho0 + x_a, 2) + pow(y_a, 2));
    double phi_a = atan2(y_a, params->rho0 + x_a);
    if (phi_a < 0) phi_a += 2.0 * M_PI;

    double rho_q = sqrt(pow(params->rho0 + x_q, 2) + pow(y_q, 2));
    double phi_q = atan2(y_q, params->rho0 + x_q);
    if (phi_q < 0) phi_q += 2.0 * M_PI;

    /* Вектор от источника к наблюдателю */
    double R_rho = rho_a * cos(phi_a - phi_q) - rho_q;
    double R_phi = rho_a * sin(phi_a - phi_q);
    double R_z = z_a - z_q;
    double R = sqrt(pow(R_rho, 2) + pow(R_phi, 2) + pow(R_z, 2));

    /* Скорость и ускорение источника */
    double v_q = rho_q * params->omega;
    double a_q = -rho_q * pow(params->omega, 2);

    /* Радиус Лиенара-Вихерта (учет запаздывания) */
    double R_star = R - (R_phi * v_q) / params->c;

    /* Защита от деления на ноль */
    if (R_star < 1e-10 || R < 1e-10) {
        *E1_rho = *E1_phi = *E1_z = 0.0;
        *E2_rho = *E2_phi = *E2_z = 0.0;
        *E_total_rho = *E_total_phi = *E_total_z = 0.0;
        return;
    }

    /* Вычисление градиентного поля E1 */
    double common_factor1 = 1.0 / pow(R_star, 2);
    double velocity_factor1 = /*1.0 +*/ (R_rho * a_q) / pow(params->c, 2) /*- pow(v_q / params->c, 2)*/;

    *E1_rho = common_factor1 * (R_rho * velocity_factor1 / R_star);
    *E1_phi = common_factor1 * (R_phi * velocity_factor1 / R_star - v_q / params->c);
    *E1_z   = common_factor1 * (R_z   * velocity_factor1 / R_star);

    /* Вычисление поля самоиндукции E2 */
    double common_factor2 = 1.0 / pow(R_star, 2);
    double velocity_factor2 = (R / R_star) * (pow(v_q / params->c, 2) - (R_rho * a_q) / pow(params->c, 2) - 1.0) + 1.0;

    *E2_rho = common_factor2 * (- a_q * R / pow(params->c, 2));
    *E2_phi = common_factor2 * (v_q / params->c * velocity_factor2);
    *E2_z = 0;

    /* Вычисление суммарного поля */
    double common_factor = 1.0 / pow(R_star, 3);
    double velocity_factor = /*1.0 +*/ (R_rho * a_q) / pow(params->c, 2) /*- pow(v_q, 2) / pow(params->c, 2)*/;

    *E_total_rho = common_factor * (R_rho * velocity_factor - (a_q * R_star * R) / pow(params->c, 2));
    *E_total_phi = common_factor * ((R_phi - (R * v_q) / params->c) * velocity_factor);
    *E_total_z = common_factor * (R_z * velocity_factor);
}


// the charge distribution
static inline cubareal rho_q (cubareal r0, cubareal q)
{
    return 3 * (q) / (4*M_PI*(r0)*(r0)*(r0));
}

// интегрирование по координатам заряда источника потенциала
// E2 - индуктивная компонента массы
// В приближении малых скоростей ${}^{v} \big / {}_{c}\ll 1$
// и малых ускорений $a{{r}_{0}}\ll {{c}^{2}}$ и при игнорировании запаздывания
// Электрическое поле самоиндукции ($z$ компонента)
/*
\[{E}_{2}=\\
\int\limits_{{{r}_{q}}}\int\limits_{{{\varphi}_{q}}}\int\limits_{{{\theta}_{q}}}\\
{\left\{ -\frac{{a_z}}{{{c}^{2}}{{{R}_{0}}}} \right\}\\
{\rho \left( {{r}_{q}} \right){{r}_{q}}^{2}\sin \left( {{\theta }_{q}} \right)}\ }d{{\theta }_{q}}d{{\varphi }_{q}}d{{r}_{q}}\]
*/
static inline void int_q (const ProblemParams *params, cubareal q, cubareal psi_a, cubareal theta_a, cubareal ra, cubareal psi_q, cubareal theta_q, cubareal rq,
    double *int_q_E1_rho, double *int_q_E1_phi, double *int_q_E1_z,
    double *int_q_E2_rho, double *int_q_E2_phi, double *int_q_E2_z,
    double *int_q_E_rho,  double *int_q_E_phi,  double *int_q_E_z

 )
{
    cubareal r0 = params->R0;
    // В приближении малых скоростей ${}^{v}/{}_{c}\ll 1$
    // и малых ускорений $a{{r}_{0}}\ll {{c}^{2}}$
    // и при игнорировании запаздывания
    cubareal R = R0 (ra, theta_a, rq, theta_q, psi_q);

    //double mean_a = -params->rho0 * pow(params->omega, 2);

    /* Вычисление электрического поля */
    double E1_rho, E1_phi, E1_z;
    double E2_rho, E2_phi, E2_z;
    double E_rho,  E_phi,  E_z;

    compute_electric_field_with_delay(ra, theta_a, psi_a,
                           rq, theta_q, psi_q, params,
                          &E1_rho, &E1_phi, &E1_z,
                          &E2_rho, &E2_phi, &E2_z,
                          &E_rho,  &E_phi,  &E_z);

//printf("R = %e _R = %e\n", R, _R);
//printf("1/R = %e E_rho = %e\n", 1.0/R, E_rho);

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

    //return mean_a * rho_q(r0, q) * Sq(rq) * sin(theta_q) / R / pow(params->c, 2);
    //return rho_q(r0, q) * Sq(rq) * sin(theta_q) / R;
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
    double m_perp_B = 2 * U;                 /* Вариация типа В (теоретическое значение) */

    /* Вычисление поперечной массы */
    //double m_perp_A = integral[0] / Gamma / (params.c * params.c);  /* Вариация типа А */
    //double m = f/Gamma/ (params->c * params->c);

    double m1_rho = f1_rho / Gamma / (params->c * params->c);
    double m1_phi = f1_phi / Gamma / (params->c * params->c);
    double m1_z   = f1_z   / Gamma / (params->c * params->c);

    double m2_rho = f2_rho / Gamma / (params->c * params->c);
    double m2_phi = f2_phi / Gamma / (params->c * params->c);
    double m2_z   = f2_z   / Gamma / (params->c * params->c);

    double m_rho = f_rho / Gamma / (params->c * params->c);
    double m_phi = f_phi / Gamma / (params->c * params->c);
    double m_z   = f_z   / Gamma / (params->c * params->c);

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

