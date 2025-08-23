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
    double R0; /* Радиус сферы */
    double v;  /* Продольная скорость */
    double a;  /* продольное ускорение */
    double c;  /* Скорость света */
} ProblemParams;

/* Функция для вычисления электрического поля по Лиенару-Вихерту для продольного случая */
static void compute_electric_field(
    double ra, double theta_a, double phi_a,
    double rq, double theta_q, double phi_q,
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
    
    /* Вектор от источника к наблюдателю */
    double R_x = x_a - x_q;
    double R_y = y_a - y_q;
    double R_z = z_a - z_q;
    double R = sqrt(R_x*R_x + R_y*R_y + R_z*R_z);
    
    /* Скорость и ускорение источника */
    double v = params->v;  /* Продольная скорость */
    double a = params->a;  /* Продольное ускорение */
    
    /* Радиус Лиенара-Вихерта (учет запаздывания) */
    double R_star = R - (R_z * v) / params->c;
    
    /* Защита от деления на ноль */
    if (R_star < 1e-10 || R < 1e-10) {
        *E1_x = *E1_y = *E1_z = 0.0;
        *E2_x = *E2_y = *E2_z = 0.0;
        *E_total_x = *E_total_y = *E_total_z = 0.0;
        return;
    }
    
    /* Вычисление градиентного поля E1 */
    double common_factor1 = 1.0 / pow(R_star, 2);
    double velocity_factor1 = 1.0 + (R_z * a) / pow(params->c, 2) - pow(v, 2) / pow(params->c, 2);
    
    *E1_x = common_factor1 * (R_x * velocity_factor1 / R_star);
    *E1_y = common_factor1 * (R_y * velocity_factor1 / R_star);
    *E1_z = common_factor1 * (R_z * velocity_factor1 / R_star - v / params->c);
    
    /* Вычисление поля самоиндукции E2 */
    double common_factor2 = 1.0 / pow(R_star, 2);
    double velocity_factor2 = (R / R_star) * (pow(v, 2) / pow(params->c, 2) - (R_z * a) / pow(params->c, 2) - 1.0) + 1.0;
    
    *E2_x = common_factor2 * (-a * R_x / pow(params->c, 2));
    *E2_y = common_factor2 * (-a * R_y / pow(params->c, 2));
    *E2_z = common_factor2 * (v / params->c * velocity_factor2 - a * R_z / pow(params->c, 2));
    
    /* Вычисление суммарного поля */
    *E_total_x = *E1_x + *E2_x;
    *E_total_y = *E1_y + *E2_y;
    *E_total_z = *E1_z + *E2_z;
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
    double *int_q_E1_x, double *int_q_E1_y, double *int_q_E1_z,
    double *int_q_E2_x, double *int_q_E2_y, double *int_q_E2_z,
    double *int_q_E_x,  double *int_q_E_y,  double *int_q_E_z

 )
{
    cubareal r0 = params->R0;
    // В приближении малых скоростей ${}^{v}/{}_{c}\ll 1$
    // и малых ускорений $a{{r}_{0}}\ll {{c}^{2}}$
    // и при игнорировании запаздывания
    cubareal R = R0 (ra, theta_a, rq, theta_q, psi_q);

    /* Вычисление электрического поля */
    double E1_x, E1_y, E1_z;
    double E2_x, E2_y, E2_z;
    double E_x,  E_y,  E_z;

    compute_electric_field/*with_delay*/(ra, theta_a, psi_a,
                           rq, theta_q, psi_q, params,
                          &E1_x, &E1_y, &E1_z,
                          &E2_x, &E2_y, &E2_z,
                          &E_x,  &E_y,  &E_z);

//printf("R = %e _R = %e\n", R, _R);
//printf("1/R = %e E_x = %e\n", 1.0/R, E_x);

    double k_q = rho_q(r0, q) * Sq(rq) * sin(theta_q);

    *int_q_E1_x = k_q * E1_x;
    *int_q_E1_y = k_q * E1_y;
    *int_q_E1_z = k_q * E1_z;
    *int_q_E2_x = k_q * E2_x;
    *int_q_E2_y = k_q * E2_y;
    *int_q_E2_z = k_q * E2_z;
    *int_q_E_x  = k_q * E_x;
    *int_q_E_y  = k_q * E_y;
    *int_q_E_z  = k_q * E_z;

    //return mean_a * rho_q(r0, q) * Sq(rq) * sin(theta_q) / R / pow(params->c, 2);
    //return rho_q(r0, q) * Sq(rq) * sin(theta_q) / R;
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
    double k_a = rho_q(r0, q) * Sq(ra) * sin(theta_a);

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
    double m_perp_B = 2 * U;                 /* Вариация типа В (теоретическое значение) */

    /* Вычисление продольной массы */
    //double m_perp_A = integral[0] / Gamma / (params.c * params.c);  /* Вариация типа А */
    //double m = f/Gamma/ (params->c * params->c);

    double m1_x = f1_x / Gamma / (params->c * params->c);
    double m1_y = f1_y / Gamma / (params->c * params->c);
    double m1_z = f1_z / Gamma / (params->c * params->c);

    double m2_x = f2_x / Gamma / (params->c * params->c);
    double m2_y = f2_y / Gamma / (params->c * params->c);
    double m2_z = f2_z / Gamma / (params->c * params->c);

    double m_x = f_x / Gamma / (params->c * params->c);
    double m_y = f_y / Gamma / (params->c * params->c);
    double m_z = f_z / Gamma / (params->c * params->c);

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

