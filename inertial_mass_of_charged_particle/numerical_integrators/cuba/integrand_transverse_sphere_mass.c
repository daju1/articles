#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "cuba.h"


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
    double *E_rho, double *E_phi, double *E_z
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
        *E_rho = 0.0;
        *E_phi = 0.0;
        *E_z = 0.0;
        return;
    }
    
    /* Основная формула Лиенара-Вихерта */
    double common_factor = 1.0 / pow(R_star, 3);
    double velocity_factor = 1.0 + (R_rho * a_q) / pow(params->c, 2) - pow(v_q, 2) / pow(params->c, 2);
    
    *E_rho = common_factor * (R_rho * velocity_factor + (a_q * R_star * R) / pow(params->c, 2));
    *E_phi = common_factor * ((R_phi - (R * v_q) / params->c) * velocity_factor);
    *E_z = common_factor * (R_z * velocity_factor);
}

/* Функция-интегранд для библиотеки Cuba */
int integrand(
    const int *ndim, 
    const double x[], 
    const int *ncomp, 
    double f[], 
    void *userdata
) {
    ProblemParams *params = (ProblemParams *)userdata;
    
    /* Преобразование координат из [0,1] в физические значения */
    double ra = params->R0 * pow(x[0], 1.0/3.0);  /* r = R0 * x^(1/3) для равномерного распределения объема */
    double theta_a = acos(1.0 - 2.0 * x[1]);
    double psi_a = 2.0 * M_PI * x[2];
    
    double rq = params->R0 * pow(x[3], 1.0/3.0);
    double theta_q = acos(1.0 - 2.0 * x[4]);
    double psi_q = 2.0 * M_PI * x[5];
    
    /* Вычисление электрического поля */
    double E_rho, E_phi, E_z;
    compute_electric_field(ra, theta_a, psi_a, rq, theta_q, psi_q, params, &E_rho, &E_phi, &E_z);
    
    /* Плотность заряда (равномерное распределение) */
    double rho_charge = 3.0 / (4.0 * M_PI * pow(params->R0, 3));  /* Для единичного заряда */
    double dq_a = rho_charge * pow(ra, 2) * sin(theta_a);  /* Заряд точки наблюдателя */
    double dq_q = rho_charge * pow(rq, 2) * sin(theta_q);  /* Заряд источника */
    
    /* Сила самодействия */
    double dF_rho = dq_a * E_rho * dq_q;
    double dF_phi = dq_a * E_phi * dq_q;
    double dF_z = dq_a * E_z * dq_q;
    
    /* Учет якобиана преобразования координат */
    double jacobian = pow(params->R0, 7) * 8.0 * pow(M_PI, 2) * pow(x[0], 2) * (1.0 - 2.0 * x[1]) * pow(x[3], 2) * (1.0 - 2.0 * x[4]);
    
    /* Сохранение результатов */
    f[0] = dF_rho * jacobian;  /* Радиальная компонента силы */
    f[1] = dF_phi * jacobian;  /* Азимутальная компонента силы */
    f[2] = dF_z * jacobian;    /* Z-компонента силы */
    
    return 0;
}