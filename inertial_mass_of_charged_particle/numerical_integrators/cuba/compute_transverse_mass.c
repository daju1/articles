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

/* Прототип функции-интегранда */
extern int integrand(const int *, const double *, const int *, double *, void *);

int main() {
    /* Параметры задачи */
    ProblemParams params;
    
    /* Параметры, соответствующие эксперименту Менде */
    params.R0 = 0.0002778;        /* радиус сферы (половина R1) */
    params.rho0 = 0.0002778;      /* расстояние от оси вращения */
    params.omega = sqrt(9.489e8 / params.rho0);  /* угловая скорость */
    params.c = 299792458.0;       /* скорость света в м/с */
    
    params.R0 = 1;        /* радиус сферы (половина R1) */
    params.rho0 = 10;      /* расстояние от оси вращения */
    params.omega = 1;  /* угловая скорость */
    params.c = 1;       /* скорость света в м/с */
    
    /* Параметры интегрирования */
    const int ndim = 6;           /* 6-мерный интеграл */
    const int ncomp = 3;          /* 3 компоненты силы: F_rho, F_phi, F_z */
    const int nvec = 1;           /* Интегрирование по одной точке за раз */
    const double epsrel = 1e-3;   /* Относительная точность */
    const double epsabs = 1e-12;  /* Абсолютная точность */
    #define VERBOSE 2
    const int flags = VERBOSE;    /* Вывод промежуточной информации */
    const int mineval = 0;        /* Минимальное количество оценок */
    const int maxeval = 1000000;  /* Максимальное количество оценок */
    const int key = 0;            /* Ключ для выбора метода (0 = автоматический выбор) */
    
    /* Массивы для результатов */
    int nregions, neval, fail;
    double integral[3], error[3], prob[3];
    
    printf("===== Вычисление поперечной электромагнитной массы =====\n");
    printf("Параметры задачи:\n");
    printf("  R0    = %.6e м\n", params.R0);
    printf("  rho0  = %.6e м\n", params.rho0);
    printf("  omega = %.6e рад/с\n", params.omega);
    printf("\n");
    printf("Параметры интегрирования:\n");
    printf("  ndim     = %d\n", ndim);
    printf("  ncomp    = %d\n", ncomp);
    printf("  epsrel   = %.1e\n", epsrel);
    printf("  epsabs   = %.1e\n", epsabs);
    printf("  maxeval  = %d\n", maxeval);
    printf("\n");
    
    printf("Запуск вычислений...\n");

#define SEED 0
#define NSTART 1000
#define NINCREASE 500
#define NBATCH 1000
#define GRIDNO 0
#define STATEFILE NULL
#define SPIN NULL

    /* Вызов алгоритма Cuhre (детерминированный метод) */
    Vegas(ndim, ncomp, integrand, &params, nvec, epsrel, epsabs, 
          flags, SEED, mineval, maxeval, NSTART, NINCREASE, NBATCH, GRIDNO, STATEFILE, SPIN, 
          &neval, &fail, 
          integral, error, prob);

    /* Вычисление ускорения */
    double Gamma = params.rho0 * pow(params.omega, 2);
    
    /* Вычисление поперечной массы */
    double m_perp_A = integral[0] / Gamma;  /* Вариация типа А */
    double U = 3.0 / (5.0 * params.R0);     /* Энергия электрического поля */
    double m_perp_B = U / (params.c * params.c);  /* Вариация типа В (теоретическое значение) */
    
    /* Вывод результатов */
    printf("===== Результаты вычислений =====\n");
    printf("Сила самодействия:\n");
    printf("  F_rho = %.6e ± %.2e Н\n", integral[0], error[0]);
    printf("  F_phi = %.6e ± %.2e Н\n", integral[1], error[1]);
    printf("  F_z   = %.6e ± %.2e Н\n", integral[2], error[2]);
    printf("\n");
    printf("Ускорение: Gamma = %.6e м/с²\n", Gamma);
    printf("Поперечная масса (вариация А): m_perp_A = %.6e кг\n", m_perp_A);
    printf("Теоретическая масса (вариация В): m_perp_B = %.6e кг\n", m_perp_B);
    printf("Отношение m_perp_A / m_perp_B = %.6f\n", m_perp_A / m_perp_B);
    printf("\n");
    printf("Статистика:\n");
    printf("  Области интегрирования: %d\n", nregions);
    printf("  Оценки подынтегральной функции: %d\n", neval);
    printf("  Код завершения: %d\n", fail);
    printf("  χ² вероятность: %.6e\n", prob[0]);
    
    /* Проверка коэффициента 4/3 */
    double expected_ratio = 4.0 / 3.0;
    double observed_ratio = m_perp_A / m_perp_B;
    double deviation = fabs(observed_ratio - expected_ratio) / expected_ratio;
    
    printf("\n");
    printf("===== Анализ коэффициента 4/3 =====\n");
    printf("Ожидаемый коэффициент: %.6f\n", expected_ratio);
    printf("Наблюдаемый коэффициент: %.6f\n", observed_ratio);
    printf("Отклонение: %.6f%%\n", deviation * 100.0);
    
    if (deviation < epsrel) {
        printf("Результат подтверждает наличие коэффициента 4/3\n");
    } else {
        printf("Результат НЕ подтверждает наличие коэффициента 4/3\n");
    }
    
    return 0;
}