// mendrive_libtest.c
#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "mendrive_det.h"
#include "mendrive_log.h"
#include "mendrive_point2d.h"
#include "mendrive_isolines.h"

typedef void (*det_init_t)(const mendrive_params_t* p);
typedef void (*det_eval_t) (mendrive_scalar_t kz, mendrive_scalar_t sz, mendrive_scalar_t *det_re, mendrive_scalar_t *det_im);
typedef void (*det_eval_old_t) (mendrive_scalar_t kz, mendrive_scalar_t sz, mendrive_scalar_t *det_re, mendrive_scalar_t *det_im);
typedef void (*det_div_diff_kz_eval_t)(mendrive_scalar_t kz, mendrive_scalar_t sz, mendrive_scalar_t *div_re, mendrive_scalar_t *div_im);
typedef void (*det_div_diff_sz_eval_t)(mendrive_scalar_t kz, mendrive_scalar_t sz, mendrive_scalar_t *div_re, mendrive_scalar_t *div_im);

typedef void (*K_E_v_eval_t)(mendrive_scalar_t kz, mendrive_scalar_t sz,
        mendrive_complex_t *K_E_vacuum);
typedef void (*K_H_v_eval_t)(mendrive_scalar_t kz, mendrive_scalar_t sz,
        mendrive_complex_t *K_H_vacuum);

typedef void (*K_E_l_eval_t)(mendrive_scalar_t kz, mendrive_scalar_t sz,
        mendrive_complex_t *K_E_left_conductor);
typedef void (*K_H_l_eval_t)(mendrive_scalar_t kz, mendrive_scalar_t sz,
        mendrive_complex_t *K_H_left_conductor);

typedef void (*K_E_r_eval_t)(mendrive_scalar_t kz, mendrive_scalar_t sz,
        mendrive_complex_t *K_E_right_conductor) ;
typedef void (*K_H_r_eval_t)(mendrive_scalar_t kz, mendrive_scalar_t sz,
        mendrive_complex_t *K_H_right_conductor);

int (*solve_newton_root_t)(double *kz_guess, double *sz_guess,
                      double *kz_sol, double *sz_sol,
                      double epsabs, int max_iter);

typedef int  (*newton_adaptive_step_t)(
    mendrive_scalar_t *kz, mendrive_scalar_t *sz,
    mendrive_scalar_t *step_re_re, mendrive_scalar_t *step_im_re,
    mendrive_scalar_t *step_re_im, mendrive_scalar_t *step_im_im,
    mendrive_scalar_t *f_abs_out,
    mendrive_scalar_t step_decrease,
    mendrive_scalar_t step_increase,
    mendrive_scalar_t delta_eps,
    mendrive_scalar_t f_abs_eps,
    int max_retries
);

typedef int (*newton_adaptive_t) (
    mendrive_scalar_t *kz,
    mendrive_scalar_t *sz
);

/**
 * Solve with verbose output
 */
typedef int (*newton_classic_solve_verbose_t) (
    mendrive_scalar_t* kz_inout,
    mendrive_scalar_t* sz_inout
);

// Прототипы функций из вашей библиотеки
typedef int (*compute_det_contours_t)(
    long double, long double, int,
    long double, long double, int,
    det_contours_result_t*,
    long double,
    int,
    long double
);

typedef void (*free_det_contours_t)(det_contours_result_t*);

// Находим острые углы (предполагаем, что функция в той же библиотеке)
typedef void (*find_sharp_corners_t)(
    const long double*, const long double*, int, long double, char*, long double *
);

// Тестовая функция для острых углов
typedef int (*test_sharp_corners_t)(const contour_line_t* line,
                                   long double cos_max_angle,
                                   long double sin_min_angle,
                                   long double sin_max_angle,
                                   int window_size,
                                   long double local_angle_staircase_threshold,  // 0.3L
                                   long double total_angle_threshold,            // 0.3L
                                   long double concentration_threshold,          // 0.4L
                                   long double local_angle_sharp_threshold,      // 0.6L
                                   long double det_threshold,
                                   const char* name,
                                   corner2d_t* sharp_corners,
                                   int max_sharp_corners
                                   );

// ============================================================================
// Тест метода Ньютона
// ============================================================================
int test_det_eval(det_eval_t det_eval_fn, det_eval_old_t det_eval_old_fn) {
    printf("\n");
    printf("\nТЕСТ det_eval (адаптивный)\n");
    printf("\n");

    // Начальное приближение (взято из вашего графического решения)
    mendrive_scalar_t kz = 2.176e-5L;
    mendrive_scalar_t sz = 34.592L;

    mendrive_scalar_t det_re;
    mendrive_scalar_t det_im;

    det_eval_old_fn(kz, sz, &det_re, &det_im);
    printf("det_old = " MPREC_LOG_FMT_SCALAR " " MPREC_LOG_FMT_SCALAR "\n", det_re, det_im);

    det_eval_fn(kz, sz, &det_re, &det_im);
    printf("det = " MPREC_LOG_FMT_SCALAR " " MPREC_LOG_FMT_SCALAR "\n", det_re, det_im);
}

// ============================================================================
// Тест метода Ньютона
// ============================================================================
int test_newton_method(newton_adaptive_t newton_adaptive_fn) {
    printf("\n");
    printf("\nТЕСТ МЕТОДА НЬЮТОНА (адаптивный)\n");
    printf("\n");

    // Начальное приближение (взято из вашего графического решения)
    mendrive_scalar_t kz = 2.176e-5L;
    mendrive_scalar_t sz = 34.592L;

    return newton_adaptive_fn(&kz, &sz);
}

int test_newton_classic_method_verbose(newton_classic_solve_verbose_t newton_classic_solve_verbose_fn) {
    printf("\n");
    printf("\nТЕСТ МЕТОДА НЬЮТОНА (classic verbose)\n");
    printf("\n");

    // Начальное приближение (взято из вашего графического решения)
    mendrive_scalar_t kz = 2.176e-5L;
    mendrive_scalar_t sz = 34.592L;

    return newton_classic_solve_verbose_fn(&kz, &sz);
}

// ============================================================================
// Тест изолиний (как в оригинале)
// ============================================================================
static int test_contours(
    det_eval_t det_eval_fn,
    compute_det_contours_t compute_fn,
    free_det_contours_t free_fn,
    test_sharp_corners_t test_corners_fn
) {
    printf("\n");
    printf("\nТЕСТ ПОСТРОЕНИЯ ИЗОЛИНИЙ\n");
    printf("\n");

    det_contours_result_t contours = {0};
    int status = compute_fn(
        -4.0L, 4.0L, 1000,   // kz_min, kz_max, nk
        -5.0L, 5.0L, 1000,   // sz_min, sz_max, ns
        &contours,
        1e300L,              // eps_nan
        0,                   // min_points_count
        1e-10L               // eps_merge
    );

    if (status != 0) {
        printf("❌ Ошибка при построении изолиний (код %d)\n", status);
        return 1;
    }

    printf("✅ Изолинии построены успешно\n");
    printf("   Re=0: %d линий\n", contours.n_re_contours);
    printf("   Im=0: %d линий\n", contours.n_im_contours);

    // Тест острых углов на нескольких линиях
    const int max_sharp_corners = 100;
    corner2d_t sharp_corners[max_sharp_corners];
    char name[128];

    // Тестируем первые 2 линии Re=0
    int n_tested_re = (contours.n_re_contours > 2) ? 2 : contours.n_re_contours;
    for (int i = 0; i < n_tested_re; ++i) {
        snprintf(name, sizeof(name), "Re=0 line=%d", i);
        int n_corners = test_corners_fn(
            &contours.re_zero[i],
            -0.94L, 0.5L, 0.8L, 10,
            0.3L, 0.3L, 0.4L, 0.6L, 1.0L,
            name, sharp_corners, max_sharp_corners
        );
        printf("   %s: найдено острых углов = %d\n", name, n_corners);
    }

    // Тестируем первые 2 линии Im=0
    int n_tested_im = (contours.n_im_contours > 2) ? 2 : contours.n_im_contours;
    for (int i = 0; i < n_tested_im; ++i) {
        snprintf(name, sizeof(name), "Im=0 line=%d", i);
        int n_corners = test_corners_fn(
            &contours.im_zero[i],
            -0.94L, 0.5L, 0.8L, 10,
            0.3L, 0.3L, 0.4L, 0.6L, 1.0L,
            name, sharp_corners, max_sharp_corners
        );
        printf("   %s: найдено острых углов = %d\n", name, n_corners);
    }

    free_fn(&contours);
    printf("\n");
    return 0;
}


int main() {
    void* h = dlopen("./mendrive_libtest.so", RTLD_NOW);
    if (!h) {
        printf("dlopen failed: %s\n", dlerror());
        return 1;
    }

    // Загружаем функции
    #define LOAD_FUNC(name, type) \
        type name##_fn = (type)dlsym(h, #name); \
        if (!name##_fn) { \
            fprintf(stderr, "❌ Не найдена функция '%s': %s\n", #name, dlerror()); \
            dlclose(h); \
            return 1; \
        }

//     det_init_t det_init_fn =
//         (det_init_t)dlsym(h, "det_init");
//     det_eval_t det_eval_fn =
//         (det_eval_t)dlsym(h, "det_eval");

//     compute_det_contours_t compute_det_contours_fn =
//         (compute_det_contours_t)dlsym(h, "compute_det_contours");
//     free_det_contours_t free_det_contours_fn =
//         (free_det_contours_t)dlsym(h, "free_det_contours");

//     test_sharp_corners_t test_sharp_corners_fn =
//         (test_sharp_corners_t)dlsym(h, "test_sharp_corners");

//     if (!compute_det_contours_fn || !free_det_contours_fn) {
//         printf("Не найдены символы: %s\n", dlerror());
//         dlclose(h);
//         return 1;
//     }

    LOAD_FUNC(det_init, det_init_t);
    LOAD_FUNC(det_eval, det_eval_t);
    LOAD_FUNC(det_eval_old, det_eval_old_t);
    LOAD_FUNC(det_div_diff_kz_eval, det_div_diff_kz_eval_t);
    LOAD_FUNC(det_div_diff_sz_eval, det_div_diff_sz_eval_t);
//    LOAD_FUNC(solve_newton_root, solve_newton_root_t);
    LOAD_FUNC(newton_adaptive_step, newton_adaptive_step_t);
    LOAD_FUNC(newton_adaptive, newton_adaptive_t);
    LOAD_FUNC(newton_classic_solve_verbose, newton_classic_solve_verbose_t);
    LOAD_FUNC(compute_det_contours, compute_det_contours_t);
    LOAD_FUNC(free_det_contours, free_det_contours_t);
    LOAD_FUNC(test_sharp_corners, test_sharp_corners_t);

    printf("Библиотека загружена успешно!\n");

    mendrive_params_t params;

    mendrive_scalar_t sigma_0_d = 9.0e9;

    params.a = 0.1;
    params.b = 0.5;
    params.m = 0;
    params.mu_l_xx = 1.0;
    params.mu_l_yy = 1.0;
    params.mu_l_zz = 1.0;
    params.mu_r_xx = 10.0;
    params.mu_r_yy = 10.0;
    params.mu_r_zz = 10.0;
    params.mu_l_yz = 0.0;
    params.mu_l_zy = 0.0;
    params.mu_r_yz = 0.0;
    params.mu_r_zy = 0.0;
    params.sigma_e_l_xx = 62500000 * sigma_0_d;
    params.sigma_e_l_yy = 62500000 * sigma_0_d;
    params.sigma_e_l_zz = 62500000 * sigma_0_d;
    params.sigma_e_r_xx = 0.000110000 * sigma_0_d;
    params.sigma_e_r_yy = 0.000110000 * sigma_0_d;
    params.sigma_e_r_zz = 0.000110000 * sigma_0_d;
    params.sigma_m_l_xx = 0.0;
    params.sigma_m_l_yy = 0.0;
    params.sigma_m_l_zz = 0.0;
    params.sigma_m_r_xx = 0.000062500 * sigma_0_d;
    params.sigma_m_r_yy = 0.000062500 * sigma_0_d;
    params.sigma_m_r_zz = 0.000062500 * sigma_0_d;
    params.eps_l_xx = 1.0;
    params.eps_l_yy = 1.0;
    params.eps_l_zz = 1.0;
    params.eps_r_xx = 16.0;
    params.eps_r_yy = 16.0;
    params.eps_r_zz = 16.0;
    params.c = 2.9979245800e10;
    params.omega = 1.0e6;
    params.mu_0 = 1.0;
    params.epsilon_0 = 1.0;

    det_init_fn(&params);
    printf("✅ Параметры инициализированы\n\n");

    // === ЗАПУСК ТЕСТОВ ===
    int exit_code = 0;

#if 1
    // Тест 1: метод Ньютона
    if (test_newton_method(
            newton_adaptive_fn
        ) != 0) {
        exit_code = 1;
    }

    test_newton_classic_method_verbose(newton_classic_solve_verbose_fn);
#endif

    test_det_eval(det_eval_fn, det_eval_old_fn);

#if 0
    // Тест 2: изолинии и острые углы
    if (test_contours(
            det_eval_fn,
            compute_det_contours_fn,
            free_det_contours_fn,
            test_sharp_corners_fn
        ) != 0) {
        exit_code = 1;
    }
#endif

#if 1
#else
    // Тестируем построение изолиний
    det_contours_result_t contours = {0};
    int status = compute_det_contours_fn(
        -4.0L, 4.0L, 1000,  // kz
        -5.0L, 5.0L, 1000,  // sz
        &contours,
        1e300L,
        0,
        1e-10
    );

    if (status != 0) {
        printf("Ошибка при построении изолиний!\n");
        dlclose(h);
        return 1;
    }

    printf("Re=0: %d линий\n", contours.n_re_contours);
    printf("Im=0: %d линий\n", contours.n_im_contours);

    char name[128];
    const int max_sharp_corners = 100;
    int  window_size = 10;
    long double cos_max_angle = -0.94;
    long double local_angle_staircase_threshold = 0.3L;
    long double total_angle_threshold           = 0.3L;
    long double concentration_threshold         = 0.4L;
    long double local_angle_sharp_threshold     = 0.6L;
    long double det_threshold = 1.0;
    corner2d_t sharp_corners[max_sharp_corners];
    // Тестируем линии Re=0
    for (int i = 0; i < contours.n_re_contours; ++i) {
        sprintf(name, "Re=0 line=%d", i);
        test_sharp_corners_fn(&contours.re_zero[i],
        cos_max_angle,
        0.5, 0.8,
        window_size,
        local_angle_staircase_threshold,  // 0.3L
        total_angle_threshold,            // 0.3L
        concentration_threshold,          // 0.4L
        local_angle_sharp_threshold,      // 0.6L
        det_threshold,
        name, sharp_corners, max_sharp_corners);
        //break;
    }

    // Тестируем линии Im=0
    for (int i = 0; i < contours.n_im_contours; ++i) {
        sprintf(name, "Im=0 line=%d", i);
        test_sharp_corners_fn(&contours.im_zero[i],
        cos_max_angle,
        0.5, 0.8,
        window_size,
        local_angle_staircase_threshold,  // 0.3L
        total_angle_threshold,            // 0.3L
        concentration_threshold,          // 0.4L
        local_angle_sharp_threshold,      // 0.6L
        det_threshold,
        name, sharp_corners, max_sharp_corners
        );
    }

    // Освобождаем память
    free_det_contours_fn(&contours);
#endif
    dlclose(h);
    printf("\nТест завершён успешно!\n");
    return 0;
}