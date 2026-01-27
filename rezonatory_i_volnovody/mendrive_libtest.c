// mendrive_libtest.c
#include <dlfcn.h>
#include <stdio.h>
#include <stdlib.h>

#include "mendrive_det.h"
#include "mendrive_point2d.h"
#include "mendrive_isolines.h"

typedef void (*det_init_t)(const mendrive_params_t* p);
typedef void (*det_eval_t) (long double kz, long double sz, long double *det_re, long double *det_im);

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

int main() {
    void* h = dlopen("./mendrive_libtest.so", RTLD_NOW);
    if (!h) {
        printf("dlopen failed: %s\n", dlerror());
        return 1;
    }

    // Загружаем функции
    det_init_t det_init_fn =
        (det_init_t)dlsym(h, "det_init");
    det_eval_t det_eval_fn =
        (det_eval_t)dlsym(h, "det_eval");

    compute_det_contours_t compute_det_contours_fn =
        (compute_det_contours_t)dlsym(h, "compute_det_contours");
    free_det_contours_t free_det_contours_fn =
        (free_det_contours_t)dlsym(h, "free_det_contours");

    test_sharp_corners_t test_sharp_corners_fn =
        (test_sharp_corners_t)dlsym(h, "test_sharp_corners");

    if (!compute_det_contours_fn || !free_det_contours_fn) {
        printf("Не найдены символы: %s\n", dlerror());
        dlclose(h);
        return 1;
    }

    printf("Библиотека загружена успешно!\n");

    mendrive_params_t params;

    long double sigma_0_d = 9.0e9;

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
    dlclose(h);
    printf("\nТест завершён успешно!\n");
    return 0;
}