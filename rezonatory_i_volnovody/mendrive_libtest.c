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
                                   const char* name,
                                   corner2d_t* sharp_corners,
                                   int max_sharp_corners);

int main() {
    void* h = dlopen("./mendrive.so", RTLD_NOW);
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

    params.a = 0.50000000000000000000000000000000000000;
    params.mu_l_xx = 1.0000000000000000000000000000000000000;
    params.mu_l_yy = 1.0000000000000000000000000000000000000;
    params.mu_l_zz = 1.0000000000000000000000000000000000000;
    params.mu_r_xx = 100.00000000000000000000000000000000000;
    params.mu_r_yy = 100.00000000000000000000000000000000000;
    params.mu_r_zz = 100.00000000000000000000000000000000000;
    params.mu_l_yz = 0.00000000000000000000000000000000000000;
    params.mu_l_zy = 0.00000000000000000000000000000000000000;
    params.mu_r_yz = 0.00000000000000000000000000000000000000;
    params.mu_r_zy = 0.00000000000000000000000000000000000000;
    params.sigma_e_l_xx = 5.6250000000000000000000000000000000000e17;
    params.sigma_e_l_yy = 5.6250000000000000000000000000000000000e17;
    params.sigma_e_l_zz = 5.6250000000000000000000000000000000000e17;
    params.sigma_e_r_xx = 9.9000000000000000000000000000000000000e12;
    params.sigma_e_r_yy = 9.9000000000000000000000000000000000000e12;
    params.sigma_e_r_zz = 9.9000000000000000000000000000000000000e12;
    params.sigma_m_l_xx = 0.00000000000000000000000000000000000000;
    params.sigma_m_l_yy = 0.00000000000000000000000000000000000000;
    params.sigma_m_l_zz = 0.00000000000000000000000000000000000000;
    params.sigma_m_r_xx = 6.4800000000000000000000000000000000000e7;
    params.sigma_m_r_yy = 6.4800000000000000000000000000000000000e7;
    params.sigma_m_r_zz = 6.4800000000000000000000000000000000000e7;
    params.eps_l_xx = 1.0000000000000000000000000000000000000;
    params.eps_l_yy = 1.0000000000000000000000000000000000000;
    params.eps_l_zz = 1.0000000000000000000000000000000000000;
    params.eps_r_xx = 16.000000000000000000000000000000000000;
    params.eps_r_yy = 16.000000000000000000000000000000000000;
    params.eps_r_zz = 16.000000000000000000000000000000000000;
    params.c = 2.9979245800000000000000000000000000000e10;
    params.omega = 1.0000000000000000000000000000000000000e7;

    det_init_fn(&params);

    // Тестируем построение изолиний
    det_contours_result_t contours = {0};
    int status = compute_det_contours_fn(
        -200.0L, 200.0L, 400,  // kz
        -10.0L,   10.0L, 400,  // sz
        &contours,
        1e300L
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
    corner2d_t sharp_corners[max_sharp_corners];
    // Тестируем линии Re=0
    //for (int i = 0; i < contours.n_re_contours; ++i) {
    for (int i = 5; i < 6/*contours.n_re_contours*/; ++i) {
        sprintf(name, "Re=0 line=%d", i);
        test_sharp_corners_fn(&contours.re_zero[i], -0.94, 0.5, 0.8, name, sharp_corners, max_sharp_corners);
        //break;
    }

    // Тестируем линии Im=0
//     for (int i = 0; i < contours.n_im_contours; ++i) {
//         sprintf(name, "Im=0 line=%d", i);
//         test_sharp_corners_fn(&contours.im_zero[i], -0.94, 0.5, 0.8, name, sharp_corners, max_sharp_corners);
//     }

    // Освобождаем память
    free_det_contours_fn(&contours);
    dlclose(h);
    printf("\nТест завершён успешно!\n");
    return 0;
}