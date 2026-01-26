#ifndef MENDRIVE_CHARACTERISTIC_ROOTS_H
#define MENDRIVE_CHARACTERISTIC_ROOTS_H

// Результат: список корней
typedef struct {
    point2d_t* roots;
    int n_roots;
} characteristic_roots_t;

// Основная функция
int find_characteristic_roots(
    long double kz_min, long double kz_max, int nk,
    long double sz_min, long double sz_max, int ns,
    det_contours_result_t* contours,
    characteristic_roots_t* result,
    long double eps_nan,
    long double eps_det,
    long double extrap_len,
    long double cos_max_angle,
    long double sin_min_angle,
    long double sin_max_angle,
    int window_size,
    int use_tracing,  // 0 — Marching Squares, 1 — трассировка
    int min_isoline_points_count,
    long double isoline_merge_segments_epsilon,
    int sharp,
    long double local_angle_staircase_threshold,  // 0.3L
    long double total_angle_threshold,            // 0.3L
    long double concentration_threshold,          // 0.4L
    long double local_angle_sharp_threshold       // 0.6L
);

// Освобождение памяти
void free_characteristic_roots(characteristic_roots_t* r);

#endif