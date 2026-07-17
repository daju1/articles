#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mendrive_log.h"
#include "mendrive_point2d.h"
#include "mendrive_isolines.h"
#include "mendrive_isolines_traced.h"
#include "mendrive_characteristic_roots.h"
#include "mendrive_contour_intersections.h"
#include "mendrive_contour_intersections_sharp.h"

int find_characteristic_roots(
    long double kz_min, long double kz_max, int nk,
    long double sz_min, long double sz_max, int ns,
    det_contours_result_t* contours,
    characteristic_roots_t* result,
    long double eps_nan,
    long double eps_det,
    long double eps_sin_r_s,
    long double extrap_len,
    long double cos_max_angle,
    long double sin_min_angle,
    long double sin_max_angle,
    int window_size,
    int use_tracing,
    int min_isoline_points_count,
    long double isoline_merge_segments_epsilon,
    int sharp,
    long double local_angle_staircase_threshold,  // 0.3L
    long double total_angle_threshold,            // 0.3L
    long double concentration_threshold,          // 0.4L
    long double local_angle_sharp_threshold,      // 0.6L
    long double det_threshold
) {
#ifdef LOGGING
MPREC_LOG_DEBUG("find_characteristic_roots\n");
#endif
    if (!result) return -1;
    result->roots = NULL;
    result->n_roots = 0;

    // Шаг 1: Построение изолиний
    int status;
    if (use_tracing) {
        status = compute_det_contours_traced(
            kz_min, kz_max, nk, sz_min, sz_max, ns, contours, eps_nan
        );
    } else {
        status = compute_det_contours(
            kz_min, kz_max, nk, sz_min, sz_max, ns, contours, eps_nan,
            min_isoline_points_count, isoline_merge_segments_epsilon
        );
    }

#ifdef LOGGING
MPREC_LOG_DEBUG("find_characteristic_roots: compute_det_contours status = %d\n", status);
#endif

    if (status != 0) {
#ifdef LOGGING
MPREC_LOG_DEBUG("find_characteristic_roots: compute_det_contours failed. return -1\n");
#endif
        return -1;
    }

    // Шаг 2: Найдём все пересечения между Re=0 и Im=0 линиями
    const int MAX_INTERSECTIONS = 10000;
    point2d_t* intersections = (point2d_t*)malloc(MAX_INTERSECTIONS * sizeof(point2d_t));
    if (!intersections) {
        free_det_contours(contours);
#ifdef LOGGING
MPREC_LOG_DEBUG("find_characteristic_roots no intersections returns -1\n");
#endif
        return -1;
    }

    int total_count = 0;

    // Перебираем все пары ломаных: каждая Re=0 с каждой Im=0
    for (int i = 0; i < contours->n_re_contours && total_count < MAX_INTERSECTIONS; ++i) {
        for (int j = 0; j < contours->n_im_contours && total_count < MAX_INTERSECTIONS; ++j) {
            const contour_line_t* cu = &contours->re_zero[i];
            const contour_line_t* cv = &contours->im_zero[j];

            if (cu->n_points < 2 || cv->n_points < 2) continue;

            // Извлекаем координаты
            long double* cu_x = (long double*)malloc(cu->n_points * sizeof(long double));
            long double* cu_y = (long double*)malloc(cu->n_points * sizeof(long double));
            long double* cv_x = (long double*)malloc(cv->n_points * sizeof(long double));
            long double* cv_y = (long double*)malloc(cv->n_points * sizeof(long double));

            for (int k = 0; k < cu->n_points; ++k) {
                cu_x[k] = cu->points[k].x;
                cu_y[k] = cu->points[k].y;
            }
            for (int k = 0; k < cv->n_points; ++k) {
                cv_x[k] = cv->points[k].x;
                cv_y[k] = cv->points[k].y;
            }

            // Находим пересечения
            int count;
            if (sharp) {
                count = find_contour_intersections_with_corners(
                    cu_x, cu_y, cu->n_points,
                    cv_x, cv_y, cv->n_points,
                    intersections + total_count,
                    MAX_INTERSECTIONS - total_count,
                    eps_det, eps_sin_r_s,
                    extrap_len,
                    cos_max_angle,
                    sin_min_angle,
                    sin_max_angle,
                    window_size,
                    local_angle_staircase_threshold,  // 0.3L
                    total_angle_threshold,            // 0.3L
                    concentration_threshold,          // 0.4L
                    local_angle_sharp_threshold,      // 0.6L
                    det_threshold
                );
            }
            else
            {
                count = find_contour_intersections(
                    cu_x, cu_y, cu->n_points,
                    cv_x, cv_y, cv->n_points,
                    intersections + total_count,
                    MAX_INTERSECTIONS - total_count,
                    eps_det, eps_sin_r_s
                );
            }

            total_count += count;

            free(cu_x); free(cu_y);
            free(cv_x); free(cv_y);

            if (total_count >= MAX_INTERSECTIONS) break;
        }
    }
#ifdef LOGGING
MPREC_LOG_DEBUG("find_characteristic_roots total_count %d\n", total_count);
#endif
    if (total_count == 0) {
        free(intersections);
        result->roots = NULL;
        result->n_roots = 0;
        return 0;
    }

    // Убираем дубликаты (опционально)
    // ... (можно добавить позже)

    result->roots = intersections;
    result->n_roots = total_count;
#ifdef LOGGING
MPREC_LOG_DEBUG("find_characteristic_roots return 0\n");
#endif
    return 0;
}

// Освобождение памяти
void free_characteristic_roots(characteristic_roots_t* r) {
    if (r && r->roots) {
        free(r->roots);
        r->roots = NULL;
        r->n_roots = 0;
    }
}