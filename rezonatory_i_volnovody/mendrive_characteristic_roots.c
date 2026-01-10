#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mendrive_point2d.h"
#include "mendrive_isolines.h"
#include "mendrive_isolines_traced.h"
#include "mendrive_characteristic_roots.h"
#include "mendrive_contour_intersections.h"
#include "mendrive_contour_intersections_sharp.h"

int find_characteristic_roots(
    long double kz_min, long double kz_max, int nk,
    long double sz_min, long double sz_max, int ns,
    characteristic_roots_t* result,
    long double eps_nan,
    long double eps_det,
    long double extrap_len,
    long double cos_max_angle,
    int use_tracing
) {
    if (!result) return -1;
    result->roots = NULL;
    result->n_roots = 0;

    det_contours_result_t contours = {0};

    // Шаг 1: Построение изолиний
    int status;
    if (use_tracing) {
        status = compute_det_contours_traced(
            kz_min, kz_max, nk, sz_min, sz_max, ns, &contours, eps_nan
        );
    } else {
        status = compute_det_contours(
            kz_min, kz_max, nk, sz_min, sz_max, ns, &contours, eps_nan
        );
    }

    if (status != 0) {
        return -1;
    }

    // Шаг 2: Найдём все пересечения между Re=0 и Im=0 линиями
    const int MAX_INTERSECTIONS = 10000;
    point2d_t* intersections = (point2d_t*)malloc(MAX_INTERSECTIONS * sizeof(point2d_t));
    if (!intersections) {
        free_det_contours(&contours);
        return -1;
    }

    int total_count = 0;

    // Перебираем все пары ломаных: каждая Re=0 с каждой Im=0
    for (int i = 0; i < contours.n_re_contours && total_count < MAX_INTERSECTIONS; ++i) {
        for (int j = 0; j < contours.n_im_contours && total_count < MAX_INTERSECTIONS; ++j) {
            const contour_line_t* cu = &contours.re_zero[i];
            const contour_line_t* cv = &contours.im_zero[j];

            if (cu->n_points < 2 || cv->n_points < 2) continue;

            // Извлекаем координаты
            long double* cu_x = (long double*)malloc(cu->n_points * sizeof(long double));
            long double* cu_y = (long double*)malloc(cu->n_points * sizeof(long double));
            long double* cv_x = (long double*)malloc(cv->n_points * sizeof(long double));
            long double* cv_y = (long double*)malloc(cv->n_points * sizeof(long double));

            for (int k = 0; k < cu->n_points; ++k) {
                cu_x[k] = cu->points[k].kz;
                cu_y[k] = cu->points[k].sz;
            }
            for (int k = 0; k < cv->n_points; ++k) {
                cv_x[k] = cv->points[k].kz;
                cv_y[k] = cv->points[k].sz;
            }

            // Находим пересечения
            int count = find_contour_intersections_with_corners(
                cu_x, cu_y, cu->n_points,
                cv_x, cv_y, cv->n_points,
                intersections + total_count,
                MAX_INTERSECTIONS - total_count,
                eps_det,
                extrap_len,
                cos_max_angle
            );

            total_count += count;

            free(cu_x); free(cu_y);
            free(cv_x); free(cv_y);

            if (total_count >= MAX_INTERSECTIONS) break;
        }
    }

    free_det_contours(&contours);

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