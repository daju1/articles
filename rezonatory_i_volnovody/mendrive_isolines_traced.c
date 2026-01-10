#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mendrive_point2d.h"
#include "mendrive_isolines.h"

#include <gsl/gsl_multiroots.h>


// Внешние функции (вы должны убедиться, что они видны при линковке)
extern void det_init(const void* params);
extern void det_eval(long double kz, long double sz, long double* det_re, long double* det_im);

// Вспомогательная функция: вычисление F и градиента
static void eval_F_and_grad(
    long double kz, long double sz,
    int use_re,
    long double* F_val,
    long double* dF_dkz,
    long double* dF_dsz
) {
    const long double h = 1e-8L;

    long double re, im;
    det_eval(kz, sz, &re, &im);
    *F_val = (use_re ? re : im);

    // Численный градиент
    long double re1, im1, re2, im2;

    // ∂/∂kz
    det_eval(kz + h, sz, &re1, &im1);
    det_eval(kz - h, sz, &re2, &im2);
    *dF_dkz = (use_re ? (re1 - re2) : (im1 - im2)) / (2.0L * h);

    // ∂/∂sz
    det_eval(kz, sz + h, &re1, &im1);
    det_eval(kz, sz - h, &re2, &im2);
    *dF_dsz = (use_re ? (re1 - re2) : (im1 - im2)) / (2.0L * h);
}

// Коррекция точки на уровень F=0 методом Ньютона (1D)
static int correct_to_isoline(
    long double* kz, long double* sz,
    int use_re,
    long double tolerance,
    int max_iter
) {
    for (int iter = 0; iter < max_iter; ++iter) {
        long double F_val, dF_dkz, dF_dsz;
        eval_F_and_grad(*kz, *sz, use_re, &F_val, &dF_dkz, &dF_dsz);

        long double grad_norm2 = dF_dkz*dF_dkz + dF_dsz*dF_dsz;
        if (grad_norm2 < 1e-30L) return -1; // градиент нулевой

        long double delta = F_val / grad_norm2;
        long double kz_new = *kz - delta * dF_dkz;
        long double sz_new = *sz - delta * dF_dsz;

        if (fabsl(*kz - kz_new) < tolerance && fabsl(*sz - sz_new) < tolerance) {
            *kz = kz_new;
            *sz = sz_new;
            return 0;
        }

        *kz = kz_new;
        *sz = sz_new;
    }
    return -1; // не сошлось
}

// Основная функция трассировки
int trace_isoline(
    long double kz_start,
    long double sz_start,
    int use_re,
    long double step_size,
    long double tolerance,
    long double kz_min, long double kz_max,
    long double sz_min, long double sz_max,
    point2d_t** points,
    int* n_points
) {
    // Корректируем стартовую точку на изолинию
    long double kz = kz_start, sz = sz_start;
    if (correct_to_isoline(&kz, &sz, use_re, tolerance, 20) != 0) {
        return -1; // не удалось скорректировать
    }

    // Выделение памяти
    int cap = 1024;
    point2d_t* pts = (point2d_t*)malloc(cap * sizeof(point2d_t));
    if (!pts) return -1;

    pts[0].kz = kz;
    pts[0].sz = sz;
    int count = 1;

    long double direction = 1.0L; // направление трассировки

    for (int step = 0; step < 100000; ++step) {
        // Вычисляем градиент
        long double F_val, dF_dkz, dF_dsz;
        eval_F_and_grad(kz, sz, use_re, &F_val, &dF_dkz, &dF_dsz);

        // Касательный вектор: (-dF/dsz, dF/dkz)
        long double tx = -dF_dsz;
        long double ty =  dF_dkz;
        long double t_norm = sqrtl(tx*tx + ty*ty);
        if (t_norm < 1e-15L) break;

        tx /= t_norm;
        ty /= t_norm;

        // Метод Рунге-Кутты 4-го порядка (упрощённый: один шаг)
        long double kz1 = kz + direction * step_size * tx;
        long double sz1 = sz + direction * step_size * ty;

        // Корректируем новую точку на изолинию
        if (correct_to_isoline(&kz1, &sz1, use_re, tolerance, 10) != 0) {
            // Попробуем уменьшить шаг
            step_size *= 0.5L;
            if (step_size < 1e-12L) break;
            continue;
        }

        // Проверка границ
        if (kz1 < kz_min || kz1 > kz_max || sz1 < sz_min || sz1 > sz_max) {
            break;
        }

        // Проверка замыкания (если вернулись близко к началу)
        if (count > 10) {
            long double dx = kz1 - pts[0].kz;
            long double dy = sz1 - pts[0].sz;
            if (dx*dx + dy*dy < 1e-10L) {
                break; // замкнули контур
            }
        }

        // Добавляем точку
        if (count >= cap) {
            cap *= 2;
            point2d_t* tmp = (point2d_t*)realloc(pts, cap * sizeof(point2d_t));
            if (!tmp) { free(pts); return -1; }
            pts = tmp;
        }
        pts[count].kz = kz1;
        pts[count].sz = sz1;
        count++;

        // Обновляем текущую точку
        kz = kz1;
        sz = sz1;

        // Адаптивный шаг (увеличиваем, если всё хорошо)
        step_size = fmin(step_size * 1.1L, 1.0L);
    }

    *points = pts;
    *n_points = count;
    return 0;
}

// Новая функция: трассировка изолиний
int compute_det_contours_traced(
    long double kz_min, long double kz_max, int nk_fine,
    long double sz_min, long double sz_max, int ns_fine,
    det_contours_result_t* result,
    long double eps_nan
) {
    if (!result || nk_fine < 2 || ns_fine < 2) return -1;

    // === Используем ГРУБУЮ сетку для поиска стартовых точек ===
    int nk = (nk_fine < 200) ? nk_fine : 100; // не более 100x100
    int ns = (ns_fine < 200) ? ns_fine : 100;

    long double* kz_grid = (long double*)malloc(nk * sizeof(long double));
    long double* sz_grid = (long double*)malloc(ns * sizeof(long double));
    long double* det_re = (long double*)malloc(nk * ns * sizeof(long double));
    long double* det_im = (long double*)malloc(nk * ns * sizeof(long double));

    if (!kz_grid || !sz_grid || !det_re || !det_im) {
        free(kz_grid); free(sz_grid); free(det_re); free(det_im);
        return -1;
    }

    for (int i = 0; i < nk; ++i) {
        kz_grid[i] = kz_min + (kz_max - kz_min) * i / (nk - 1);
    }
    for (int j = 0; j < ns; ++j) {
        sz_grid[j] = sz_min + (sz_max - sz_min) * j / (ns - 1);
    }

    for (int j = 0; j < ns; ++j) {
        for (int i = 0; i < nk; ++i) {
            det_eval(kz_grid[i], sz_grid[j], &det_re[i + j*nk], &det_im[i + j*nk]);
        }
    }

    // === Параметры ограничения ===
    const int MAX_LINES = 50;
    const long double TRACE_STEP = 0.3L; // больше шаг → быстрее

    // === Re=0 ===
    contour_line_t* re_lines = NULL;
    int n_re = 0, cap_re = 16;
    re_lines = (contour_line_t*)malloc(cap_re * sizeof(contour_line_t));
    if (!re_lines) goto cleanup_alloc;

    char* visited = (char*)calloc(nk * ns, 1);
    if (!visited) goto cleanup_alloc;

    for (int j = 0; j < ns-1 && n_re < MAX_LINES; ++j) {
        for (int i = 0; i < nk-1 && n_re < MAX_LINES; ++i) {
            if (visited[i + j*nk]) continue;

            int found = 0;
            long double kz0 = 0.0L, sz0 = 0.0L;

            if (det_re[i + j*nk] * det_re[i+1 + j*nk] < 0) {
                kz0 = 0.5L * (kz_grid[i] + kz_grid[i+1]);
                sz0 = sz_grid[j];
                found = 1;
            }
            else if (det_re[i+1 + j*nk] * det_re[i+1 + (j+1)*nk] < 0) {
                kz0 = kz_grid[i+1];
                sz0 = 0.5L * (sz_grid[j] + sz_grid[j+1]);
                found = 1;
            }
            else if (det_re[i+1 + (j+1)*nk] * det_re[i + (j+1)*nk] < 0) {
                kz0 = 0.5L * (kz_grid[i] + kz_grid[i+1]);
                sz0 = sz_grid[j+1];
                found = 1;
            }
            else if (det_re[i + (j+1)*nk] * det_re[i + j*nk] < 0) {
                kz0 = kz_grid[i];
                sz0 = 0.5L * (sz_grid[j] + sz_grid[j+1]);
                found = 1;
            }

            if (found) {
                point2d_t* points = NULL;
                int n_points = 0;
                int status = trace_isoline(
                    kz0, sz0, 1,
                    TRACE_STEP, 1e-10L, // чуть менее строгая точность
                    kz_min, kz_max, sz_min, sz_max,
                    &points, &n_points
                );

                if (status == 0 && n_points > 1 && n_points < 5000) {
                    visited[i + j*nk] = 1;
                    if (n_re >= cap_re) {
                        cap_re *= 2;
                        contour_line_t* tmp = (contour_line_t*)realloc(re_lines, cap_re * sizeof(contour_line_t));
                        if (!tmp) goto cleanup_re;
                        re_lines = tmp;
                    }
                    re_lines[n_re].points = points;
                    re_lines[n_re].n_points = n_points;
                    n_re++;
                } else {
                    free(points);
                }
            }
        }
    }

    free(visited);

    // === Im=0 ===
    contour_line_t* im_lines = NULL;
    int n_im = 0, cap_im = 16;
    im_lines = (contour_line_t*)malloc(cap_im * sizeof(contour_line_t));
    if (!im_lines) goto cleanup_re;

    visited = (char*)calloc(nk * ns, 1);
    if (!visited) goto cleanup_im;

    for (int j = 0; j < ns-1 && n_im < MAX_LINES; ++j) {
        for (int i = 0; i < nk-1 && n_im < MAX_LINES; ++i) {
            if (visited[i + j*nk]) continue;

            int found = 0;
            long double kz0 = 0.0L, sz0 = 0.0L;

            if (det_im[i + j*nk] * det_im[i+1 + j*nk] < 0) {
                kz0 = 0.5L * (kz_grid[i] + kz_grid[i+1]);
                sz0 = sz_grid[j];
                found = 1;
            }
            else if (det_im[i+1 + j*nk] * det_im[i+1 + (j+1)*nk] < 0) {
                kz0 = kz_grid[i+1];
                sz0 = 0.5L * (sz_grid[j] + sz_grid[j+1]);
                found = 1;
            }
            else if (det_im[i+1 + (j+1)*nk] * det_im[i + (j+1)*nk] < 0) {
                kz0 = 0.5L * (kz_grid[i] + kz_grid[i+1]);
                sz0 = sz_grid[j+1];
                found = 1;
            }
            else if (det_im[i + (j+1)*nk] * det_im[i + j*nk] < 0) {
                kz0 = kz_grid[i];
                sz0 = 0.5L * (sz_grid[j] + sz_grid[j+1]);
                found = 1;
            }

            if (found) {
                point2d_t* points = NULL;
                int n_points = 0;
                int status = trace_isoline(
                    kz0, sz0, 0,
                    TRACE_STEP, 1e-10L,
                    kz_min, kz_max, sz_min, sz_max,
                    &points, &n_points
                );

                if (status == 0 && n_points > 1 && n_points < 5000) {
                    visited[i + j*nk] = 1;
                    if (n_im >= cap_im) {
                        cap_im *= 2;
                        contour_line_t* tmp = (contour_line_t*)realloc(im_lines, cap_im * sizeof(contour_line_t));
                        if (!tmp) goto cleanup_full;
                        im_lines = tmp;
                    }
                    im_lines[n_im].points = points;
                    im_lines[n_im].n_points = n_points;
                    n_im++;
                } else {
                    free(points);
                }
            }
        }
    }

    free(visited);
    free(kz_grid); free(sz_grid); free(det_re); free(det_im);

    result->re_zero = re_lines;
    result->n_re_contours = n_re;
    result->im_zero = im_lines;
    result->n_im_contours = n_im;
    return 0;

cleanup_full:
    for (int i = 0; i < n_im; ++i) free(im_lines[i].points);
    free(im_lines);
cleanup_im:
    free(im_lines);
cleanup_re:
    for (int i = 0; i < n_re; ++i) free(re_lines[i].points);
    free(re_lines);
cleanup_alloc:
    free(kz_grid); free(sz_grid); free(det_re); free(det_im);
    if (visited) free(visited);
    return -1;
}