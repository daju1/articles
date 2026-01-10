#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mendrive_point2d.h"
#include "mendrive_isolines.h"
#include "mendrive_contour_intersections.h"

#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>
#include <gsl/gsl_roots.h>

// === Внешние функции (вы должны убедиться, что они видны при линковке) ===
extern void det_init(const void* params);
extern void det_eval(long double kz, long double sz, long double* det_re, long double* det_im);

// === Вспомогательные структуры для уточнения ===
typedef struct {
    long double k0, k1, s0, s1;
    int use_re; // 1 — Re, 0 — Im
} edge_params_t;

static double F_scalar(double t, void *p) {
    edge_params_t *ep = (edge_params_t*)p;
    long double kz = ep->k0 + t * (ep->k1 - ep->k0);
    long double sz = ep->s0 + t * (ep->s1 - ep->s0);
    long double re, im;
    det_eval(kz, sz, &re, &im);
    return (double)(ep->use_re ? re : im);
}

// Уточнение пересечения ребра с уровнем 0
static int refine_edge_intersection(
    long double k0, long double s0,
    long double k1, long double s1,
    int use_re,
    long double t_guess,
    long double *t_refined,
    long double eps_abs
) {
    if (t_guess < 0.0L || t_guess > 1.0L) {
        *t_refined = t_guess;
        return 0;
    }

    edge_params_t params = {k0, k1, s0, s1, use_re};
    gsl_function F = {&F_scalar, &params};

    const gsl_root_fsolver_type *T = gsl_root_fsolver_brent;
    gsl_root_fsolver *s = gsl_root_fsolver_alloc(T);
    if (!s) return -1;

    // Расширяем интервал немного, но оставляем в [0,1]
    double t_low = fmax(0.0, t_guess - 0.1);
    double t_high = fmin(1.0, t_guess + 0.1);

    double f_low = F_scalar(t_low, &params);
    double f_high = F_scalar(t_high, &params);

    // Если знаки одинаковые — расширяем до [0,1]
    if (f_low * f_high > 0.0) {
        t_low = 0.0;
        t_high = 1.0;
        f_low = F_scalar(0.0, &params);
        f_high = F_scalar(1.0, &params);
        if (f_low * f_high > 0.0) {
            gsl_root_fsolver_free(s);
            *t_refined = t_guess;
            return -1; // нет корня
        }
    }

    int status = gsl_root_fsolver_set(s, &F, t_low, t_high);
    if (status != GSL_SUCCESS) {
        gsl_root_fsolver_free(s);
        *t_refined = t_guess;
        return -1;
    }

    int iter = 0, max_iter = 50;
    do {
        iter++;
        status = gsl_root_fsolver_iterate(s);
        if (status != GSL_SUCCESS) break;

        double t_mid = gsl_root_fsolver_root(s);
        double f_mid = F_scalar(t_mid, &params);
        status = gsl_root_test_residual(f_mid, eps_abs);
        if (status == GSL_SUCCESS) {
            *t_refined = (long double)t_mid;
            gsl_root_fsolver_free(s);
            return 0;
        }
    } while (status == GSL_CONTINUE && iter < max_iter);

    *t_refined = (long double)gsl_root_fsolver_root(s);
    gsl_root_fsolver_free(s);
    return 0;
}

// Макросы
#define IS_INVALID(x, eps) (isnan(x) || isinf(x) || fabsl(x) > (eps))

// Таблица Marching Squares
static const int edge_table[16] = {
    0x0, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00
};

// Сравнение точек для сборки ломаных
typedef struct {
    long double x, y;
} pt_ld;

static int compare_pts(const void* a, const void* b) {
    pt_ld* pa = (pt_ld*)a;
    pt_ld* pb = (pt_ld*)b;
    if (fabsl(pa->x - pb->x) > 1e-15L) return (pa->x > pb->x) ? 1 : -1;
    if (fabsl(pa->y - pb->y) > 1e-15L) return (pa->y > pb->y) ? 1 : -1;
    return 0;
}

// Основная функция построения изолиний
static int extract_contours_from_grid(
    const long double* F, int nk, int ns,
    const long double* kz_grid, const long double* sz_grid,
    long double level,
    long double eps_nan,
    int use_re, // ← добавлено: нужно для уточнения
    contour_line_t** out_contours,
    int* n_contours
) {
    *out_contours = NULL;
    *n_contours = 0;

    if (nk < 2 || ns < 2) return 0;

    pt_ld* segments = NULL;
    size_t seg_count = 0, seg_cap = 1024;
    segments = (pt_ld*)malloc(seg_cap * 2 * sizeof(pt_ld));
    if (!segments) return -1;

    for (int j = 0; j < ns - 1; ++j) {
        for (int i = 0; i < nk - 1; ++i) {
            int idx[4] = { i + j*nk, (i+1) + j*nk, (i+1) + (j+1)*nk, i + (j+1)*nk };
            long double f[4] = { F[idx[0]], F[idx[1]], F[idx[2]], F[idx[3]] };
            int valid[4] = {
                !IS_INVALID(f[0], eps_nan),
                !IS_INVALID(f[1], eps_nan),
                !IS_INVALID(f[2], eps_nan),
                !IS_INVALID(f[3], eps_nan)
            };

            int cubeindex = 0;
            if (valid[0] && f[0] < level) cubeindex |= 1;
            if (valid[1] && f[1] < level) cubeindex |= 2;
            if (valid[2] && f[2] < level) cubeindex |= 4;
            if (valid[3] && f[3] < level) cubeindex |= 8;

            int edges = edge_table[cubeindex];
            if (edges == 0) continue;

            pt_ld verts[4];
            int vcount = 0;

            // Нижнее ребро
            if (edges & 1 && valid[0] && valid[1]) {
                long double t_lin = (level - f[0]) / (f[1] - f[0]);
                long double t_refined = t_lin;
                refine_edge_intersection(
                    kz_grid[i], sz_grid[j],
                    kz_grid[i+1], sz_grid[j],
                    use_re,
                    (double)t_lin,
                    &t_refined,
                    1e-12
                );
                verts[vcount].x = kz_grid[i] + t_refined * (kz_grid[i+1] - kz_grid[i]);
                verts[vcount].y = sz_grid[j];
                vcount++;
            }
            // Правое ребро
            if (edges & 2 && valid[1] && valid[2]) {
                long double t_lin = (level - f[1]) / (f[2] - f[1]);
                long double t_refined = t_lin;
                refine_edge_intersection(
                    kz_grid[i+1], sz_grid[j],
                    kz_grid[i+1], sz_grid[j+1],
                    use_re,
                    (double)t_lin,
                    &t_refined,
                    1e-12
                );
                verts[vcount].x = kz_grid[i+1];
                verts[vcount].y = sz_grid[j] + t_refined * (sz_grid[j+1] - sz_grid[j]);
                vcount++;
            }
            // Верхнее ребро
            if (edges & 4 && valid[2] && valid[3]) {
                long double t_lin = (level - f[2]) / (f[3] - f[2]);
                long double t_refined = t_lin;
                refine_edge_intersection(
                    kz_grid[i+1], sz_grid[j+1],
                    kz_grid[i], sz_grid[j+1],
                    use_re,
                    (double)t_lin,
                    &t_refined,
                    1e-12
                );
                verts[vcount].x = kz_grid[i+1] + t_refined * (kz_grid[i] - kz_grid[i+1]);
                verts[vcount].y = sz_grid[j+1];
                vcount++;
            }
            // Левое ребро
            if (edges & 8 && valid[3] && valid[0]) {
                long double t_lin = (level - f[3]) / (f[0] - f[3]);
                long double t_refined = t_lin;
                refine_edge_intersection(
                    kz_grid[i], sz_grid[j+1],
                    kz_grid[i], sz_grid[j],
                    use_re,
                    (double)t_lin,
                    &t_refined,
                    1e-12
                );
                verts[vcount].x = kz_grid[i];
                verts[vcount].y = sz_grid[j+1] + t_refined * (sz_grid[j] - sz_grid[j+1]);
                vcount++;
            }

            for (int vi = 0; vi < vcount - 1; vi += 2) {
                if (seg_count >= seg_cap) {
                    seg_cap *= 2;
                    pt_ld* tmp = (pt_ld*)realloc(segments, seg_cap * 2 * sizeof(pt_ld));
                    if (!tmp) { free(segments); return -1; }
                    segments = tmp;
                }
                segments[2*seg_count] = verts[vi];
                segments[2*seg_count + 1] = verts[vi+1];
                seg_count++;
            }
        }
    }

    // === СБОРКА ЛОМАНЫХ С ДОПУСКОМ ===
    const long double eps_merge = 1e-10L;

    if (seg_count == 0) {
        free(segments);
        return 0;
    }

    typedef struct {
        pt_ld p;
        size_t seg_index;
        int is_end;
    } endpoint_t;

    endpoint_t* endpoints = (endpoint_t*)malloc(2 * seg_count * sizeof(endpoint_t));
    if (!endpoints) { free(segments); return -1; }

    for (size_t i = 0; i < seg_count; ++i) {
        endpoints[2*i].p = segments[2*i];
        endpoints[2*i].seg_index = i;
        endpoints[2*i].is_end = 0;

        endpoints[2*i + 1].p = segments[2*i + 1];
        endpoints[2*i + 1].seg_index = i;
        endpoints[2*i + 1].is_end = 1;
    }

    char* used_seg = (char*)calloc(seg_count, 1);
    if (!used_seg) { free(segments); free(endpoints); return -1; }

    contour_line_t* contours = NULL;
    int n_c = 0, cap_c = 16;
    contours = (contour_line_t*)malloc(cap_c * sizeof(contour_line_t));

    for (size_t start_i = 0; start_i < seg_count; ++start_i) {
        if (used_seg[start_i]) continue;

        point2d_t* pts = NULL;
        int pcount = 0, pcap = 64;
        pts = (point2d_t*)malloc(pcap * sizeof(point2d_t));
        if (!pts) goto cleanup;

        pt_ld current_start = segments[2*start_i];
        pt_ld current_end   = segments[2*start_i + 1];

        pts[pcount].kz = current_start.x;
        pts[pcount].sz = current_start.y;
        pcount++;
        pts[pcount].kz = current_end.x;
        pts[pcount].sz = current_end.y;
        pcount++;

        used_seg[start_i] = 1;

        // Продолжаем вперёд
        pt_ld tail = current_end;
        int extended;
        do {
            extended = 0;
            for (size_t j = 0; j < seg_count; ++j) {
                if (used_seg[j]) continue;
                pt_ld a = segments[2*j], b = segments[2*j + 1];
                if (fabsl(a.x - tail.x) < eps_merge && fabsl(a.y - tail.y) < eps_merge) {
                    if (pcount >= pcap) {
                        pcap *= 2;
                        point2d_t* tmp = (point2d_t*)realloc(pts, pcap * sizeof(point2d_t));
                        if (!tmp) goto cleanup;
                        pts = tmp;
                    }
                    pts[pcount].kz = b.x;
                    pts[pcount].sz = b.y;
                    pcount++;
                    tail = b;
                    used_seg[j] = 1;
                    extended = 1;
                    break;
                }
                if (fabsl(b.x - tail.x) < eps_merge && fabsl(b.y - tail.y) < eps_merge) {
                    if (pcount >= pcap) {
                        pcap *= 2;
                        point2d_t* tmp = (point2d_t*)realloc(pts, pcap * sizeof(point2d_t));
                        if (!tmp) goto cleanup;
                        pts = tmp;
                    }
                    pts[pcount].kz = a.x;
                    pts[pcount].sz = a.y;
                    pcount++;
                    tail = a;
                    used_seg[j] = 1;
                    extended = 1;
                    break;
                }
            }
        } while (extended);

        // Теперь пробуем "нарастить" назад
        pt_ld head = current_start;
        do {
            extended = 0;
            for (size_t j = 0; j < seg_count; ++j) {
                if (used_seg[j]) continue;
                pt_ld a = segments[2*j], b = segments[2*j + 1];
                if (fabsl(b.x - head.x) < eps_merge && fabsl(b.y - head.y) < eps_merge) {
                    if (pcount >= pcap) {
                        pcap *= 2;
                        point2d_t* tmp = (point2d_t*)realloc(pts, pcap * sizeof(point2d_t));
                        if (!tmp) goto cleanup;
                        pts = tmp;
                    }
                    memmove(&pts[1], &pts[0], pcount * sizeof(point2d_t));
                    pts[0].kz = a.x;
                    pts[0].sz = a.y;
                    pcount++;
                    head = a;
                    used_seg[j] = 1;
                    extended = 1;
                    break;
                }
                if (fabsl(a.x - head.x) < eps_merge && fabsl(a.y - head.y) < eps_merge) {
                    if (pcount >= pcap) {
                        pcap *= 2;
                        point2d_t* tmp = (point2d_t*)realloc(pts, pcap * sizeof(point2d_t));
                        if (!tmp) goto cleanup;
                        pts = tmp;
                    }
                    memmove(&pts[1], &pts[0], pcount * sizeof(point2d_t));
                    pts[0].kz = b.x;
                    pts[0].sz = b.y;
                    pcount++;
                    head = b;
                    used_seg[j] = 1;
                    extended = 1;
                    break;
                }
            }
        } while (extended);

        if (n_c >= cap_c) {
            cap_c *= 2;
            contour_line_t* tmp = (contour_line_t*)realloc(contours, cap_c * sizeof(contour_line_t));
            if (!tmp) goto cleanup;
            contours = tmp;
        }
        contours[n_c].points = pts;
        contours[n_c].n_points = pcount;
        n_c++;
        continue;

    cleanup:
        free(pts);
        free(contours);
        free(segments);
        free(endpoints);
        free(used_seg);
        return -1;
    }

    free(segments);
    free(endpoints);
    free(used_seg);

    *out_contours = contours;
    *n_contours = n_c;
    return 0;
}

// === Основная экспортная функция ===
int compute_det_contours(
    long double kz_min, long double kz_max, int nk,
    long double sz_min, long double sz_max, int ns,
    det_contours_result_t* result,
    long double eps_nan
) {
    if (!result || nk < 2 || ns < 2) return -1;

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

    int err1 = extract_contours_from_grid(det_re, nk, ns, kz_grid, sz_grid, 0.0L, eps_nan, 1, &result->re_zero, &result->n_re_contours);
    int err2 = extract_contours_from_grid(det_im, nk, ns, kz_grid, sz_grid, 0.0L, eps_nan, 0, &result->im_zero, &result->n_im_contours);

    free(kz_grid);
    free(sz_grid);
    free(det_re);
    free(det_im);

    if (err1 != 0 || err2 != 0) {
        free_det_contours(result);
        return -1;
    }

    return 0;
}

void free_det_contours(det_contours_result_t* r) {
    if (!r) return;
    if (r->re_zero) {
        for (int i = 0; i < r->n_re_contours; ++i) {
            free(r->re_zero[i].points);
        }
        free(r->re_zero);
        r->re_zero = NULL;
    }
    if (r->im_zero) {
        for (int i = 0; i < r->n_im_contours; ++i) {
            free(r->im_zero[i].points);
        }
        free(r->im_zero);
        r->im_zero = NULL;
    }
    r->n_re_contours = r->n_im_contours = 0;
}