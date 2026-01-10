#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include "mendrive_point2d.h"
#include "mendrive_isolines.h"

static int extract_contours_from_grid(
    const long double* F, int nk, int ns,
    const long double* kz_grid, const long double* sz_grid,
    long double level,
    long double eps_nan,
    contour_line_t** out_contours,
    int* n_contours
);

// Внешние функции (вы должны убедиться, что они видны при линковке)
extern void det_init(const void* params);
extern void det_eval(long double kz, long double sz, long double* det_re, long double* det_im);

// Макросы для удобства
#define IS_INVALID(x, eps) (isnan(x) || isinf(x) || fabsl(x) > (eps))

// Внутренняя структура ячейки Marching Squares
static const int edge_table[16] = {
    0x0, 0x109, 0x203, 0x30a, 0x406, 0x50f, 0x605, 0x70c,
    0x80c, 0x905, 0xa0f, 0xb06, 0xc0a, 0xd03, 0xe09, 0xf00
};

// Сборка ломаных из отрезков
typedef struct {
    long double x, y;
} pt_ld;

static int compare_pts(const void* a, const void* b) {
    pt_ld* pa = (pt_ld*)a;
    pt_ld* pb = (pt_ld*)b;
    if (fabsl(pa->x - pb->x) > 1e-18L) return (pa->x > pb->x) ? 1 : -1;
    if (fabsl(pa->y - pb->y) > 1e-18L) return (pa->y > pb->y) ? 1 : -1;
    return 0;
}

// Основная функция построения изолиний
static int extract_contours_from_grid(
    const long double* F, int nk, int ns,
    const long double* kz_grid, const long double* sz_grid,
    long double level,
    long double eps_nan,
    contour_line_t** out_contours,
    int* n_contours
) {
    *out_contours = NULL;
    *n_contours = 0;

    if (nk < 2 || ns < 2) return 0;

    // Список всех отрезков
    pt_ld* segments = NULL;
    size_t seg_count = 0, seg_cap = 1024;
    segments = (pt_ld*)malloc(seg_cap * 2 * sizeof(pt_ld));
    if (!segments) return -1;

    // Обход всех ячеек
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

            // Формируем битовую маску
            int cubeindex = 0;
            if (valid[0] && f[0] < level) cubeindex |= 1;
            if (valid[1] && f[1] < level) cubeindex |= 2;
            if (valid[2] && f[2] < level) cubeindex |= 4;
            if (valid[3] && f[3] < level) cubeindex |= 8;

            int edges = edge_table[cubeindex];
            if (edges == 0) continue;

            // Интерполируем точки пересечения
            pt_ld verts[4];
            int vcount = 0;

            // Нижнее ребро
            if (edges & 1 && valid[0] && valid[1]) {
                long double t = (level - f[0]) / (f[1] - f[0]);
                verts[vcount].x = kz_grid[i] + t * (kz_grid[i+1] - kz_grid[i]);
                verts[vcount].y = sz_grid[j];
                vcount++;
            }
            // Правое ребро
            if (edges & 2 && valid[1] && valid[2]) {
                long double t = (level - f[1]) / (f[2] - f[1]);
                verts[vcount].x = kz_grid[i+1];
                verts[vcount].y = sz_grid[j] + t * (sz_grid[j+1] - sz_grid[j]);
                vcount++;
            }
            // Верхнее ребро
            if (edges & 4 && valid[2] && valid[3]) {
                long double t = (level - f[2]) / (f[3] - f[2]);
                verts[vcount].x = kz_grid[i+1] + t * (kz_grid[i] - kz_grid[i+1]);
                verts[vcount].y = sz_grid[j+1];
                vcount++;
            }
            // Левое ребро
            if (edges & 8 && valid[3] && valid[0]) {
                long double t = (level - f[3]) / (f[0] - f[3]);
                verts[vcount].x = kz_grid[i];
                verts[vcount].y = sz_grid[j+1] + t * (sz_grid[j] - sz_grid[j+1]);
                vcount++;
            }

            // Добавляем отрезки
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

    // === Сборка ломаных из отрезков ===
    if (seg_count == 0) {
        free(segments);
        return 0;
    }

    // Создаём хэш-подобный список концов
    pt_ld* all_pts = (pt_ld*)malloc(2 * seg_count * sizeof(pt_ld));
    if (!all_pts) { free(segments); return -1; }
    memcpy(all_pts, segments, 2 * seg_count * sizeof(pt_ld));
    qsort(all_pts, 2 * seg_count, sizeof(pt_ld), compare_pts);

    // Простая сборка: ищем цепочки
    char* used = (char*)calloc(seg_count, 1);
    if (!used) { free(segments); free(all_pts); return -1; }

    contour_line_t* contours = NULL;
    int n_c = 0, cap_c = 16;
    contours = (contour_line_t*)malloc(cap_c * sizeof(contour_line_t));

    for (size_t si = 0; si < seg_count; ++si) {
        if (used[si]) continue;

        // Начинаем новую ломаную
        point2d_t* pts = NULL;
        int pcount = 0, pcap = 64;
        pts = (point2d_t*)malloc(pcap * sizeof(point2d_t));
        if (!pts) goto cleanup;

        // Добавляем первый отрезок
        pts[pcount].kz = segments[2*si].x;
        pts[pcount].sz = segments[2*si].y;
        pcount++;
        pts[pcount].kz = segments[2*si + 1].x;
        pts[pcount].sz = segments[2*si + 1].y;
        pcount++;
        used[si] = 1;

        // Продолжаем вперёд
        pt_ld current = segments[2*si + 1];
        int extended;
        do {
            extended = 0;
            for (size_t sj = 0; sj < seg_count; ++sj) {
                if (used[sj]) continue;
                pt_ld a = segments[2*sj], b = segments[2*sj + 1];
                if (fabsl(a.x - current.x) < 1e-15L && fabsl(a.y - current.y) < 1e-15L) {
                    if (pcount >= pcap) {
                        pcap *= 2;
                        point2d_t* tmp = (point2d_t*)realloc(pts, pcap * sizeof(point2d_t));
                        if (!tmp) goto cleanup;
                        pts = tmp;
                    }
                    pts[pcount].kz = b.x;
                    pts[pcount].sz = b.y;
                    pcount++;
                    current = b;
                    used[sj] = 1;
                    extended = 1;
                    break;
                }
                if (fabsl(b.x - current.x) < 1e-15L && fabsl(b.y - current.y) < 1e-15L) {
                    if (pcount >= pcap) {
                        pcap *= 2;
                        point2d_t* tmp = (point2d_t*)realloc(pts, pcap * sizeof(point2d_t));
                        if (!tmp) goto cleanup;
                        pts = tmp;
                    }
                    pts[pcount].kz = a.x;
                    pts[pcount].sz = a.y;
                    pcount++;
                    current = a;
                    used[sj] = 1;
                    extended = 1;
                    break;
                }
            }
        } while (extended);

        // Сохраняем ломаную
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
        free(all_pts);
        free(used);
        return -1;
    }

    free(segments);
    free(all_pts);
    free(used);

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

    // Выделение памяти под сетку
    long double* kz_grid = (long double*)malloc(nk * sizeof(long double));
    long double* sz_grid = (long double*)malloc(ns * sizeof(long double));
    long double* det_re = (long double*)malloc(nk * ns * sizeof(long double));
    long double* det_im = (long double*)malloc(nk * ns * sizeof(long double));

    if (!kz_grid || !sz_grid || !det_re || !det_im) {
        free(kz_grid); free(sz_grid); free(det_re); free(det_im);
        return -1;
    }

    // Заполнение сетки
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

    // Извлечение контуров
    int err1 = extract_contours_from_grid(det_re, nk, ns, kz_grid, sz_grid, 0.0L, eps_nan, &result->re_zero, &result->n_re_contours);
    int err2 = extract_contours_from_grid(det_im, nk, ns, kz_grid, sz_grid, 0.0L, eps_nan, &result->im_zero, &result->n_im_contours);

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

// === Функции для ctypes ===
void get_contour_kz(const contour_line_t* c, int idx, long double** x, int* n) {
    if (!c || idx < 0 || idx >= 1) return; // ctypes вызывает по одной ломаной за раз
    *x = (long double*)malloc(c->n_points * sizeof(long double));
    for (int i = 0; i < c->n_points; ++i) {
        (*x)[i] = c->points[i].kz;
    }
    *n = c->n_points;
}

void get_contour_sz(const contour_line_t* c, int idx, long double** y, int* n) {
    if (!c || idx < 0 || idx >= 1) return;
    *y = (long double*)malloc(c->n_points * sizeof(long double));
    for (int i = 0; i < c->n_points; ++i) {
        (*y)[i] = c->points[i].sz;
    }
    *n = c->n_points;
}