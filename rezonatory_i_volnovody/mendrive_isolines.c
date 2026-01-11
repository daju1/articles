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
    if (fabsl(pa->x - pb->x) > 1e-16L) return (pa->x > pb->x) ? 1 : -1;
    if (fabsl(pa->y - pb->y) > 1e-16L) return (pa->y > pb->y) ? 1 : -1;
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
#if 1
    // === Аналитический поиск пересечений с уровнем 0 через билинейную модель ===
    for (int j = 0; j < ns - 1; ++j) {
        for (int i = 0; i < nk - 1; ++i) {
            // Значения в углах ячейки
            long double f00 = F[i + j*nk];
            long double f10 = F[(i+1) + j*nk];
            long double f11 = F[(i+1) + (j+1)*nk];
            long double f01 = F[i + (j+1)*nk];

            // Проверка на NaN/Inf
            if (IS_INVALID(f00, eps_nan) || IS_INVALID(f10, eps_nan) ||
                IS_INVALID(f11, eps_nan) || IS_INVALID(f01, eps_nan)) continue;

            // Координаты углов
            long double k0 = kz_grid[i],     k1 = kz_grid[i+1];
            long double s0 = sz_grid[j],     s1 = sz_grid[j+1];

            // Билинейная модель: F(k,s) = A*k*s + B*k + C*s + D
            long double denom = (k1 - k0) * (s1 - s0);
            if (fabsl(denom) < 1e-18L) continue;

            long double A = (f00 - f10 - f01 + f11) / denom;
            long double B = (f10 - f00) / (k1 - k0) - A * s0;
            long double C = (f01 - f00) / (s1 - s0) - A * k0;
            long double D = f00 - B*k0 - C*s0 - A*k0*s0;

            // Сдвигаем уровень: ищем F = level → F - level = 0
            D -= level;

            pt_ld verts[4];
            int vcount = 0;

            // --- Нижнее ребро: s = s0 ---
            {
                long double denom_edge = A * s0 + B;
                if (fabsl(denom_edge) > 1e-18L) {
                    long double k = -(C * s0 + D) / denom_edge;
                    if (k >= k0 && k <= k1) {
                        verts[vcount].x = k;
                        verts[vcount].y = s0;
                        vcount++;
                    }
                }
            }

            // --- Правое ребро: k = k1 ---
            {
                long double denom_edge = A * k1 + C;
                if (fabsl(denom_edge) > 1e-18L) {
                    long double s = -(B * k1 + D) / denom_edge;
                    if (s >= s0 && s <= s1) {
                        verts[vcount].x = k1;
                        verts[vcount].y = s;
                        vcount++;
                    }
                }
            }

            // --- Верхнее ребро: s = s1 ---
            {
                long double denom_edge = A * s1 + B;
                if (fabsl(denom_edge) > 1e-18L) {
                    long double k = -(C * s1 + D) / denom_edge;
                    if (k >= k0 && k <= k1) {
                        verts[vcount].x = k;
                        verts[vcount].y = s1;
                        vcount++;
                    }
                }
            }

            // --- Левое ребро: k = k0 ---
            {
                long double denom_edge = A * k0 + C;
                if (fabsl(denom_edge) > 1e-18L) {
                    long double s = -(B * k0 + D) / denom_edge;
                    if (s >= s0 && s <= s1) {
                        verts[vcount].x = k0;
                        verts[vcount].y = s;
                        vcount++;
                    }
                }
            }
#if 0
            // Добавляем отрезки (по 2 точки)
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
#else
            // === КРИТИЧЕСКОЕ ИСПРАВЛЕНИЕ: правильно соединяем точки ===
            // В билинейной модели изолиния — гипербола, которая может пересекать
            // ячейку в 0, 2 или 4 точках. Но никогда в 1 или 3!
            // Однако из-за численных ошибок может быть 3 точки.
            // Поэтому мы сортируем точки по углам ячейки и соединяем соседние.

            if (vcount == 0 || vcount > 4) continue;

            // Сортируем точки по порядку обхода ячейки (нижнее -> правое -> верхнее -> левое)
            // Для простоты — просто добавляем все возможные отрезки между парами
            for (int a = 0; a < vcount; ++a) {
                for (int b = a + 1; b < vcount; ++b) {
                    // Добавляем отрезок между a и b
                    if (seg_count >= seg_cap) {
                        seg_cap *= 2;
                        pt_ld* tmp = (pt_ld*)realloc(segments, seg_cap * 2 * sizeof(pt_ld));
                        if (!tmp) { free(segments); return -1; }
                        segments = tmp;
                    }
                    segments[2*seg_count] = verts[a];
                    segments[2*seg_count + 1] = verts[b];
                    seg_count++;
                }
            }
#endif
        }
    }

#else
    #if 1
    // Обход всех ячеек
    for (int j = 0; j < ns - 1; ++j) {
        for (int i = 0; i < nk - 1; ++i) {
            // Значения в углах
            long double f00 = F[i + j*nk];
            long double f10 = F[(i+1) + j*nk];
            long double f11 = F[(i+1) + (j+1)*nk];
            long double f01 = F[i + (j+1)*nk];

            // Проверка валидности
            if (IS_INVALID(f00, eps_nan) || IS_INVALID(f10, eps_nan) ||
                IS_INVALID(f11, eps_nan) || IS_INVALID(f01, eps_nan)) continue;

            // Координаты углов
            long double k0 = kz_grid[i], k1 = kz_grid[i+1];
            long double s0 = sz_grid[j], s1 = sz_grid[j+1];

            pt_ld verts[4];
            int vcount = 0;

            // --- Нижнее ребро (между (k0,s0) и (k1,s0)) ---
            if ((f00 < level && level < f10 ) || ( f10 < level && level < f00 )) {
                long double t = (level - f00) / (f10 - f00);
                // Учитываем направление сетки
                long double k_min = k0 > k1 ? k1 : k0;
                long double k_max = k0 < k1 ? k1 : k0;
                long double k = k0 + t * (k1 - k0);
                if (k > k_min && k < k_max) {
                    verts[vcount].x = k;
                    verts[vcount].y = s0;
                    vcount++;
                }
            }

            // --- Правое ребро (между (k1,s0) и (k1,s1)) ---
            if ((f10 < level && f11 > level) || (f10 > level && f11 < level)) {
                long double t = (level - f10) / (f11 - f10);
                //printf("right t=%Lf f10=%Le f11=%Le\n", t, f10, f11); fflush(stdout);
                long double s_min = s0 > s1 ? s1 : s0;
                long double s_max = s0 < s1 ? s1 : s0;
                long double s = s0 + t * (s1 - s0);
                if (s > s_min && s < s_max) {
                    verts[vcount].x = k1;
                    verts[vcount].y = s;
                    vcount++;
                }
            }

            // --- Верхнее ребро (между (k1,s1) и (k0,s1)) ---
            if (( f11 < level && level < f01 ) || ( f01 < level && level < f11 )) {
                long double t = (level - f11) / (f01 - f11);
                long double k_min = k0 > k1 ? k1 : k0;
                long double k_max = k0 < k1 ? k1 : k0;
                long double k = k1 + t * (k0 - k1);
                if (k > k_min && k < k_max) {
                    verts[vcount].x = k;
                    verts[vcount].y = s1;
                    vcount++;
                }
            }

            // --- Левое ребро (между (k0,s1) и (k0,s0)) ---
            if ((f01 < level && level < f00 ) || (f00 < level && level < f01)) {
                long double t = (level - f00) / (f01 - f00);
                //printf("left  t=%Lf f00=%Le f01=%Le\n", t, f00, f01); fflush(stdout);
                long double s_min = s0 > s1 ? s1 : s0;
                long double s_max = s0 < s1 ? s1 : s0;
                long double s = s0 + t * (s1 - s0);
                if (s > s_min && s < s_max) {
                    verts[vcount].x = k0;
                    verts[vcount].y = s;
                    vcount++;
                }
            }

            // === КРИТИЧЕСКОЕ ИСПРАВЛЕНИЕ: правильно соединяем точки ===
            // В простом случае должно быть 0 или 2 точки, но из-за ошибок может быть 4
            if (vcount == 2) {
                // Добавляем один отрезок
                if (seg_count >= seg_cap) {
                    seg_cap *= 2;
                    pt_ld* tmp = (pt_ld*)realloc(segments, seg_cap * 2 * sizeof(pt_ld));
                    if (!tmp) { free(segments); return -1; }
                    segments = tmp;
                }
                segments[2*seg_count] = verts[0];
                segments[2*seg_count + 1] = verts[1];
                seg_count++;
            }
            /*else if (vcount == 4) {
                // Добавляем два отрезка: (0-1) и (2-3)
                // Но сначала отсортируем точки по порядку обхода ячейки
                // Для простоты — просто добавим все возможные пары
                for (int a = 0; a < 4; a += 2) {
                    if (seg_count >= seg_cap) {
                        seg_cap *= 2;
                        pt_ld* tmp = (pt_ld*)realloc(segments, seg_cap * 2 * sizeof(pt_ld));
                        if (!tmp) { free(segments); return -1; }
                        segments = tmp;
                    }
                    segments[2*seg_count] = verts[a];
                    segments[2*seg_count + 1] = verts[a+1];
                    seg_count++;
                }
            }*/
            // Если vcount == 1 или 3 — игнорируем (численная ошибка)
        }
    }
    #else
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
    #endif
#endif
    // === Сборка ломаных из отрезков ===
    if (seg_count == 0) {
        free(segments);
        return 0;
    }

#if 1
    // === УЛУЧШЕННАЯ СБОРКА ЛОМАНЫХ С ДОПУСКОМ ===
    const long double eps_merge = 1e-10L; // ← подберите под масштаб ваших данных

    // Создаём список всех концов отрезков
    typedef struct {
        pt_ld p;
        int seg_index;
        int is_end; // 0 = начало, 1 = конец
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

    // Основной цикл: для каждого неиспользованного отрезка — строим цепочку
    for (size_t start_i = 0; start_i < seg_count; ++start_i) {
        if (used_seg[start_i]) continue;

        point2d_t* pts = NULL;
        int pcount = 0, pcap = 64;
        pts = (point2d_t*)malloc(pcap * sizeof(point2d_t));
        if (!pts) goto cleanup;

        // Начинаем с этого отрезка
        pt_ld current_start = segments[2*start_i];
        pt_ld current_end   = segments[2*start_i + 1];

        pts[pcount].kz = current_start.x;
        pts[pcount].sz = current_start.y;
        pcount++;
        pts[pcount].kz = current_end.x;
        pts[pcount].sz = current_end.y;
        pcount++;

        used_seg[start_i] = 1;

        // Продолжаем вперёд от current_end
        pt_ld tail = current_end;
        int extended;
        do {
            extended = 0;
            for (size_t j = 0; j < seg_count; ++j) {
                if (used_seg[j]) continue;
                pt_ld a = segments[2*j], b = segments[2*j + 1];
                if (fabsl(a.x - tail.x) < eps_merge && fabsl(a.y - tail.y) < eps_merge) {
                    // Присоединяем b
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
                    // Присоединяем a
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

        // Теперь пробуем "нарастить" назад от current_start
        pt_ld head = current_start;
        do {
            extended = 0;
            for (size_t j = 0; j < seg_count; ++j) {
                if (used_seg[j]) continue;
                pt_ld a = segments[2*j], b = segments[2*j + 1];
                if (fabsl(b.x - head.x) < eps_merge && fabsl(b.y - head.y) < eps_merge) {
                    // Присоединяем a В НАЧАЛО
                    if (pcount >= pcap) {
                        pcap *= 2;
                        point2d_t* tmp = (point2d_t*)realloc(pts, pcap * sizeof(point2d_t));
                        if (!tmp) goto cleanup;
                        pts = tmp;
                    }
                    // Сдвигаем все точки вправо
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
                    // Присоединяем b В НАЧАЛО
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
#else
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
#endif
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
