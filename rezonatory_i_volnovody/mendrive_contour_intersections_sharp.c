#include <math.h>
#include <complex.h>
#include "mendrive_point2d.h"
#include "mendrive_contour_intersections_sharp.h"
#include <stdio.h> // This line includes the standard input/output library

// Возвращает true, если точка i в ломаной является "острым" изломом
static int is_sharp_corner(
    const long double* x, const long double* y, int n, int i,
    long double cos_max_angle  // например, -0.866L для угла > 150°
) {
    int m = 10;
    if (i <= m || i >= n - m - 1) return 0; // концевые точки не считаются

    long double v1x = x[i] - x[i-1-m];
    long double v1y = y[i] - y[i-1-m];
    long double v2x = x[i+1+m] - x[i];
    long double v2y = y[i+1+m] - y[i];

    long double len1 = sqrtl(v1x*v1x + v1y*v1y);
    long double len2 = sqrtl(v2x*v2x + v2y*v2y);

    if (len1 < 1e-18L || len2 < 1e-18L) return 0;

    long double dot = v1x * v2x + v1y * v2y;
    long double cos_angle = dot / (len1 * len2);

    // Угол "острый", если он сильно отличается от 180° (cos ≈ -1)
    // Например: если cos_angle > -0.9 → угол < ~154°, но мы хотим "излом" → лучше использовать отклонение от прямой
    // Альтернатива: смотреть на |cross| / (len1*len2) — синус угла
    // Но проще: считать изломом, если угол < threshold_deg (например, 160°)
    // Тогда: cos_angle > cos(threshold_rad)

    return (cos_angle > cos_max_angle); // например, cos_max_angle = -0.94L (~160°)
}

// Экстраполирует одно плечо угла назад или вперёд
// dir = -1 → экстраполяция "назад" (через [i-1] -> [i])
// dir = +1 → экстраполяция "вперёд" (через [i] -> [i+1])
static void extrapolate_segment(
    const long double* x, const long double* y, int i, int dir,
    long double extrap_len,
    long double* out_p_x, long double* out_p_y,
    long double* out_r_x, long double* out_r_y
) {
    long double dx, dy, len;
    if (dir == -1) {
        dx = x[i] - x[i-1];
        dy = y[i] - y[i-1];
    } else {
        dx = x[i+1] - x[i];
        dy = y[i+1] - y[i];
    }
    len = sqrtl(dx*dx + dy*dy);
    if (len < 1e-18L) {
        dx = dy = 0.0L;
        len = 1.0L;
    }

    // Единичный вектор направления
    long double ux = dx / len;
    long double uy = dy / len;

    // Точка начала экстраполированного отрезка — это угловая точка
    *out_p_x = x[i];
    *out_p_y = y[i];

    // Направление — продолжение плеча
    *out_r_x = ux * extrap_len;
    *out_r_y = uy * extrap_len;
}

// Вспомогательная: пересечение двух отрезков
// Возвращает 1, если пересекаются, и заполняет (out_x, out_y)
// Использует параметрическое представление:
//   p + t * r,    q + u * s,    t,u ∈ [0,1]
static int segment_intersection(
    long double p_x, long double p_y,
    long double r_x, long double r_y,
    long double q_x, long double q_y,
    long double s_x, long double s_y,
    long double* out_x, long double* out_y
) {
    long double r_cross_s = r_x * s_y - r_y * s_x;
    if ( fabsl(r_cross_s ) < 1e-18L) {
        // printf("parallel\n");
        return 0; // параллельны или совпадают
    }

    long double q_minus_p_x = q_x - p_x;
    long double q_minus_p_y = q_y - p_y;

    long double t = (q_minus_p_x * s_y - q_minus_p_y * s_x) / r_cross_s;
    long double u = (q_minus_p_x * r_y - q_minus_p_y * r_x) / r_cross_s;

    if (t < 0.0L || t > 1.0L || u < 0.0L || u > 1.0L)
    //if (t < -1.0L || t > 2.0L || u < -1.0L || u > 2.0L)
    {
        // printf("t = %Lf, u = %Lf\n", t, u);
        return 0;
    }

    *out_x = p_x + t * r_x;
    *out_y = p_y + t * r_y;
    return 1;
}

// Возвращаем 2, если отрезки почти параллельны и близки — возможно, касание
/*static int segment_intersection_or_tangent(
    long double p_x, long double p_y,
    long double r_x, long double r_y,
    long double q_x, long double q_y,
    long double s_x, long double s_y,
    long double* out_x, long double* out_y,
    long double eps_det,
    long double eps_dist   // например, 1e-4L — порог близости
) {
    long double r_cross_s = r_x * s_y - r_y * s_x;

    if (fabsl(r_cross_s) < eps_det) {
        // Параллельны — проверим расстояние от точки q до прямой p + t*r

        // Расстояние от точки q до прямой через p с направлением r
        long double dist_num = fabsl( (q_x - p_x)*r_y - (q_y - p_y)*r_x );
        long double len_r = sqrtl(r_x*r_x + r_y*r_y);
        if (len_r < eps_det) return 0;

        long double dist = dist_num / len_r;

        if (dist > eps_dist) return 0; // слишком далеко — не касание

        // Проверим, находится ли проекция точки q на прямую p + t*r в пределах отрезка [0,1]
        long double proj_t = ((q_x - p_x)*r_x + (q_y - p_y)*r_y) / (r_x*r_x + r_y*r_y);
        if (proj_t < 0.0L || proj_t > 1.0L) return 0;

        // Точка касания — проекция q на прямую p + t*r
        *out_x = p_x + proj_t * r_x;
        *out_y = p_y + proj_t * r_y;

        return 2; // "касание"
    }

    // Обычное пересечение
    long double q_minus_p_x = q_x - p_x;
    long double q_minus_p_y = q_y - p_y;

    long double t = (q_minus_p_x * s_y - q_minus_p_y * s_x) / r_cross_s;
    long double u = (q_minus_p_x * r_y - q_minus_p_y * r_x) / r_cross_s;

    if (t < 0.0L || t > 1.0L || u < 0.0L || u > 1.0L)
        return 0;

    *out_x = p_x + t * r_x;
    *out_y = p_y + t * r_y;
    return 1;
}*/

int find_contour_intersections_with_corners(
    const long double* cu_x, const long double* cu_y, int cu_n,
    const long double* cv_x, const long double* cv_y, int cv_n,
    point2d_t* intersections, int max_intersections,
    long double eps_det,
    long double extrap_len,          // длина экстраполяции (например, 1e-3)
    long double cos_max_angle        // порог остроты угла, например -0.94L
) {
    if (cu_n < 2 || cv_n < 2 || !intersections || max_intersections <= 0) {
        return 0;
    }

    int count = 0;

    // 1. Обычные пересечения отрезок-отрезок
    for (int i = 0; i < cu_n - 1; ++i) {
        long double p_x = cu_x[i], p_y = cu_y[i];
        long double r_x = cu_x[i+1] - p_x;
        long double r_y = cu_y[i+1] - p_y;

        for (int j = 0; j < cv_n - 1; ++j) {
            long double q_x = cv_x[j], q_y = cv_y[j];
            long double s_x = cv_x[j+1] - q_x;
            long double s_y = cv_y[j+1] - q_y;

            long double x, y;
 #if 0
            int res = segment_intersection_or_tangent(
                p_x, p_y, r_x, r_y,
                q_x, q_y, s_x, s_y,
                &x, &y,
                eps_det,
                1e-4L  // eps_dist — подберите под масштаб
            );
 #else
             int res = segment_intersection(
                p_x, p_y, r_x, r_y,
                q_x, q_y, s_x, s_y,
                &x, &y
            );
#endif

            /*if (res == 1 || res == 2) {  // 1 — пересечение, 2 — касание
                if (count >= max_intersections) goto overflow;
                intersections[count].kz = x;
                intersections[count].sz = y;
                count++;
            }*/
        }
    }

    // 2. Обработка острых изломов на cu → экстраполяция → пересечение с cv
    for (int i = 1; i < cu_n - 1; ++i) {
        if (!is_sharp_corner(cu_x, cu_y, cu_n, i, cos_max_angle)) continue;

        // Два плеча: назад и вперёд
        for (int dir = -1; dir <= 1; dir += 2) {
            long double p_x, p_y, r_x, r_y;
            extrapolate_segment(cu_x, cu_y, i, dir, extrap_len, &p_x, &p_y, &r_x, &r_y);

            // Пересекаем этот экстраполированный отрезок со всеми отрезками cv
            for (int j = 0; j < cv_n - 1; ++j) {
                long double q_x = cv_x[j], q_y = cv_y[j];
                long double s_x = cv_x[j+1] - q_x;
                long double s_y = cv_y[j+1] - q_y;

                long double x, y;
                // Важно: используем eps_det!
                long double r_cross_s = r_x * s_y - r_y * s_x;
                if (fabsl(r_cross_s) < eps_det) continue;

                long double qm_px = q_x - p_x;
                long double qm_py = q_y - p_y;
                long double t = (qm_px * s_y - qm_py * s_x) / r_cross_s;
                long double u = (qm_px * r_y - qm_py * r_x) / r_cross_s;

                // Экстраполированный отрезок: t ∈ [0, 1] (мы его задали длиной extrap_len)
                // Отрезок cv: u ∈ [0, 1]
                //if (t >= 0.0L && t <= 1.0L && u >= 0.0L && u <= 1.0L) {
                if (t >= -0.5L && t <= 1.5L && u >= -0.5L && u <= 1.5L) {
                    x = p_x + t * r_x;
                    y = p_y + t * r_y;
                    if (count >= max_intersections) goto overflow;
                    intersections[count].kz = x;
                    intersections[count].sz = y;
                    count++;
                }
            }
        }
    }

    // 3. Обработка острых изломов на cv → экстраполяция → пересечение с cu
    for (int j = 1; j < cv_n - 1; ++j) {
        if (!is_sharp_corner(cv_x, cv_y, cv_n, j, cos_max_angle)) continue;

        for (int dir = -1; dir <= 1; dir += 2) {
            long double p_x, p_y, r_x, r_y;
            extrapolate_segment(cv_x, cv_y, j, dir, extrap_len, &p_x, &p_y, &r_x, &r_y);

            for (int i = 0; i < cu_n - 1; ++i) {
                long double q_x = cu_x[i], q_y = cu_y[i];
                long double s_x = cu_x[i+1] - q_x;
                long double s_y = cu_y[i+1] - q_y;

                long double x, y;
                long double r_cross_s = r_x * s_y - r_y * s_x;
                if (fabsl(r_cross_s) < eps_det) continue;

                long double qm_px = q_x - p_x;
                long double qm_py = q_y - p_y;
                long double t = (qm_px * s_y - qm_py * s_x) / r_cross_s;
                long double u = (qm_px * r_y - qm_py * r_x) / r_cross_s;

                //if (t >= 0.0L && t <= 1.0L && u >= 0.0L && u <= 1.0L) {
                if (t >= -0.5L && t <= 1.5L && u >= -0.5L && u <= 1.5L) {
                    x = p_x + t * r_x;
                    y = p_y + t * r_y;
                    if (count >= max_intersections) goto overflow;
                    intersections[count].kz = x;
                    intersections[count].sz = y;
                    count++;
                }
            }
        }
    }

    return count;

overflow:
    fprintf(stderr, "find_contour_intersections_with_corners: max_intersections exceeded\n");
    return count;
}