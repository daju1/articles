#include <math.h>
#include <complex.h>
#include "mendrive_contour_intersections.h"
#include <stdio.h> // This line includes the standard input/output library


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
    {
        // printf("t = %Lf, u = %Lf\n", t, u);
        return 0;
    }

    *out_x = p_x + t * r_x;
    *out_y = p_y + t * r_y;
    return 1;
}

// Функция: найти все пересечения двух ломаных (контурных линий)
int find_contour_intersections(
    const long double* cu_x, const long double* cu_y, int cu_n,
    const long double* cv_x, const long double* cv_y, int cv_n,
    point2d_t* intersections, int max_intersections,
    long double eps_det
) {
    if (cu_n < 2 || cv_n < 2 || !intersections || max_intersections <= 0) {
        // printf("find_contour_intersections wrong params\n");
        return 0;
    }

    int count = 0;

    // Проходим по всем парам отрезков
    for (int i = 0; i < cu_n - 1; ++i) {
        long double p_x = cu_x[i],   p_y = cu_y[i];
        long double r_x = cu_x[i+1] - p_x;
        long double r_y = cu_y[i+1] - p_y;

        for (int j = 0; j < cv_n - 1; ++j) {
            long double q_x = cv_x[j],   q_y = cv_y[j];
            long double s_x = cv_x[j+1] - q_x;
            long double s_y = cv_y[j+1] - q_y;

            long double x, y;
            if (!segment_intersection(p_x, p_y, r_x, r_y, q_x, q_y, s_x, s_y, &x, &y))
                continue;

            if (count >= max_intersections) {
                fprintf(stderr, "find_contour_intersections: max_intersections exceeded\n");
                return count;
            }

            intersections[count].kz = x;
            intersections[count].sz = y;
            count++;
        }
    }

    return count;
}
