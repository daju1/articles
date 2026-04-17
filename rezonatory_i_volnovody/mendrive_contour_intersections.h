#ifndef MENDDRIVE_CONTOUR_INTERSECTION_H
#define MENDDRIVE_CONTOUR_INTERSECTION_H

// Вспомогательная: пересечение двух отрезков
// Возвращает 1, если пересекаются, и заполняет (out_x, out_y)
// Использует параметрическое представление:
//   p + t * r,    q + u * s,    t,u ∈ [0,1]
int segment_intersection(
    long double p_x, long double p_y,
    long double r_x, long double r_y,
    long double q_x, long double q_y,
    long double s_x, long double s_y,
    long double* out_x, long double* out_y,
    long double* sin_r_s, long double eps_sin_r_s
);

// Находит все пересечения ломаных cu и cv
// Возвращает число найденных точек, записывает их в intersections (память выделяет вызывающая сторона)
// eps_det — порог |det| для фильтрации (например, 1e-6)
int find_contour_intersections(
    const long double* cu_x, const long double* cu_y, int cu_n,
    const long double* cv_x, const long double* cv_y, int cv_n,
    point2d_t* intersections, int max_intersections,
    long double eps_det,
    long double eps_sin_r_s
);

#endif
