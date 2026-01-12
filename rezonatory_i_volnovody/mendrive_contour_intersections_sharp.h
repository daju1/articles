#ifndef MENDDRIVE_SHARP_CONTOUR_INTERSECTION_H
#define MENDDRIVE_SHARP_CONTOUR_INTERSECTION_H

// Адаптивное определение острых углов
void find_sharp_corners_adaptive(
    const long double* x, const long double* y, int n,
    char* is_sharp  // выходной массив: 1 если излом, 0 иначе
);

// Находит все пересечения ломаных cu и cv
// Возвращает число найденных точек, записывает их в intersections (память выделяет вызывающая сторона)
// eps_det — порог |det| для фильтрации (например, 1e-6)
int find_contour_intersections_with_corners(
    const long double* cu_x, const long double* cu_y, int cu_n,
    const long double* cv_x, const long double* cv_y, int cv_n,
    point2d_t* intersections, int max_intersections,
    long double eps_det,
    long double extrap_len,          // длина экстраполяции (например, 1e-3)
    long double cos_max_angle        // порог остроты угла, например -0.94L
);

#endif