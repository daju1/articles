#ifndef MENDDRIVE_SHARP_CONTOUR_INTERSECTION_H
#define MENDDRIVE_SHARP_CONTOUR_INTERSECTION_H

// Адаптивное определение острых углов
void find_sharp_corners(
    const long double* x, const long double* y, int n,
    long double cos_max_angle,
    long double sin_min_angle,
    long double sin_max_angle,
    char* is_sharp_cos,  // выходной массив: 1 если излом, 0 иначе
    char* is_sharp_sin,  // выходной массив: 1 если излом, 0 иначе
    long double* cosines,
    long double* sines
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
    long double cos_max_angle,       // порог остроты угла, например -0.94L
    long double sin_min_angle,       // порог остроты угла, например 0.5L
    long double sin_max_angle,       // порог остроты угла, например 0.5L
    int window_size,
    long double local_angle_staircase_threshold,  // 0.3L
    long double total_angle_threshold,            // 0.3L
    long double concentration_threshold,          // 0.4L
    long double local_angle_sharp_threshold,      // 0.6L
    long double det_threshold
);

// Тестовая функция для острых углов
int test_sharp_corners(const contour_line_t* line,
                       long double cos_max_angle,
                       long double sin_min_angle,
                       long double sin_max_angle,
                       int window_size,
                       long double local_angle_staircase_threshold,  // 0.3L
                       long double total_angle_threshold,            // 0.3L
                       long double concentration_threshold,          // 0.4L
                       long double local_angle_sharp_threshold,      // 0.6L
                       long double det_threshold,
                       const char* name,
                       corner2d_t* sharp_corners,
                       int max_sharp_corners
                       );

#endif