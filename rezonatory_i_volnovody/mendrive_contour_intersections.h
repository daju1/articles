#ifndef MENDDRIVE_CONTOUR_INTERSECTION_H
#define MENDDRIVE_CONTOUR_INTERSECTION_H

typedef struct {
    long double kz;
    long double sz;
} point2d_t;

// Находит все пересечения ломаных cu и cv
// Возвращает число найденных точек, записывает их в intersections (память выделяет вызывающая сторона)
// eps_det — порог |det| для фильтрации (например, 1e-6)
int find_contour_intersections(
    const long double* cu_x, const long double* cu_y, int cu_n,
    const long double* cv_x, const long double* cv_y, int cv_n,
    point2d_t* intersections, int max_intersections,
    long double eps_det
);


#endif
