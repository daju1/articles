#ifndef MENDRIVE_ISOLINES_TRACED_H
#define MENDRIVE_ISOLINES_TRACED_H

// Трассировка одной изолинии
int trace_isoline(
    long double kz_start,
    long double sz_start,
    int use_re, // 1 — Re(det), 0 — Im(det)
    long double step_size,      // начальный шаг (например, 0.1)
    long double tolerance,      // точность коррекции (например, 1e-12)
    long double kz_min, long double kz_max,
    long double sz_min, long double sz_max,
    point2d_t** points,
    int* n_points
);

int compute_det_contours_traced(
    long double kz_min, long double kz_max, int nk,
    long double sz_min, long double sz_max, int ns,
    det_contours_result_t* result,
    long double eps_nan
);

#endif