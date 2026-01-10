#ifndef MENDRIVE_ISOLINES_H
#define MENDRIVE_ISOLINES_H

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

// Одна ломаная
typedef struct {
    point2d_t* points;
    int n_points;
} contour_line_t;

// Результат: Re=0 и Im=0 контуры
typedef struct {
    contour_line_t* re_zero;   // det_re == 0
    int n_re_contours;
    contour_line_t* im_zero;   // det_im == 0
    int n_im_contours;
} det_contours_result_t;

// Основная функция
int compute_det_contours(
    long double kz_min, long double kz_max, int nk,
    long double sz_min, long double sz_max, int ns,
    det_contours_result_t* result,
    long double eps_nan      // порог для NaN/Inf: если |val| > eps_nan → считаем invalid
);

// Освобождение памяти
void free_det_contours(det_contours_result_t* r);

#endif