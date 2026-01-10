// mendrive_newton.h — совместимо с вашим mendrive_det.h
#ifndef MENDDRIVE_NEWTON_H
#define MENDDRIVE_NEWTON_H

#include "mendrive_det.h"

// Аналитические производные (если у вас уже сгенерированы)
void det_div_diff_kz_eval(long double kz, long double sz,
                          long double *re, long double *im);
void det_div_diff_sz_eval(long double kz, long double sz,
                          long double *re, long double *im);

// Адаптивный шаг Ньютона (полностью на long double)
int newton_adaptive_step(
    long double *kz, long double *sz,
    long double *step_re, long double *step_im,
    long double *f_abs_out,
    long double abs_m,
    long double step_m,
    long double delta_eps,
    long double f_abs_eps
);

#endif