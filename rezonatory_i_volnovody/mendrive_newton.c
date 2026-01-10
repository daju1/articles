#include "mendrive_newton.h"
#include <math.h>
#include <complex.h>
// #define I (0.0L + 1.0L * I)

static inline long double complex cdiv(long double complex a, long double complex b) {
    long double d = creall(b) * creall(b) + cimagl(b) * cimagl(b);
    return (a * conjl(b)) / d;
}

int newton_adaptive_step(
    long double *kz, long double *sz,
    long double *step_re, long double *step_im,
    long double *f_abs_out,
    long double abs_m,
    long double step_m,
    long double delta_eps,
    long double f_abs_eps
) {
    long double f_re, f_im;
    det_eval(*kz, *sz, &f_re, &f_im);
    long double f_abs = sqrtl(f_re * f_re + f_im * f_im);
    *f_abs_out = f_abs;

    if (f_abs < f_abs_eps) return 2; // converged by |f|

    long double dfdkz_re, dfdkz_im, dfdsz_re, dfdsz_im;
    det_div_diff_kz_eval(*kz, *sz, &dfdkz_re, &dfdkz_im);
    det_div_diff_sz_eval(*kz, *sz, &dfdsz_re, &dfdsz_im);

    long double dfdkz_abs = sqrtl(dfdkz_re * dfdkz_re + dfdkz_im * dfdkz_im);
    long double dfdsz_abs = sqrtl(dfdsz_re * dfdsz_re + dfdsz_im * dfdsz_im);

    if (dfdkz_abs < 1e-30L || dfdsz_abs < 1e-30L) return -1;

    long double complex f_c = f_re + I * f_im;
    long double complex dfdkz_c = dfdkz_re + I * dfdkz_im;
    long double complex dfdsz_c = dfdsz_re + I * dfdsz_im;

    long double complex delta_kz_c = cdiv(f_c, dfdkz_c);
    long double complex delta_sz_c = cdiv(f_c, dfdsz_c);

    long double delta_re_re = creall(delta_kz_c);
    long double delta_re_im = cimagl(delta_kz_c);
    long double delta_im_re = creall(delta_sz_c);
    long double delta_im_im = cimagl(delta_sz_c);

    // Конвергенция по delta
    if (fabsl(delta_re_re) < delta_eps && fabsl(delta_re_im) < delta_eps &&
        fabsl(delta_im_re) < delta_eps && fabsl(delta_im_im) < delta_eps)
        return 1;

    // ——— адаптивный демпфирующий шаг ———
    // Re(kz)
    long double kz_new = *kz - (*step_re) * delta_re_re;
    long double f_re_new, f_im_new;
    det_eval(kz_new, *sz, &f_re_new, &f_im_new);
    long double f_abs_new = sqrtl(f_re_new * f_re_new + f_im_new * f_im_new);

    if (f_abs_new > abs_m * f_abs) {
        *step_re *= 0.1L;
        kz_new = *kz - (*step_re) * delta_re_re;
        det_eval(kz_new, *sz, &f_re_new, &f_im_new);
        f_abs_new = sqrtl(f_re_new * f_re_new + f_im_new * f_im_new);
    }
    if (*step_re < 0.9L) *step_re /= step_m;

    // Im(sz)
    long double sz_new = *sz - (*step_im) * delta_im_im;
    det_eval(*kz, sz_new, &f_re_new, &f_im_new);
    f_abs_new = sqrtl(f_re_new * f_re_new + f_im_new * f_im_new);

    if (f_abs_new > abs_m * f_abs) {
        *step_im *= 0.1L;
        sz_new = *sz - (*step_im) * delta_im_im;
        det_eval(*kz, sz_new, &f_re_new, &f_im_new);
        f_abs_new = sqrtl(f_re_new * f_re_new + f_im_new * f_im_new);
    }
    if (*step_im < 0.9L) *step_im /= step_m;

    *kz = kz_new;
    *sz = sz_new;

    return 0; // continue
}