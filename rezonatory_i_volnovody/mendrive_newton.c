#include "mendrive_newton.h"
#include <math.h>
#include <complex.h>

/**
 * Исправленная версия newton_adaptive_step
 * Полностью воспроизводит логику Python класса newton_prec.find_newton_complex_root()
 *
 * Ключевые изменения:
 * 1. Теперь 4 шага вместо 2: re_d_re, im_d_re, re_d_im, im_d_im
 * 2. Каждый шаг имеет свой адаптивный множитель
 * 3. Внутренний цикл повторов (max_retries) при неудачном шаге
 * 4. f_abs обновляется после каждого успешного шага
 * 5. Координата обновляется только при успешном уменьшении |f|
 */

int newton_adaptive_step(
    long double *kz, long double *sz,
    long double *step_re_re, long double *step_im_re,
    long double *step_re_im, long double *step_im_im,
    long double *f_abs_out,
    long double abs_m,
    long double step_decrease,   // было step__m = 0.1 в Python
    long double step_increase,   // было step_m = 0.2 в Python (делим на него)
    long double delta_eps,
    long double f_abs_eps,
    int max_retries              // N = 10 в Python
) {
    // Вычисляем текущее значение функции
    long double f_re, f_im;
    det_eval(*kz, *sz, &f_re, &f_im);
    long double f_abs = sqrtl(f_re * f_re + f_im * f_im);
    *f_abs_out = f_abs;

    // Проверка сходимости по |f|
    if (f_abs < f_abs_eps) return 2; // converged by |f|

    // Вычисляем f/df для обоих направлений
    long double div_kz_re, div_kz_im, div_sz_re, div_sz_im;
    det_div_diff_kz_eval(*kz, *sz, &div_kz_re, &div_kz_im);
    det_div_diff_sz_eval(*kz, *sz, &div_sz_re, &div_sz_im);

    // Проверка на вырожденность производных
    long double div_kz_abs = sqrtl(div_kz_re * div_kz_re + div_kz_im * div_kz_im);
    long double div_sz_abs = sqrtl(div_sz_re * div_sz_re + div_sz_im * div_sz_im);
    if (div_kz_abs < 1e-30L || div_sz_abs < 1e-30L) return -1; // degenerate

    // Сохраняем дельты для проверки сходимости
    // delta = f / df (это то, что вычитаем из координаты)
    long double delta_re_re = div_kz_re;  // Re(f/df_kz) -> применяется к kz
    long double delta_im_re = div_kz_im;  // Im(f/df_kz) -> применяется к kz
    long double delta_re_im = div_sz_re;  // Re(f/df_sz) -> применяется к sz
    long double delta_im_im = div_sz_im;  // Im(f/df_sz) -> применяется к sz

    long double f_re_new, f_im_new, f_abs_new;
    long double kz_new, sz_new;
    int n;

    // ===== ШАГ 1: re_d_re - обновление kz по Re(f/df_kz) =====
    n = max_retries;
    while (n > 0) {
        n--;
        kz_new = *kz - (*step_re_re) * delta_re_re;
        det_eval(kz_new, *sz, &f_re_new, &f_im_new);
        f_abs_new = sqrtl(f_re_new * f_re_new + f_im_new * f_im_new);

        if (f_abs_new > abs_m * f_abs) {
            *step_re_re *= step_decrease;
            continue;
        } else {
            *kz = kz_new;
            f_abs = f_abs_new;
            if (*step_re_re < 0.9L) {
                *step_re_re /= step_increase;
            }
            break;
        }
    }

    // ===== ШАГ 2: im_d_re - обновление kz по Im(f/df_kz) =====
    // Пересчитываем f/df_kz с новым kz
    det_div_diff_kz_eval(*kz, *sz, &div_kz_re, &div_kz_im);
    delta_im_re = div_kz_im;

    n = max_retries;
    while (n > 0) {
        n--;
        kz_new = *kz - (*step_im_re) * delta_im_re;
        det_eval(kz_new, *sz, &f_re_new, &f_im_new);
        f_abs_new = sqrtl(f_re_new * f_re_new + f_im_new * f_im_new);

        if (f_abs_new > abs_m * f_abs) {
            *step_im_re *= step_decrease;
            continue;
        } else {
            *kz = kz_new;
            f_abs = f_abs_new;
            if (*step_im_re < 0.9L) {
                *step_im_re /= step_increase;
            }
            break;
        }
    }

    // ===== ШАГ 3: re_d_im - обновление sz по Re(f/df_sz) =====
    // Пересчитываем f/df_sz с новым kz
    det_div_diff_sz_eval(*kz, *sz, &div_sz_re, &div_sz_im);
    delta_re_im = div_sz_re;

    n = max_retries;
    while (n > 0) {
        n--;
        sz_new = *sz - (*step_re_im) * delta_re_im;
        det_eval(*kz, sz_new, &f_re_new, &f_im_new);
        f_abs_new = sqrtl(f_re_new * f_re_new + f_im_new * f_im_new);

        if (f_abs_new > abs_m * f_abs) {
            *step_re_im *= step_decrease;
            continue;
        } else {
            *sz = sz_new;
            f_abs = f_abs_new;
            if (*step_re_im < 0.9L) {
                *step_re_im /= step_increase;
            }
            break;
        }
    }

    // ===== ШАГ 4: im_d_im - обновление sz по Im(f/df_sz) =====
    // Пересчитываем f/df_sz с новым sz
    det_div_diff_sz_eval(*kz, *sz, &div_sz_re, &div_sz_im);
    delta_im_im = div_sz_im;

    n = max_retries;
    while (n > 0) {
        n--;
        sz_new = *sz - (*step_im_im) * delta_im_im;
        det_eval(*kz, sz_new, &f_re_new, &f_im_new);
        f_abs_new = sqrtl(f_re_new * f_re_new + f_im_new * f_im_new);

        if (f_abs_new > abs_m * f_abs) {
            *step_im_im *= step_decrease;
            continue;
        } else {
            *sz = sz_new;
            f_abs = f_abs_new;
            if (*step_im_im < 0.9L) {
                *step_im_im /= step_increase;
            }
            break;
        }
    }

    // Обновляем выходное значение |f|
    *f_abs_out = f_abs;

    // Проверка сходимости по дельтам
    if (fabsl(delta_re_re) < delta_eps &&
        fabsl(delta_im_re) < delta_eps &&
        fabsl(delta_re_im) < delta_eps &&
        fabsl(delta_im_im) < delta_eps) {
        return 1; // converged by delta
    }

    // Проверка сходимости по |f| после всех шагов
    if (f_abs < f_abs_eps) {
        return 2; // converged by |f|
    }

    return 0; // continue iterations
}
