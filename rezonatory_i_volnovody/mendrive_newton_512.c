#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <mpfr.h>
#include <gmp.h>

// Точность: 512 бит (~154 десятичных знака)
#define MPFR_PREC 512

// Прототипы внешних функций (должны быть реализованы в вашей библиотеке)
// Принимают и возвращают значения в формате строк для передачи в MPFR
extern void det_eval_mpfr_str(
    const char* kz_str, const char* sz_str,
    char* f_re_out, char* f_im_out,
    size_t buf_size
);

extern void det_jacobian_mpfr_str(
    const char* kz_str, const char* sz_str,
    char* dRe_dkz_out, char* dRe_dsz_out,
    char* dIm_dkz_out, char* dIm_dsz_out,
    size_t buf_size
);

// Вспомогательная функция: строка → MPFR
static void str_to_mpfr(mpfr_t x, const char* str) {
    mpfr_set_str(x, str, 10, MPFR_RNDN);
}

// Вспомогательная функция: MPFR → строка
static void mpfr_to_str(char* buf, size_t size, mpfr_t x) {
    mpfr_snprintf(buf, size, "%.150Rf", x);
}

// Классический метод Ньютона для системы 2×2
int newton_mpfr_512(
    char* kz_str,      // вход/выход: начальное приближение и результат
    char* sz_str,      // вход/выход
    int max_iter,
    mpfr_t f_abs_out   // выход: |f| на последней итерации
) {
    // Инициализация переменных
    mpfr_t kz, sz, kz_new, sz_new;
    mpfr_t f_re, f_im, f_abs;
    mpfr_t dRe_dkz, dRe_dsz, dIm_dkz, dIm_dsz;
    mpfr_t detJ, delta_kz, delta_sz;
    mpfr_t alpha, f_abs_new, kz_trial, sz_trial;

    mpfr_inits2(MPFR_PREC, kz, sz, kz_new, sz_new,
                f_re, f_im, f_abs,
                dRe_dkz, dRe_dsz, dIm_dkz, dIm_dsz,
                detJ, delta_kz, delta_sz,
                alpha, f_abs_new, kz_trial, sz_trial, (mpfr_ptr)0);

    // Загрузка начальных значений
    str_to_mpfr(kz, kz_str);
    str_to_mpfr(sz, sz_str);

    mpfr_set_ui(alpha, 1, MPFR_RNDN);  // начальный шаг линейного поиска

    char buf[512];
    int converged = 0;

    for (int iter = 0; iter < max_iter; iter++) {
        // 1. Вычисляем значение функции f(kz, sz)
        mpfr_to_str(buf, sizeof(buf), kz);
        char sz_buf[512];
        mpfr_to_str(sz_buf, sizeof(sz_buf), sz);

        char f_re_str[512], f_im_str[512];
        det_eval_mpfr_str(buf, sz_buf, f_re_str, f_im_str, sizeof(f_re_str));

        str_to_mpfr(f_re, f_re_str);
        str_to_mpfr(f_im, f_im_str);

        // |f| = sqrt(Re² + Im²)
        mpfr_sqr(f_abs, f_re, MPFR_RNDN);
        mpfr_fma(f_abs, f_im, f_im, f_abs, MPFR_RNDN);  // f_abs += f_im²
        mpfr_sqrt(f_abs, f_abs, MPFR_RNDN);

        // Проверка сходимости
        mpfr_t eps;
        mpfr_init2(eps, MPFR_PREC);
        mpfr_set_str(eps, "1e-100", 10, MPFR_RNDN);

        if (mpfr_cmp(f_abs, eps) < 0) {
            converged = 1;
            mpfr_set(f_abs_out, f_abs, MPFR_RNDN);
            break;
        }

        // 2. Вычисляем якобиан
        char dRe_dkz_str[512], dRe_dsz_str[512];
        char dIm_dkz_str[512], dIm_dsz_str[512];

        det_jacobian_mpfr_str(buf, sz_buf,
            dRe_dkz_str, dRe_dsz_str, dIm_dkz_str, dIm_dsz_str,
            sizeof(dRe_dkz_str));

        str_to_mpfr(dRe_dkz, dRe_dkz_str);
        str_to_mpfr(dRe_dsz, dRe_dsz_str);
        str_to_mpfr(dIm_dkz, dIm_dkz_str);
        str_to_mpfr(dIm_dsz, dIm_dsz_str);

        // 3. Вычисляем определитель якобиана
        // detJ = dRe/dkz * dIm/dsz - dRe/dsz * dIm/dkz
        mpfr_mul(detJ, dRe_dkz, dIm_dsz, MPFR_RNDN);
        mpfr_fms(detJ, dRe_dsz, dIm_dkz, detJ, MPFR_RNDN);  // detJ -= dRe_dsz * dIm_dkz

        // Проверка на вырожденность
        mpfr_abs(detJ, detJ, MPFR_RNDN);
        mpfr_set_str(eps, "1e-120", 10, MPFR_RNDN);
        if (mpfr_cmp(detJ, eps) < 0) {
            fprintf(stderr, "⚠️ Якобиан вырожден на итерации %d\n", iter);
            mpfr_clears(eps, (mpfr_ptr)0);
            mpfr_clears(kz, sz, kz_new, sz_new, f_re, f_im, f_abs,
                        dRe_dkz, dRe_dsz, dIm_dkz, dIm_dsz, detJ,
                        delta_kz, delta_sz, alpha, f_abs_new,
                        kz_trial, sz_trial, (mpfr_ptr)0);
            return -1;
        }

        // 4. Решаем систему J·Δx = -f
        // Δkz = (-f_re * dIm/dsz + f_im * dRe/dsz) / detJ
        mpfr_mul(delta_kz, f_re, dIm_dsz, MPFR_RNDN);
        mpfr_neg(delta_kz, delta_kz, MPFR_RNDN);
        mpfr_fma(delta_kz, f_im, dRe_dsz, delta_kz, MPFR_RNDN);
        mpfr_div(delta_kz, delta_kz, detJ, MPFR_RNDN);

        // Δsz = ( f_re * dIm/dkz - f_im * dRe/dkz) / detJ
        mpfr_mul(delta_sz, f_re, dIm_dkz, MPFR_RNDN);
        mpfr_fms(delta_sz, f_im, dRe_dkz, delta_sz, MPFR_RNDN);
        mpfr_div(delta_sz, delta_sz, detJ, MPFR_RNDN);

        // 5. Линейный поиск для устойчивости
        mpfr_set(alpha, mpfr_const_one, MPFR_RNDN);  // alpha = 1.0

        for (int ls_iter = 0; ls_iter < 20; ls_iter++) {
            // kz_trial = kz - alpha * delta_kz
            mpfr_mul(kz_trial, alpha, delta_kz, MPFR_RNDN);
            mpfr_sub(kz_trial, kz, kz_trial, MPFR_RNDN);

            // sz_trial = sz - alpha * delta_sz
            mpfr_mul(sz_trial, alpha, delta_sz, MPFR_RNDN);
            mpfr_sub(sz_trial, sz, sz_trial, MPFR_RNDN);

            // Вычисляем |f| в пробной точке
            mpfr_to_str(buf, sizeof(buf), kz_trial);
            mpfr_to_str(sz_buf, sizeof(sz_buf), sz_trial);
            det_eval_mpfr_str(buf, sz_buf, f_re_str, f_im_str, sizeof(f_re_str));

            mpfr_t f_re_t, f_im_t, f_abs_t;
            mpfr_inits2(MPFR_PREC, f_re_t, f_im_t, f_abs_t, (mpfr_ptr)0);
            str_to_mpfr(f_re_t, f_re_str);
            str_to_mpfr(f_im_t, f_im_str);
            mpfr_sqr(f_abs_t, f_re_t, MPFR_RNDN);
            mpfr_fma(f_abs_t, f_im_t, f_im_t, f_abs_t, MPFR_RNDN);
            mpfr_sqrt(f_abs_t, f_abs_t, MPFR_RNDN);

            // Принимаем шаг, если |f| уменьшилась
            if (mpfr_cmp(f_abs_t, f_abs) < 0 || mpfr_cmp_ui(alpha, 1e-10) < 0) {
                mpfr_set(kz, kz_trial, MPFR_RNDN);
                mpfr_set(sz, sz_trial, MPFR_RNDN);
                mpfr_set(f_abs, f_abs_t, MPFR_RNDN);
                mpfr_clears(f_re_t, f_im_t, f_abs_t, (mpfr_ptr)0);
                break;
            }

            // Уменьшаем шаг
            mpfr_mul_ui(alpha, alpha, 5, MPFR_RNDN);
            mpfr_div_ui(alpha, alpha, 10, MPFR_RNDN);  // alpha *= 0.5

            mpfr_clears(f_re_t, f_im_t, f_abs_t, (mpfr_ptr)0);
        }

        mpfr_clear(eps);
    }

    // Сохраняем результат
    mpfr_to_str(kz_str, 512, kz);
    mpfr_to_str(sz_str, 512, sz);
    mpfr_set(f_abs_out, f_abs, MPFR_RNDN);

    // Очистка
    mpfr_clears(kz, sz, kz_new, sz_new, f_re, f_im, f_abs,
                dRe_dkz, dRe_dsz, dIm_dkz, dIm_dsz, detJ,
                delta_kz, delta_sz, alpha, f_abs_new,
                kz_trial, sz_trial, (mpfr_ptr)0);

    return converged ? 1 : 0;
}

// Python-совместимая обёртка
int newton_mpfr_wrapper(
    char* kz_str,      // in/out, длина >= 512
    char* sz_str,      // in/out, длина >= 512
    char* f_abs_str,   // out, длина >= 512
    int max_iter
) {
    mpfr_t f_abs;
    mpfr_init2(f_abs, MPFR_PREC);

    int ret = newton_mpfr_512(kz_str, sz_str, max_iter, f_abs);

    if (ret >= 0) {
        mpfr_snprintf(f_abs_str, 512, "%.150Rf", f_abs);
    }

    mpfr_clear(f_abs);
    return ret;
}