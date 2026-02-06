#include <stdio.h>
#include "mendrive_det.h"
#include "mendrive_newton.h"
#include "mendrive_log.h"
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
    mendrive_scalar_t *kz, mendrive_scalar_t *sz,
    mendrive_scalar_t *step_re_re, mendrive_scalar_t *step_im_re,
    mendrive_scalar_t *step_re_im, mendrive_scalar_t *step_im_im,
    mendrive_scalar_t *f_abs_out,
    mendrive_scalar_t step_decrease,   // было step__m = 0.1 в Python
    mendrive_scalar_t step_increase,   // было step_m = 0.2 в Python (делим на него)
    mendrive_scalar_t delta_eps,
    mendrive_scalar_t f_abs_eps,
    int max_retries              // N = 10 в Python
) {
    // Вычисляем текущее значение функции
    mendrive_scalar_t f_re, f_im;
    det_eval(*kz, *sz, &f_re, &f_im);

    mendrive_scalar_t f_abs = MENDRIVE_SQRT(f_re * f_re + f_im * f_im);
    *f_abs_out = f_abs;

    MPREC_LOG_DEBUG("Начало итерации Ньютона: |f|=" MPREC_LOG_FMT_SCALAR " kz=" MPREC_LOG_FMT_SCALAR
                    ", sz=" MPREC_LOG_FMT_SCALAR "", f_abs,
                    *kz, *sz);

    // Проверка сходимости по |f|
    // if (f_abs < f_abs_eps) return 2; // converged by |f|

    // Вычисляем f/df для обоих направлений
    mendrive_scalar_t div_kz_re, div_kz_im;
    det_div_diff_kz_eval(*kz, *sz, &div_kz_re, &div_kz_im);
    MPREC_LOG_TRACE("f/df_kz = (" MPREC_LOG_FMT_SCALAR ", " MPREC_LOG_FMT_SCALAR ")",
                    div_kz_re, div_kz_im);

    // Проверка на вырожденность производных
    mendrive_scalar_t div_kz_abs = MENDRIVE_SQRT(div_kz_re * div_kz_re + div_kz_im * div_kz_im);
    if (div_kz_abs < 1e-64L) {
        MPREC_LOG_INFO("degenerate derivative div_kz_abs " MPREC_LOG_FMT_SCALAR "",
            div_kz_abs );
       return -1; // degenerate
    }

    // Сохраняем дельты для проверки сходимости
    // delta = f / df (это то, что вычитаем из координаты)

    mendrive_scalar_t f_re_new, f_im_new, f_abs_new;
    mendrive_scalar_t kz_new, sz_new;
    int n, accepted = 0;

    // ===== ШАГ 1: re_d_re - обновление kz по Re(f/df_kz) =====
    n = max_retries;
    while (n > 0) {
        n--;
        kz_new = *kz - (*step_re_re) * div_kz_re;
        det_eval(kz_new, *sz, &f_re_new, &f_im_new);
        f_abs_new = MENDRIVE_SQRT(f_re_new * f_re_new + f_im_new * f_im_new);

        MPREC_LOG_TRACE("Шаг 1: div_kz_re=" MPREC_LOG_FMT_SCALAR
                         " step" MPREC_LOG_FMT_SCALAR
                         " Δ sz=" MPREC_LOG_FMT_SCALAR,
                        div_kz_re, (*step_re_re), (*step_re_re) * div_kz_re);

        if (f_abs_new > f_abs) {
            step_decrease = f_abs / f_abs_new;
            *step_re_re *= step_decrease;
            MPREC_LOG_TRACE("Шаг 1 отклонён: |f|_new=" MPREC_LOG_FMT_SCALAR " > " MPREC_LOG_FMT_SCALAR " step_decrease=" MPREC_LOG_FMT_SCALAR " step_re_re=" MPREC_LOG_FMT_SCALAR "", 
                           f_abs_new, (f_abs), step_decrease, *step_re_re);
            continue;
        } else {
            *kz = kz_new;
            MPREC_LOG_DEBUG("Шаг 1 принят: |f|=" MPREC_LOG_FMT_SCALAR "->" MPREC_LOG_FMT_SCALAR " Δ|f|=" MPREC_LOG_FMT_SCALAR " step=" MPREC_LOG_FMT_SCALAR "kz=" MPREC_LOG_FMT_SCALAR ", ",
                            f_abs, f_abs_new, f_abs - f_abs_new, *step_re_re, *kz);
            f_abs = f_abs_new;

            *step_re_re /= step_increase;
            *step_im_re /= step_increase;
            *step_re_im /= step_increase;
            *step_im_im /= step_increase;

            accepted += 1;
            break;
        }
    }

    // ===== ШАГ 2: im_d_re - обновление kz по Im(f/df_kz) =====
    // Пересчитываем f/df_kz с новым kz
    det_div_diff_kz_eval(*kz, *sz, &div_kz_re, &div_kz_im);

    n = max_retries;
    while (n > 0) {
        n--;
        kz_new = *kz - (*step_im_re) * div_kz_im;
        det_eval(kz_new, *sz, &f_re_new, &f_im_new);
        f_abs_new = MENDRIVE_SQRT(f_re_new * f_re_new + f_im_new * f_im_new);

        MPREC_LOG_TRACE("Шаг 2: div_kz_im=" MPREC_LOG_FMT_SCALAR
                         " step" MPREC_LOG_FMT_SCALAR
                         " Δ sz=" MPREC_LOG_FMT_SCALAR,
                        div_kz_im, (*step_im_re), (*step_im_re) * div_kz_im);

        if (f_abs_new > f_abs) {
            step_decrease = f_abs / f_abs_new;
            *step_im_re *= step_decrease;
            MPREC_LOG_TRACE("Шаг 2 отклонён: |f|_new=" MPREC_LOG_FMT_SCALAR " > " MPREC_LOG_FMT_SCALAR " step_decrease=" MPREC_LOG_FMT_SCALAR " step_im_re=" MPREC_LOG_FMT_SCALAR "", 
                           f_abs_new, (f_abs), step_decrease, *step_im_re);
            continue;
        } else {
            *kz = kz_new;
            MPREC_LOG_DEBUG("Шаг 2 принят: |f|=" MPREC_LOG_FMT_SCALAR "->" MPREC_LOG_FMT_SCALAR " Δ|f|=" MPREC_LOG_FMT_SCALAR " step=" MPREC_LOG_FMT_SCALAR ", kz=" MPREC_LOG_FMT_SCALAR "",
                            f_abs, f_abs_new, f_abs - f_abs_new, *step_im_re, *kz);
            f_abs = f_abs_new;

            *step_re_re /= step_increase;
            *step_im_re /= step_increase;
            *step_re_im /= step_increase;
            *step_im_im /= step_increase;

            accepted += 1;
            break;
        }
    }

    // ===== ШАГ 3: re_d_im - обновление sz по Re(f/df_sz) =====
    // Пересчитываем f/df_sz с новым kz
    mendrive_scalar_t div_sz_re, div_sz_im;
    det_div_diff_sz_eval(*kz, *sz, &div_sz_re, &div_sz_im);
    MPREC_LOG_TRACE("f/df_sz = (" MPREC_LOG_FMT_SCALAR ", " MPREC_LOG_FMT_SCALAR ")",
                    div_sz_re, div_sz_im);

    // Проверка на вырожденность производных
    mendrive_scalar_t div_sz_abs = MENDRIVE_SQRT(div_sz_re * div_sz_re + div_sz_im * div_sz_im);
    if (div_sz_abs < 1e-64L) {
        MPREC_LOG_INFO("degenerate derivative div_sz_abs " MPREC_LOG_FMT_SCALAR "",
            div_sz_abs );
       return -1; // degenerate
    }

    n = max_retries;
    while (n > 0) {
        n--;
        sz_new = *sz - (*step_re_im) * div_sz_re;
        det_eval(*kz, sz_new, &f_re_new, &f_im_new);
        f_abs_new = MENDRIVE_SQRT(f_re_new * f_re_new + f_im_new * f_im_new);

        MPREC_LOG_TRACE("Шаг 3: div_sz_re=" MPREC_LOG_FMT_SCALAR
                         " step=" MPREC_LOG_FMT_SCALAR
                         " Δsz=" MPREC_LOG_FMT_SCALAR,
                        div_sz_re, (*step_re_im), (*step_re_im) * div_sz_re);

        if (f_abs_new > f_abs) {
            step_decrease = f_abs / f_abs_new;
            *step_re_im *= step_decrease;
            MPREC_LOG_TRACE("Шаг 3 отклонён: |f|_new=" MPREC_LOG_FMT_SCALAR " > " MPREC_LOG_FMT_SCALAR " step_decrease=" MPREC_LOG_FMT_SCALAR " step_re_im=" MPREC_LOG_FMT_SCALAR "", 
                           f_abs_new, (f_abs), step_decrease, *step_re_im);
            continue;
        } else {
            *sz = sz_new;
            MPREC_LOG_DEBUG("Шаг 3 принят: |f|=" MPREC_LOG_FMT_SCALAR "->" MPREC_LOG_FMT_SCALAR " Δ|f|=" MPREC_LOG_FMT_SCALAR " step=" MPREC_LOG_FMT_SCALAR ", sz=" MPREC_LOG_FMT_SCALAR "",
                            f_abs, f_abs_new, f_abs - f_abs_new, *step_re_im, *sz);
            f_abs = f_abs_new;

            *step_re_re /= step_increase;
            *step_im_re /= step_increase;
            *step_re_im /= step_increase;
            *step_im_im /= step_increase;

            accepted += 1;
            break;
        }
    }

    // ===== ШАГ 4: im_d_im - обновление sz по Im(f/df_sz) =====
    // Пересчитываем f/df_sz с новым sz
    det_div_diff_sz_eval(*kz, *sz, &div_sz_re, &div_sz_im);

    n = max_retries;
    while (n > 0) {
        n--;
        sz_new = *sz - (*step_im_im) * div_sz_im;
        det_eval(*kz, sz_new, &f_re_new, &f_im_new);
        f_abs_new = MENDRIVE_SQRT(f_re_new * f_re_new + f_im_new * f_im_new);

        MPREC_LOG_TRACE("Шаг 4: div_sz_im=" MPREC_LOG_FMT_SCALAR
                         " step" MPREC_LOG_FMT_SCALAR
                         " Δ sz=" MPREC_LOG_FMT_SCALAR,
                        div_sz_im, (*step_im_im), (*step_im_im) * div_sz_im);

        if (f_abs_new > f_abs) {
            step_decrease = f_abs / f_abs_new;
            *step_im_im *= step_decrease;
            MPREC_LOG_TRACE("Шаг 4 отклонён: |f|_new=" MPREC_LOG_FMT_SCALAR " > " MPREC_LOG_FMT_SCALAR " step_decrease=" MPREC_LOG_FMT_SCALAR " step_im_im=" MPREC_LOG_FMT_SCALAR "", 
                           f_abs_new, (f_abs), step_decrease, *step_im_im);
            continue;
        } else {
            *sz = sz_new;
            MPREC_LOG_DEBUG("Шаг 4 принят: |f|=" MPREC_LOG_FMT_SCALAR "->" MPREC_LOG_FMT_SCALAR " Δ|f|=" MPREC_LOG_FMT_SCALAR " step=" MPREC_LOG_FMT_SCALAR ", sz=" MPREC_LOG_FMT_SCALAR "",
                            f_abs, f_abs_new, f_abs - f_abs_new, *step_re_im, *sz);
            f_abs = f_abs_new;

            *step_re_re /= step_increase;
            *step_im_re /= step_increase;
            *step_re_im /= step_increase;
            *step_im_im /= step_increase;

            accepted += 1;
            break;
        }
    }

    // Обновляем выходное значение |f|
    det_eval(*kz, *sz, &f_re, &f_im);
    f_abs = MENDRIVE_SQRT(f_re * f_re + f_im * f_im);
    *f_abs_out = f_abs;

    // Проверка сходимости по дельтам
    if (MENDRIVE_FABS(div_kz_re) < delta_eps &&
        MENDRIVE_FABS(div_kz_im) < delta_eps &&
        MENDRIVE_FABS(div_sz_re) < delta_eps &&
        MENDRIVE_FABS(div_sz_im) < delta_eps) {
        MPREC_LOG_INFO("Сошлось по дельтам: |Δ| < " MPREC_LOG_FMT_SCALAR "", delta_eps);
        return 1; // converged by delta
    }

    // Проверка сходимости по |f| после всех шагов
    if (f_abs < f_abs_eps) {
        MPREC_LOG_INFO("Сошлось по |f|: |f| = " MPREC_LOG_FMT_SCALAR " < " MPREC_LOG_FMT_SCALAR "", f_abs, f_abs_eps);
        return 2; // converged by |f|
    }

    if (0 == accepted) {
        MPREC_LOG_INFO("no accepted steps");
        return -2;
    }

    MPREC_LOG_DEBUG("Итерация завершена: |f| = " MPREC_LOG_FMT_SCALAR "\n", f_abs);
    return 0; // continue iterations
}


// ============================================================================
// Тест метода Ньютона
// ============================================================================
int newton_adaptive(
    mendrive_scalar_t *kz,
    mendrive_scalar_t *sz
) {
    // Адаптивные шаги (как в Python-версии)
    mendrive_scalar_t step_re_re = 0.95L;
    mendrive_scalar_t step_im_re = 0.95L;
    mendrive_scalar_t step_re_im = 0.95L;
    mendrive_scalar_t step_im_im = 0.95L;

    // Параметры алгоритма
    const mendrive_scalar_t step_decrease = 0.999L;     // уменьшение шага при отказе
    const mendrive_scalar_t step_increase = 0.999999L;  // увеличение шага при успехе
    const mendrive_scalar_t delta_eps = 1e-64L;         // порог сходимости по дельте
    const mendrive_scalar_t f_abs_eps = 1e-64L;         // порог сходимости по |f|
    const int max_retries = 5;                          // макс. попыток на шаг
    const int max_iter = 50000;                         // макс. итераций

    // Вычисляем начальное значение детерминанта
    mendrive_scalar_t det_re, det_im, f_abs;
    det_eval(*kz, *sz, &det_re, &det_im);
    f_abs = MENDRIVE_COMPLEX_ABS(det_re, det_im);

    MPREC_LOG_INFO("Начальное приближение:\n");
    MPREC_LOG_INFO("  kz = " MPREC_LOG_FMT_SCALAR "\n", *kz);
    MPREC_LOG_INFO("  sz = " MPREC_LOG_FMT_SCALAR "\n", *sz);
    MPREC_LOG_INFO("  det " MPREC_LOG_FMT_SCALAR " " MPREC_LOG_FMT_SCALAR "", det_re, det_im);
    MPREC_LOG_INFO("  |det| = " MPREC_LOG_FMT_SCALAR "\n", f_abs);
    MPREC_LOG_INFO("\n");

    // Итерации Ньютона
    int iter;
    int converged = 0;
    clock_t start = clock();

    for (iter = 0; iter < max_iter; iter++) {
        mendrive_scalar_t f_abs_prev = f_abs;

        // Выполняем один адаптивный шаг Ньютона
        int ret = newton_adaptive_step(
            kz, sz,
            &step_re_re, &step_im_re, &step_re_im, &step_im_im,
            &f_abs,
            step_decrease, step_increase,
            delta_eps, f_abs_eps, max_retries
        );

        // Пересчитываем детерминант для логгирования
        det_eval(*kz, *sz, &det_re, &det_im);

        // Логгирование итерации
        MPREC_LOG_INFO("Итерация %2d: kz = " MPREC_LOG_FMT_SCALAR ", sz = " MPREC_LOG_FMT_SCALAR ", |det| = " MPREC_LOG_FMT_SCALAR ", шаги = [" MPREC_LOG_FMT_SCALAR ", " MPREC_LOG_FMT_SCALAR ", " MPREC_LOG_FMT_SCALAR ", " MPREC_LOG_FMT_SCALAR "]",
               iter, *kz, *sz, f_abs, step_re_re, step_im_re, step_re_im, step_im_im);

        if (ret == 1) {
            MPREC_LOG_INFO(" ← СОШЛОСЬ (дельта)\n");
            converged = 1;
            break;
        } else if (ret == 2) {
            MPREC_LOG_INFO(" ← СОШЛОСЬ (|f|)\n");
            converged = 1;
            break;
        } else if ( -1 == ret ) {
            MPREC_LOG_INFO(" ← ОШИБКА (производная ≈ 0)\n");
            break;
        } else {
            MPREC_LOG_INFO("\n");
        }

        mendrive_scalar_t delta_f = MENDRIVE_FABS(f_abs - f_abs_prev);
        // Защита от зацикливания
        if (delta_f < 1e-128L && iter > 5 && -2 != ret) {
            MPREC_LOG_INFO("  ⚠️  Сходимость остановилась (|Δf| " MPREC_LOG_FMT_SCALAR ")\n", delta_f);
            break;
        }
    }

    clock_t end = clock();
    double elapsed = (double)(end - start) / CLOCKS_PER_SEC * 1000.0;

    MPREC_LOG_INFO("\n");
    MPREC_LOG_INFO("\nРЕЗУЛЬТАТЫ ТЕСТА НЬЮТОНА:\n");
    MPREC_LOG_INFO("\n");
    MPREC_LOG_INFO("Статус: %s\n", converged ? "✅ СОШЁЛСЯ" : "⚠️  НЕ СОШЁЛСЯ (достигнут лимит итераций)");
    MPREC_LOG_INFO("Итераций выполнено: %d / %d\n", iter + 1, max_iter);
    MPREC_LOG_INFO("Время выполнения: %.2f мс\n", elapsed);
    MPREC_LOG_INFO("Финальное решение:\n");
    MPREC_LOG_INFO("  kz = " MPREC_LOG_FMT_SCALAR "\n", *kz);
    MPREC_LOG_INFO("  sz = " MPREC_LOG_FMT_SCALAR "\n", *sz);
    MPREC_LOG_INFO("  det " MPREC_LOG_FMT_SCALAR " " MPREC_LOG_FMT_SCALAR "", det_re, det_im);
    MPREC_LOG_INFO("  |det| = " MPREC_LOG_FMT_SCALAR "\n", f_abs);
    MPREC_LOG_INFO("\n");

    // Диагностика качества решения
    if (f_abs < 1e-10L) {
        printf("✅ Качество решения: ВЫСОКОЕ (|det| < 1e-10)\n");
    } else if (f_abs < 1e-6L) {
        printf("⚠️  Качество решения: СРЕДНЕЕ (1e-10 ≤ |det| < 1e-6)\n");
    } else {
        printf("❌ Качество решения: НИЗКОЕ (|det| ≥ 1e-6)\n");
        printf("   Возможные причины:\n");
        printf("   - Неправильный знак мнимой части волнового вектора (затухание/рост)\n");
        printf("   - Плохая обусловленность Якобиана вблизи корня\n");
        printf("   - Система уравнений физически несовместна (переопределённая 8×7)\n");
    }
    printf("\n");

    return converged ? 0 : 1;
}

