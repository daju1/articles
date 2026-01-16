#ifndef MENDRIVE_NEWTON_H
#define MENDRIVE_NEWTON_H

/**
 * Внешние функции, которые должны быть предоставлены (генерируются из SageMath):
 * - det_eval: вычисление детерминанта f(kz, sz)
 * - det_div_diff_kz_eval: вычисление f(kz,sz) / (df/dkz)
 * - det_div_diff_sz_eval: вычисление f(kz,sz) / (df/dsz)
 */
void det_eval(long double kz, long double sz,
              long double *re, long double *im);

void det_div_diff_kz_eval(long double kz, long double sz,
                          long double *re, long double *im);

void det_div_diff_sz_eval(long double kz, long double sz,
                          long double *re, long double *im);

/**
 * newton_adaptive_step - один шаг адаптивного метода Ньютона
 *
 * Параметры:
 *   kz, sz          - указатели на текущие координаты (обновляются)
 *   step_re_re      - адаптивный шаг для Re(f/df_kz) -> kz
 *   step_im_re      - адаптивный шаг для Im(f/df_kz) -> kz
 *   step_re_im      - адаптивный шаг для Re(f/df_sz) -> sz
 *   step_im_im      - адаптивный шаг для Im(f/df_sz) -> sz
 *   f_abs_out       - выходное значение |f| после шага
 *   abs_m           - множитель для проверки улучшения (обычно 1.0)
 *   step_decrease   - коэффициент уменьшения шага при неудаче (обычно 0.1)
 *   step_increase   - коэффициент увеличения шага (делим на него, обычно 0.2)
 *   delta_eps       - порог сходимости по дельтам
 *   f_abs_eps       - порог сходимости по |f|
 *   max_retries     - максимум повторов на каждом подшаге (обычно 10)
 *
 * Возвращает:
 *   -1 : ошибка (вырожденная производная)
 *    0 : продолжать итерации
 *    1 : сходимость по дельтам
 *    2 : сходимость по |f|
 */
int newton_adaptive_step(
    long double *kz, long double *sz,
    long double *step_re_re, long double *step_im_re,
    long double *step_re_im, long double *step_im_im,
    long double *f_abs_out,
    long double abs_m,
    long double step_decrease,
    long double step_increase,
    long double delta_eps,
    long double f_abs_eps,
    int max_retries
);

#endif /* MENDRIVE_NEWTON_H */
