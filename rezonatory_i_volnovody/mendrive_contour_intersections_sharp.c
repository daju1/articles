#include <stdlib.h>
#include <math.h>
#include <complex.h>
#include "mendrive_point2d.h"
#include "mendrive_isolines.h"
#include "mendrive_contour_intersections_sharp.h"
#include <stdio.h> // This line includes the standard input/output library

#define USE_ADAPTIVE_SHARP
#define LOGGING 0

// Вспомогательная функция: вычисление косинуса угла в точке i
static long double compute_cosine_at_point(
    const long double* x,
    const long double* y,
    int n, int i
) {
    if (i <= 0 || i >= n - 1) return -2.0L; // недопустимое значение

    long double v1x = x[i] - x[i-1];
    long double v1y = y[i] - y[i-1];
    long double v2x = x[i+1] - x[i];
    long double v2y = y[i+1] - y[i];

    long double len1 = sqrtl(v1x*v1x + v1y*v1y);
    long double len2 = sqrtl(v2x*v2x + v2y*v2y);

    if (len1 < 1e-18L || len2 < 1e-18L) return -2.0L;

    long double dot = v1x * v2x + v1y * v2y;
    long double cos_angle = dot / (len1 * len2);

    // Ограничиваем диапазон [-1, 1]
    if (cos_angle < -1.0L) cos_angle = -1.0L;
    if (cos_angle > 1.0L) cos_angle = 1.0L;

    return cos_angle;
}

static long double compute_sin_at_point(
    const long double* x,
    const long double* y,
    int n, int i
) {
    if (i <= 0 || i >= n - 1) return -2.0L;

    long double v1x = x[i] - x[i-1];
    long double v1y = y[i] - y[i-1];
    long double v2x = x[i+1] - x[i];
    long double v2y = y[i+1] - y[i];

    long double len1 = sqrtl(v1x*v1x + v1y*v1y);
    long double len2 = sqrtl(v2x*v2x + v2y*v2y);

    if (len1 < 1e-18L || len2 < 1e-18L) return -2.0L;

    long double cross = fabsl(v1x * v2y - v1y * v2x);
    long double sin_angle = cross / (len1 * len2);

    return sin_angle;
}

// Вспомогательная функция сравнения для qsort
static int compare_long_double(const void* a, const void* b) {
    long double da = *(const long double*)a;
    long double db = *(const long double*)b;
    return (da > db) - (da < db);
}

// Адаптивное определение острых углов
void find_sharp_corners_sin(
    const long double* x, const long double* y,
    int n,
    long double sin_min_angle,
    char* is_sharp,  // выходной массив: 1 если излом, 0 иначе
    long double* sines
) {
    if (n < 3) {
        for (int i = 0; i < n; ++i) is_sharp[i] = 0;
        return;
    }

    // Вычисляем синусы
    int valid_count = 0;
    for (int i = 1; i < n - 1; ++i) {
        long double sin_val = compute_sin_at_point(x, y, n, i);
        if (sin_val >= 0.0L) { // валидное значение
            sines[valid_count++] = sin_val;
        }
    }

    if (valid_count == 0) {
        for (int i = 0; i < n; ++i) is_sharp[i] = 0;
        return;
    }

    // Сортируем
    qsort(sines, valid_count, sizeof(long double), compare_long_double);

    // Находим порог
    long double threshold = sin_min_angle; // значение по умолчанию
    int idx = (int)(valid_count * 0.9); // 90-й процентиль
    if (idx >= valid_count) idx = valid_count - 1;
    threshold = sines[idx];
    if (threshold < 0.5L) threshold = 0.5L; // минимальный порог

    // Помечаем острые углы
    for (int i = 0; i < n; ++i) is_sharp[i] = 0;
    for (int i = 1; i < n - 1; ++i) {
        long double sin_val = compute_sin_at_point(x, y, n, i);
        if (sin_val >= threshold) {
            is_sharp[i] = 1;
        }
    }
}

// Возвращает true, если точка i в ломаной является "острым" изломом
static int is_sharp_corner(
    const long double* x, const long double* y, int n, int i,
    long double cos_max_angle  // например, -0.866L для угла > 150°
) {
    long double cos_val = compute_cosine_at_point(x, y, n, i);

    // Угол "острый", если он сильно отличается от 180° (cos ≈ -1)
    // Например: если cos_angle > -0.9 → угол < ~154°, но мы хотим "излом" → лучше использовать отклонение от прямой
    // Альтернатива: смотреть на |cross| / (len1*len2) — синус угла
    // Но проще: считать изломом, если угол < threshold_deg (например, 160°)
    // Тогда: cos_angle > cos(threshold_rad)

    return (cos_val < cos_max_angle); // например, cos_max_angle = -0.94L (~160°)
}

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
) {
    for (int i = 0; i < n; ++i) is_sharp_cos[i] = 0;
    for (int i = 0; i < n; ++i) is_sharp_sin[i] = 0;

    if (n < 3) {
        return;
    }

    _Bool logging = LOGGING;

    // Шаг 1: вычисляем все косинусы
    int valid_count = 0;

    for (int i = 1; i < n - 1; ++i) {
        long double cos_val = compute_cosine_at_point(x, y, n, i);
        long double sin_val = compute_sin_at_point(x, y, n, i);
        if (cos_val > -1.5L && sin_val > -1.5L) { // валидное значение
            cosines[valid_count] = cos_val;
            sines[valid_count++] = sin_val;
            if (logging) printf("%d %Lf %Lf cos_val=%Lf sin_val=%Lf\n", i, x[i], y[i], cos_val, sin_val);
        }
    }

    if (valid_count == 0) {
        return;
    }

    // Шаг 3: ищем разрыв в распределении
    long double cos_threshold = cos_max_angle; // значение по умолчанию
    long double sin_threshold = sin_min_angle; // значение по умолчанию

#ifdef USE_ADAPTIVE_SHARP
    // Шаг 2: сортируем косинусы
    qsort(cosines, valid_count, sizeof(long double), 
          (int(*)(const void*, const void*))compare_long_double);
          
    qsort(sines, valid_count, sizeof(long double), 
          (int(*)(const void*, const void*))compare_long_double);

    const int HIST_BINS = 45;
    int cos_hist[HIST_BINS];
    int sin_hist[HIST_BINS];
    for (int i = 0; i < HIST_BINS ; ++i) {cos_hist[i] = 0;}
    for (int i = 0; i < HIST_BINS ; ++i) {sin_hist[i] = 0;}

    // Ищем первый бин с малым количеством точек после плотной области
    long double min_cos = -1.0L;
    long double max_cos = 1.0L;
    long double min_sin = -1.0L;
    long double max_sin = 1.0L;
    long double cos_bin_width = (max_cos - min_cos) / HIST_BINS;
    long double sin_bin_width = (max_sin - min_sin) / HIST_BINS;

    for (int i = 0; i < valid_count; ++i) {
        int cos_bin = (int)((cosines[i] - min_cos) / cos_bin_width);
        int sin_bin = (int)((  sines[i] - min_sin) / sin_bin_width);
        if (cos_bin < 0) cos_bin = 0;
        if (sin_bin < 0) sin_bin = 0;
        if (cos_bin >= HIST_BINS) cos_bin = HIST_BINS - 1;
        if (sin_bin >= HIST_BINS) sin_bin = HIST_BINS - 1;
        cos_hist[cos_bin]++;
        sin_hist[sin_bin]++;
    }
    
    // Находим плотную область слева (плавные повороты)
    int left_dense_end = -1;
    for (int i = 0; i < HIST_BINS; ++i) {
        if (logging) printf("%Lf cos_hist[%d]=%d\n", min_cos + i*cos_bin_width, i, cos_hist[i]);
        if (cos_hist[i] > 4) {
            break; // нашли конец плотной области
        }
        else if (cos_hist[i] > 0 && i < HIST_BINS / 2) { // порог плотности
            left_dense_end = i;
        }
    }

    if (left_dense_end != -1 && left_dense_end < HIST_BINS - 1) {
        // Берём начало разрыва как порог
        cos_threshold = min_cos + (left_dense_end + 1) * cos_bin_width;

        if (logging) printf("cos_threshold = %Lf\n", cos_threshold);
        // Но не менее -0.8 (угол 143°)
        if (cos_threshold > cos_max_angle) cos_threshold = cos_max_angle;
    }
    else {
        if (logging) printf("cos left_dense_end = %d\n", left_dense_end);
    }
    
    //////////////////////////////////////////////////////////////////////////////
    
    left_dense_end = -1;
    for (int i = 0; i < HIST_BINS; ++i) {
        if (logging) printf("%Lf sin_hist[%d]=%d\n", min_sin + i*sin_bin_width, i, sin_hist[i]);
        if (sin_hist[i] > 4) {
            break; // нашли конец плотной области
        }
        else if (sin_hist[i] > 0 && i < HIST_BINS / 2) { // порог плотности
            left_dense_end = i;
        }
    }

    if (left_dense_end != -1 && left_dense_end < HIST_BINS - 1) {
        // Берём начало разрыва как порог
        sin_threshold = min_sin + (left_dense_end + 1) * sin_bin_width;

        if (logging) printf("sin_threshold = %Lf\n", sin_threshold);
        // Но не менее -0.8 (угол 143°)
        if (sin_threshold < sin_min_angle) sin_threshold = sin_min_angle;
    }
    else {
        if (logging) printf("sin left_dense_end = %d\n", left_dense_end);
    }
    
    // Шаг 4: ищем максимальный разрыв между соседними бинами
#endif

    // Шаг 5: помечаем острые углы
    for (int i = 0; i < n; ++i) is_sharp_sin[i] = 0;
    for (int i = 0; i < n; ++i) is_sharp_cos[i] = 0;

    for (int i = 1; i < n - 1; ++i) {
        long double cos_val = compute_cosine_at_point(x, y, n, i);
        long double sin_val = compute_sin_at_point(x, y, n, i);
        cosines[i-1] = cos_val;
        sines[i-1]   = sin_val;
        if (cos_val <= cos_threshold) {
            is_sharp_cos[i] = 1;
        }
        if (sin_val > sin_threshold && sin_val < sin_max_angle) {
            is_sharp_sin[i] = 1;
        }
    }
}

// Экстраполирует одно плечо угла назад или вперёд
// dir = -1 → экстраполяция "назад" (через [i-1] -> [i])
// dir = +1 → экстраполяция "вперёд" (через [i] -> [i+1])
static void extrapolate_segment(
    const long double* x, const long double* y, int i, int dir,
    long double extrap_len,
    long double* out_p_x, long double* out_p_y,
    long double* out_r_x, long double* out_r_y
) {
    long double dx, dy, len;
    if (dir == -1) {
        dx = x[i] - x[i-1];
        dy = y[i] - y[i-1];
    } else {
        dx = x[i+1] - x[i];
        dy = y[i+1] - y[i];
    }
    len = sqrtl(dx*dx + dy*dy);
    if (len < 1e-18L) {
        dx = dy = 0.0L;
        len = 1.0L;
    }

    // Единичный вектор направления
    long double ux = dx / len;
    long double uy = dy / len;

    // Точка начала экстраполированного отрезка — это угловая точка
    *out_p_x = x[i];
    *out_p_y = y[i];

    // Направление — продолжение плеча
    *out_r_x = ux * extrap_len;
    *out_r_y = uy * extrap_len;
}

// Вспомогательная: пересечение двух отрезков
// Возвращает 1, если пересекаются, и заполняет (out_x, out_y)
// Использует параметрическое представление:
//   p + t * r,    q + u * s,    t,u ∈ [0,1]
static int segment_intersection(
    long double p_x, long double p_y,
    long double r_x, long double r_y,
    long double q_x, long double q_y,
    long double s_x, long double s_y,
    long double* out_x, long double* out_y
) {
    long double r_cross_s = r_x * s_y - r_y * s_x;
    if ( fabsl(r_cross_s ) < 1e-18L) {
        // printf("parallel\n");
        return 0; // параллельны или совпадают
    }

    long double q_minus_p_x = q_x - p_x;
    long double q_minus_p_y = q_y - p_y;

    long double t = (q_minus_p_x * s_y - q_minus_p_y * s_x) / r_cross_s;
    long double u = (q_minus_p_x * r_y - q_minus_p_y * r_x) / r_cross_s;

    if (t < 0.0L || t > 1.0L || u < 0.0L || u > 1.0L)
    //if (t < -1.0L || t > 2.0L || u < -1.0L || u > 2.0L)
    {
        // printf("t = %Lf, u = %Lf\n", t, u);
        return 0;
    }

    *out_x = p_x + t * r_x;
    *out_y = p_y + t * r_y;
    return 1;
}

// Возвращаем 2, если отрезки почти параллельны и близки — возможно, касание
/*static int segment_intersection_or_tangent(
    long double p_x, long double p_y,
    long double r_x, long double r_y,
    long double q_x, long double q_y,
    long double s_x, long double s_y,
    long double* out_x, long double* out_y,
    long double eps_det,
    long double eps_dist   // например, 1e-4L — порог близости
) {
    long double r_cross_s = r_x * s_y - r_y * s_x;

    if (fabsl(r_cross_s) < eps_det) {
        // Параллельны — проверим расстояние от точки q до прямой p + t*r

        // Расстояние от точки q до прямой через p с направлением r
        long double dist_num = fabsl( (q_x - p_x)*r_y - (q_y - p_y)*r_x );
        long double len_r = sqrtl(r_x*r_x + r_y*r_y);
        if (len_r < eps_det) return 0;

        long double dist = dist_num / len_r;

        if (dist > eps_dist) return 0; // слишком далеко — не касание

        // Проверим, находится ли проекция точки q на прямую p + t*r в пределах отрезка [0,1]
        long double proj_t = ((q_x - p_x)*r_x + (q_y - p_y)*r_y) / (r_x*r_x + r_y*r_y);
        if (proj_t < 0.0L || proj_t > 1.0L) return 0;

        // Точка касания — проекция q на прямую p + t*r
        *out_x = p_x + proj_t * r_x;
        *out_y = p_y + proj_t * r_y;

        return 2; // "касание"
    }

    // Обычное пересечение
    long double q_minus_p_x = q_x - p_x;
    long double q_minus_p_y = q_y - p_y;

    long double t = (q_minus_p_x * s_y - q_minus_p_y * s_x) / r_cross_s;
    long double u = (q_minus_p_x * r_y - q_minus_p_y * r_x) / r_cross_s;

    if (t < 0.0L || t > 1.0L || u < 0.0L || u > 1.0L)
        return 0;

    *out_x = p_x + t * r_x;
    *out_y = p_y + t * r_y;
    return 1;
}*/

// Возвращает длину самого длинного гладкого участка, центрированного в точке center_idx
// Гладкость определяется по разнице между соседними значениями (cosines или sines)
static int smooth_segment_length(
    const long double* values,   // массив значений (cosines или sines)
    int n_vals,                  // размер массива (обычно n-2)
    int center_idx,              // индекс в массиве values (соответствует точке i+1 в исходных данных)
    long double eps_smooth       // порог изменения, например 0.05L
) {
    if (center_idx < 0 || center_idx >= n_vals) return 0;

    int max_len = 1; // минимум — сама точка

    // Ищем влево
    int left = center_idx;
    while (left > 0) {
        long double diff = fabsl(values[left] - values[left - 1]);
        if (diff > eps_smooth) break;
        left--;
    }

    // Ищем вправо
    int right = center_idx;
    while (right < n_vals - 1) {
        long double diff = fabsl(values[right + 1] - values[right]);
        if (diff > eps_smooth) break;
        right++;
    }

    return right - left + 1;
}


// Уточняет позицию острой вершины, выбирая ту из трёх (i-1, i, i+1),
// которая лежит на самом длинном гладком участке (по косинусам и/или синусам)
// Возвращает индекс лучшей точки (в пределах [1, n-2])
static int refine_sharp_corner_by_smoothness(
    const long double* cosines, int n_cos,
    const long double* sines,   int n_sin,
    int i,                      // текущая формально острая вершина (индекс в исходной ломаной)
    int n_points,               // общее число точек ломаной
    long double eps_smooth      // порог гладкости, например 0.05L
) {
    if (n_points < 5 || i < 2 || i >= n_points - 2) {
        return i; // недостаточно точек для анализа
    }

    int candidates[3] = {i-1, i, i+1};
    int smooth_lengths[3] = {0, 0, 0};

    for (int k = 0; k < 3; ++k) {
        int cand = candidates[k];
        if (cand < 1 || cand >= n_points - 1) {
            smooth_lengths[k] = 0;
            continue;
        }

        // индекс в массивах cosines/sines: cand - 1
        int cos_idx = cand - 1;
        int sin_idx = cand - 1;

        if (cos_idx < 0 || cos_idx >= n_cos) {
            smooth_lengths[k] = 0;
            continue;
        }
        if (sin_idx < 0 || sin_idx >= n_sin) {
            smooth_lengths[k] = 0;
            continue;
        }

        // Вычисляем длину гладкого участка по косинусам
        int len_cos = smooth_segment_length(cosines, n_cos, cos_idx, eps_smooth);

        // Вычисляем длину гладкого участка по синусам
        int len_sin = smooth_segment_length(  sines, n_sin, sin_idx, eps_smooth);

        // Берём максимум из двух — можно также усреднять или взвешивать
        smooth_lengths[k] = len_cos > len_sin ? len_cos : len_sin;
    }

    // Находим кандидата с максимальной длиной гладкого участка
    int best_k = 0;
    for (int k = 1; k < 3; ++k) {
        if (smooth_lengths[k] > smooth_lengths[best_k]) {
            best_k = k;
        }
    }

    return candidates[best_k];
}

int find_contour_intersections_with_corners(
    const long double* cu_x, const long double* cu_y, int cu_n,
    const long double* cv_x, const long double* cv_y, int cv_n,
    point2d_t* intersections, int max_intersections,
    long double eps_det,
    long double extrap_len,           // длина экстраполяции (например, 1e-3)
    long double cos_max_angle,        // порог остроты угла, например -0.94L
    long double sin_min_angle,        // порог остроты угла, например 0.5L
    long double sin_max_angle         // порог остроты угла, например 0.8L
) {
    if (cu_n < 2 || cv_n < 2 || !intersections || max_intersections <= 0) {
        return 0;
    }

    int count = 0;

    // 1. Обычные пересечения отрезок-отрезок
    for (int i = 0; i < cu_n - 1; ++i) {
        long double p_x = cu_x[i], p_y = cu_y[i];
        long double r_x = cu_x[i+1] - p_x;
        long double r_y = cu_y[i+1] - p_y;

        for (int j = 0; j < cv_n - 1; ++j) {
            long double q_x = cv_x[j], q_y = cv_y[j];
            long double s_x = cv_x[j+1] - q_x;
            long double s_y = cv_y[j+1] - q_y;

            long double x, y;
 #if 0
            int res = segment_intersection_or_tangent(
                p_x, p_y, r_x, r_y,
                q_x, q_y, s_x, s_y,
                &x, &y,
                eps_det,
                1e-4L  // eps_dist — подберите под масштаб
            );
 #else
             int res = segment_intersection(
                p_x, p_y, r_x, r_y,
                q_x, q_y, s_x, s_y,
                &x, &y
            );
#endif

            /*if (res == 1 || res == 2) {  // 1 — пересечение, 2 — касание
                if (count >= max_intersections) goto overflow;
                intersections[count].kz = x;
                intersections[count].sz = y;
                count++;
            }*/
        }
    }

    // Вместо цикла с is_sharp_corner
    char* sharp_mask_cos = (char*)calloc(cu_n, sizeof(char));
    char* sharp_mask_sin = (char*)calloc(cu_n, sizeof(char));
    long double* cosines = (long double*)malloc((cu_n - 2) * sizeof(long double));
    long double* sines = (long double*)malloc((cu_n - 2) * sizeof(long double));

    find_sharp_corners(cu_x, cu_y, cu_n,
                       cos_max_angle,
                       sin_min_angle,
                       sin_max_angle,
                       sharp_mask_cos,
                       sharp_mask_sin,
                       cosines, sines);

    // 2. Обработка острых изломов на cu → экстраполяция → пересечение с cv
    for (int i = 1; i < cu_n - 1; ++i) {
        if (!sharp_mask_cos[i] && !sharp_mask_sin[i]) continue;

        // Уточняем позицию излома по гладкости
        int refined_i = refine_sharp_corner_by_smoothness(
            cosines, cu_n - 2,
            sines,   cu_n - 2,
            i, cu_n,
            0.05L  // eps_smooth — подберите под вашу задачу
        );

        // Проверяем, что refined_i валиден
        if (refined_i < 1 || refined_i >= cu_n - 1) {
            refined_i = i; // fallback
        }

        // Теперь работаем с refined_i
        for (int dir = -1; dir <= 1; dir += 2) {
            long double p_x, p_y, r_x, r_y;

            extrapolate_segment(cu_x, cu_y, refined_i, dir, extrap_len, &p_x, &p_y, &r_x, &r_y);

        // Два плеча: назад и вперёд
        //for (int dir = -1; dir <= 1; dir += 2) {
        //    long double p_x, p_y, r_x, r_y;
        //    extrapolate_segment(cu_x, cu_y, i, dir, extrap_len, &p_x, &p_y, &r_x, &r_y);

            // Пересекаем этот экстраполированный отрезок со всеми отрезками cv
            for (int j = 0; j < cv_n - 1; ++j) {
                long double q_x = cv_x[j], q_y = cv_y[j];
                long double s_x = cv_x[j+1] - q_x;
                long double s_y = cv_y[j+1] - q_y;

                long double x, y;
                // Важно: используем eps_det!
                long double r_cross_s = r_x * s_y - r_y * s_x;
                if (fabsl(r_cross_s) < eps_det) continue;

                long double qm_px = q_x - p_x;
                long double qm_py = q_y - p_y;
                long double t = (qm_px * s_y - qm_py * s_x) / r_cross_s;
                long double u = (qm_px * r_y - qm_py * r_x) / r_cross_s;

                // Экстраполированный отрезок: t ∈ [0, 1] (мы его задали длиной extrap_len)
                // Отрезок cv: u ∈ [0, 1]
                if (t >= 0.0L && t <= 1.0L && u >= 0.0L && u <= 1.0L) {
                //if (t >= -0.5L && t <= 1.5L && u >= -0.5L && u <= 1.5L) {
                    x = p_x + t * r_x;
                    y = p_y + t * r_y;
                    if (count >= max_intersections)
                    {
                        free(sines);
                        free(cosines);
                        free(sharp_mask_cos);
                        free(sharp_mask_sin);
                        goto overflow;
                    }
                    intersections[count].kz = x;
                    intersections[count].sz = y;
                    count++;
                }
            }
        }
    }

    free(sines);
    free(cosines);
    free(sharp_mask_cos);
    free(sharp_mask_sin);

    // Вместо цикла с is_sharp_corner
    sharp_mask_cos = (char*)calloc(cv_n, sizeof(char));
    sharp_mask_sin = (char*)calloc(cv_n, sizeof(char));
    cosines = (long double*)malloc((cv_n - 2) * sizeof(long double));
    sines = (long double*)malloc((cv_n - 2) * sizeof(long double));

    find_sharp_corners(cv_x, cv_y, cv_n,
                       cos_max_angle,
                       sin_min_angle,
                       sin_max_angle,
                       sharp_mask_cos,
                       sharp_mask_sin,
                       cosines, sines);

    // 3. Обработка острых изломов на cv → экстраполяция → пересечение с cu
    for (int j = 1; j < cv_n - 1; ++j) {
        if (!sharp_mask_cos[j] && !sharp_mask_sin[j]) continue;

        // Уточняем позицию излома по гладкости
        int refined_j = refine_sharp_corner_by_smoothness(
            cosines, cv_n - 2,
            sines,   cv_n - 2,
            j, cv_n,
            0.05L  // eps_smooth
        );

        if (refined_j < 1 || refined_j >= cv_n - 1) {
            refined_j = j;
        }

        for (int dir = -1; dir <= 1; dir += 2) {
            long double p_x, p_y, r_x, r_y;
            extrapolate_segment(cv_x, cv_y, refined_j, dir, extrap_len, &p_x, &p_y, &r_x, &r_y);

        //for (int dir = -1; dir <= 1; dir += 2) {
        //    long double p_x, p_y, r_x, r_y;
        //    extrapolate_segment(cv_x, cv_y, j, dir, extrap_len, &p_x, &p_y, &r_x, &r_y);

            for (int i = 0; i < cu_n - 1; ++i) {
                long double q_x = cu_x[i], q_y = cu_y[i];
                long double s_x = cu_x[i+1] - q_x;
                long double s_y = cu_y[i+1] - q_y;

                long double x, y;
                long double r_cross_s = r_x * s_y - r_y * s_x;
                if (fabsl(r_cross_s) < eps_det) continue;

                long double qm_px = q_x - p_x;
                long double qm_py = q_y - p_y;
                long double t = (qm_px * s_y - qm_py * s_x) / r_cross_s;
                long double u = (qm_px * r_y - qm_py * r_x) / r_cross_s;

                if (t >= 0.0L && t <= 1.0L && u >= 0.0L && u <= 1.0L) {
                //if (t >= -0.5L && t <= 1.5L && u >= -0.5L && u <= 1.5L) {
                    x = p_x + t * r_x;
                    y = p_y + t * r_y;
                    if (count >= max_intersections)
                    {
                        free(sines);
                        free(cosines);
                        free(sharp_mask_cos);
                        free(sharp_mask_sin);
                        goto overflow;
                    }
                    intersections[count].kz = x;
                    intersections[count].sz = y;
                    count++;
                }
            }
        }
    }
    free(sines);
    free(cosines);
    free(sharp_mask_cos);
    free(sharp_mask_sin);

    return count;

overflow:
    fprintf(stderr, "find_contour_intersections_with_corners: max_intersections exceeded\n");
    return count;
}

// Тестовая функция для острых углов
int test_sharp_corners(const contour_line_t* line,
                       long double cos_max_angle,
                       long double sin_min_angle,
                       long double sin_max_angle,
                       const char* name,
                       corner2d_t* sharp_corners,
                       int max_sharp_corners) {
    if (!line || line->n_points < 3) return 0;

    _Bool logging = LOGGING;

    // Выделяем память для маски
    char* sharp_mask_cos = (char*)calloc(line->n_points, sizeof(char));
    char* sharp_mask_sin = (char*)calloc(line->n_points, sizeof(char));
    long double* cosines = (long double*)malloc((line->n_points - 2) * sizeof(long double));
    long double* sines = (long double*)malloc((line->n_points - 2) * sizeof(long double));
    if (!sharp_mask_cos || !sharp_mask_sin ) {
        free(sharp_mask_cos);
        free(sharp_mask_sin);
        return 0;
    }
    if (!cosines)
    {
        free(sharp_mask_cos);
        free(sharp_mask_sin);
        return 0;
    }
    if (!sines)
    {
        free(cosines);
        free(sharp_mask_cos);
        free(sharp_mask_sin);
        return 0;
    }
    // Извлекаем координаты
    long double* x = (long double*)malloc(line->n_points * sizeof(long double));
    long double* y = (long double*)malloc(line->n_points * sizeof(long double));
    if (!x || !y) {
        free(sharp_mask_cos);
        free(sharp_mask_sin);
        free(cosines); free(x); free(y);
        return 0;
    }

    for (int i = 0; i < line->n_points; ++i) {
        x[i] = line->points[i].kz;
        y[i] = line->points[i].sz;
    }

    find_sharp_corners(x, y, line->n_points,
                       cos_max_angle,
                       sin_min_angle,
                       sin_max_angle,
                       sharp_mask_cos,
                       sharp_mask_sin,
                       cosines, sines);

    // Выводим результаты

    if (logging) printf("\n=== Острые углы в %s n_points=%d === \n", name, line->n_points);
    int sharp_count = 0;
    for (int i = 1; i < line->n_points - 1; ++i) {
        if (sharp_mask_cos[i] || sharp_mask_sin[i]) {
            if (sharp_count >= max_sharp_corners)
            {
                free(sharp_mask_cos);
                free(sharp_mask_sin);
                free(cosines);
                free(sines);
                free(x);
                free(y);
                goto overflow;
            }

            // Уточняем позицию излома по гладкости
            int refined_i = refine_sharp_corner_by_smoothness(
                cosines, line->n_points - 2,
                sines,   line->n_points - 2,
                i, line->n_points,
                0.05L  // eps_smooth — подберите под вашу задачу
            );

            // Проверяем, что refined_i валиден
            if (refined_i < 1 || refined_i >= line->n_points - 1) {
                refined_i = i; // fallback
            }

            // Теперь работаем с refined_i
            //for (int dir = -1; dir <= 1; dir += 2) {
            //    long double p_x, p_y, r_x, r_y;
            //    extrapolate_segment(cu_x, cu_y, refined_i, dir, extrap_len, &p_x, &p_y, &r_x, &r_y);

            if (logging)
                printf("  Точка %d: (%.6Lf, %.6Lf) cos = %Lf sin = %Lf mask_cos=%d mask_sin=%d refined_i=%d refined_pt: (%.6Lf, %.6Lf)\n",
                       i, x[i], y[i],
                       cosines[i-1], sines[i-1],
                       sharp_mask_cos[i], sharp_mask_sin[i],
                       refined_i, x[refined_i], y[refined_i]);

            sharp_corners[sharp_count].kz = x[i];
            sharp_corners[sharp_count].sz = y[i];
            sharp_corners[sharp_count].cosine = cosines[i-1];
            sharp_corners[sharp_count].sine   = sines[i-1];
            if (sharp_mask_cos[i] && sharp_mask_sin[i])
                sharp_corners[sharp_count].type = 1;
            else if (sharp_mask_cos[i])
                sharp_corners[sharp_count].type = 2;
            else if (sharp_mask_sin[i])
                sharp_corners[sharp_count].type = 3;
            sharp_corners[sharp_count].i = i;
            sharp_corners[sharp_count].refined_i;
            sharp_corners[sharp_count].refined_pt.kz = x[refined_i];
            sharp_corners[sharp_count].refined_pt.sz = y[refined_i];

            sharp_count++;
        }
    }
    if (logging) printf("Найдено острых углов: %d\n", sharp_count);

    free(sharp_mask_cos);
    free(sharp_mask_sin);
    free(cosines);
    free(sines);
    free(x);
    free(y);
    return sharp_count;

overflow:
    fprintf(stderr, "test_sharp_corners: max_sharp_corners exceeded\n");
    return sharp_count;
}
