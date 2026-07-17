from sage.all import RealField, ComplexField

from .mendrivelibc.common import get_sign_K_fast

# Create a SageMath RealField with high precision
complex128 = ComplexField(128)
real128    = RealField(128)

import numpy as np

N_re = 1000
N_im = 1000


use_phase_y = False


def calc_extrap_len(kz_range, sz_range, nk, ns):
    extrap_len = np.sqrt(
        ((kz_range[1] - kz_range[0])/nk)**2 +
        ((sz_range[1] - sz_range[0])/ns)**2)
    return extrap_len

# deprecated
def init_det_plot2d(kz_range, sz_range, nk, ns):
    kz_min, kz_max = kz_range
    sz_min, sz_max = sz_range

    kz_linspace = np.linspace(real128(kz_min), real128(kz_max), nk)
    sz_linspace = np.linspace(real128(sz_min), real128(sz_max), ns)

    kz_list = kz_linspace.tolist()
    sz_list = sz_linspace.tolist()

    kz_grid, sz_grid = np.meshgrid(kz_linspace, sz_linspace)

    u = kz_grid * np.nan
    v = sz_grid * np.nan

    shl = kz_grid * np.nan
    sel = sz_grid * np.nan
    shr = kz_grid * np.nan
    ser = sz_grid * np.nan

    return kz_grid, sz_grid, kz_list, sz_list, u, v, shl, sel, shr, ser

def plotDetContoursResult(contours):
    from sage.plot.all import Graphics
    from sage.plot.all import list_plot

    # === Извлечение ломаных ===
    cu = []
    cv = []
    pl = Graphics()

    for i_contour in range(contours.n_re_contours):
        line = contours.re_zero[i_contour]
        seg_u = []
        for i in range(line.n_points):
            pt = line.points[i]
            seg_u.append((float(pt.re), float(pt.im)))
        cu.append(seg_u)
        pl += list_plot(seg_u, size=5)

    for i_contour in range(contours.n_im_contours):
        line = contours.im_zero[i_contour]
        seg_v = []
        for i in range(line.n_points):
            pt = line.points[i]
            seg_v.append((float(pt.re), float(pt.im)))
        cv.append(seg_v)
        pl += list_plot(seg_v, color="red", size=4)

    return cu, cv, pl

def find_matplotlib_contour_intersections_c(lib, seg_u, seg_v,
                                            eps_det=1e-6,
                                            eps_sin_r_s=1e-3,
                                            extrap_len = 1e-3,
                                            cos_max_angle = 0.3,
                                            max_n=100):
    """
    cu, cv — объекты matplotlib.contour.QuadContourSet
    cu_idx, cv_idx — индексы нулевых уровней (например, cu_index_of_zero_level)
    """
    import numpy as np
    from ctypes import Structure, c_longdouble, POINTER, c_int

    class Point2D(Structure):
        _fields_ = [("x", c_longdouble),
                    ("y", c_longdouble),
                    ("sin_r_s", c_longdouble)]

    intersections = (Point2D * int(max_n))()

    # Берём первую замкнутую линию на каждом уровне (если их несколько — расширьте цикл)
    #seg_u = cu.allsegs[cu_idx][0]  # (N, 2)
    #seg_v = cv.allsegs[cv_idx][0]  # (M, 2)
    #seg_u = cu.allsegs[cu_idx][i]
    #seg_v = cv.allsegs[cv_idx][j]

    cu_x = np.ascontiguousarray(seg_u[:, 0], dtype=np.longdouble)
    cu_y = np.ascontiguousarray(seg_u[:, 1], dtype=np.longdouble)
    cv_x = np.ascontiguousarray(seg_v[:, 0], dtype=np.longdouble)
    cv_y = np.ascontiguousarray(seg_v[:, 1], dtype=np.longdouble)

    # Объявляем
    lib.find_contour_intersections.argtypes = [
        POINTER(c_longdouble), POINTER(c_longdouble), c_int,
        POINTER(c_longdouble), POINTER(c_longdouble), c_int,
        POINTER(Point2D),
        c_int,              # int max_intersections,
        c_longdouble,       # long double eps_det,
        c_longdouble,       # long double eps_sin_r_s
    ]
    lib.find_contour_intersections.restype = c_int

    n_found = lib.find_contour_intersections(
        cu_x.ctypes.data_as(POINTER(c_longdouble)),
        cu_y.ctypes.data_as(POINTER(c_longdouble)), len(cu_x),
        cv_x.ctypes.data_as(POINTER(c_longdouble)),
        cv_y.ctypes.data_as(POINTER(c_longdouble)), len(cv_y),
        intersections, max_n,
        c_longdouble(eps_det),
        c_longdouble(eps_sin_r_s)
    )

    result = []
    for i in range(n_found):
        pt = intersections[i]
        result.append((float(pt.re), float(pt.im)))
    return result

def find_resonance_roots_fast(lib, omega_re_range, omega_im_range, n_re=100, n_im=100, make_plot=False, debug_plot=False):
    """
    Быстрый поиск резонансных корней с использованием C-библиотеки.

    Returns:
        Список кортежей (kz, sz) - приближённые положения корней
    """
    import numpy as np
    from .mendrivelibc.common import det_eval_fast

    omega_re_min, omega_re_max = omega_re_range
    omega_im_min, omega_im_max = omega_im_range

    omegare_linspace = np.linspace(float(omega_re_min), float(omega_re_max), n_re)
    gamma_linspace   = np.linspace(float(omega_im_min), float(omega_im_max), n_im)

    omegare_list = omegare_linspace.tolist()
    gamma_list   = gamma_linspace.tolist()

    omegare_grid, gamma_grid = np.meshgrid(omegare_linspace, gamma_linspace)

    u = omegare_grid * np.nan
    v = gamma_grid * np.nan

    for i, omegare_val in enumerate(omegare_linspace):
        for j, gamma_val in enumerate(gamma_linspace):
            re, im = det_eval_fast(lib, omegare_val, gamma_val)
            # print(i, j, omegare_val, gamma_val, re, im)
            u[j, i] = float(re)
            v[j, i] = float(im)

    # Находим пересечения контуров Re=0 и Im=0
    import matplotlib.pyplot as plt
    fig, ax = plt.subplots()
    cu = ax.contour(omegare_linspace, gamma_linspace, u, levels=[0])
    cv = ax.contour(omegare_linspace, gamma_linspace, v, levels=[0])
    if True:
        plt.close(fig)
    else:
        ax.set_title('M8_subs_det_re')
        ax.set_xlabel('kz')
        ax.set_ylabel('sz')
        plt.show()

    # Извлекаем точки пересечения
    k_z_graphic_solutions_ = []
    if len(cu.allsegs) > 0 and len(cv.allsegs) > 0:
        for seg_u in cu.allsegs[0]:
            for seg_v in cv.allsegs[0]:
                # Находим ближайшие точки между контурами

                if debug_plot:
                    from sage.plot.all import Graphics
                    from sage.plot.all import list_plot
                    pl = Graphics()
                    pl += list_plot(seg_u, size = 5)
                    pl += list_plot(seg_v, color="red", size = 1)
                    pl.show()

                result = find_matplotlib_contour_intersections_c(lib, seg_u, seg_v)

                if None == result:
                    continue

                for found in result:
                    k_z_graphic_solutions_ += [found]

    if make_plot:
        # Define your lists of x and y coordinates
        x_coords = [k_z[0] for k_z in k_z_graphic_solutions_]
        y_coords = [k_z[0] for k_z in k_z_graphic_solutions_]

        # Create a scatter plot
        plt.scatter(x_coords, y_coords)
        plt.show()

    return k_z_graphic_solutions_

def find_resonance_roots_fast_с(lib, re_range, im_range,
                                n_re=N_re, n_im=N_im,
                                make_plot=True,
                                use_tracing=0,
                                min_isoline_points_count=0,
                                isoline_merge_segments_epsilon=1e-10,
                                sharp = 0,
                                eps_sin_r_s = 1e-3,
                                cos_max_angle = -0.94,
                                window_size = 5,
                                local_angle_staircase_threshold = 0.3,
                                total_angle_threshold           = 0.3,
                                concentration_threshold         = 0.4,
                                local_angle_sharp_threshold     = 0.6,
                                det_threshold = 1.0
                               ):
    import ctypes
    from ctypes import Structure, c_longdouble, POINTER, c_int
    from sage.plot.all import list_plot

    class Point2D(Structure):
        _fields_ = [("x", c_longdouble),
                    ("y", c_longdouble),
                    ("sin_r_s", c_longdouble)
                   ]

    # Определите структуру
    class CharacteristicRoots(Structure):
        _fields_ = [("roots", POINTER(Point2D)),
                    ("n_roots", c_int)]

    class ContourLine(Structure):
        _fields_ = [("points", POINTER(Point2D)),
                    ("n_points", c_int)]

    class DetContoursResult(Structure):
        _fields_ = [
            ("re_zero", POINTER(ContourLine)),
            ("n_re_contours", c_int),
            ("im_zero", POINTER(ContourLine)),
            ("n_im_contours", c_int)
        ]

    # Настройка
    lib.find_characteristic_roots.argtypes = [
        ctypes.c_longdouble, ctypes.c_longdouble, ctypes.c_int,
        ctypes.c_longdouble, ctypes.c_longdouble, ctypes.c_int,
        ctypes.POINTER(DetContoursResult),
        ctypes.POINTER(CharacteristicRoots),
        ctypes.c_longdouble,  # eps_nan
        ctypes.c_longdouble,  # eps_det
        ctypes.c_longdouble,  # eps_sin_r_s
        ctypes.c_longdouble,  # extrap_len
        ctypes.c_longdouble,  # cos_max_angle
        ctypes.c_longdouble,  # sin_min_angle
        ctypes.c_longdouble,  # sin_max_angle
        ctypes.c_int,         # window_size
        ctypes.c_int,         # use_tracing
        ctypes.c_int,         # min_isoline_points_count
        ctypes.c_longdouble,  # isoline_merge_segments_epsilon
        ctypes.c_int,         # sharp
        ctypes.c_longdouble,  # local_angle_staircase_threshold,
        ctypes.c_longdouble,  # total_angle_threshold
        ctypes.c_longdouble,  # concentration_threshold
        ctypes.c_longdouble,  # local_angle_sharp_threshold
        ctypes.c_longdouble,  # det_threshold
    ]
    lib.find_characteristic_roots.restype = ctypes.c_int

    lib.free_det_contours.argtypes = [ POINTER(DetContoursResult) ]
    lib.free_characteristic_roots.argtypes = [ POINTER(CharacteristicRoots) ]

    extrap_len = calc_extrap_len(re_range, im_range, n_re, n_im)

    # print("extrap_len", extrap_len)

    # Вызов
    roots = CharacteristicRoots()
    contours = DetContoursResult()
    ret = lib.find_characteristic_roots(
        ctypes.c_longdouble(re_range[0]), ctypes.c_longdouble(re_range[1]), ctypes.c_int(n_re),
        ctypes.c_longdouble(im_range[0]), ctypes.c_longdouble(im_range[1]), ctypes.c_int(n_im),
        ctypes.byref(contours),
        ctypes.byref(roots),
        ctypes.c_longdouble(1e300),    # eps_nan
        ctypes.c_longdouble(1e-12),    # eps_det
        ctypes.c_longdouble(eps_sin_r_s),
        ctypes.c_longdouble(extrap_len),     # extrap_len
        ctypes.c_longdouble(cos_max_angle),    # cos_max_angle (~160°)
        ctypes.c_longdouble(+0.5),     # sin_min_angle (~160°)
        ctypes.c_longdouble(+0.8),     # sin_max_angle (~160°)
        ctypes.c_int(window_size),
        ctypes.c_int(use_tracing), # use_tracing = 0 → Marching Squares
        ctypes.c_int(min_isoline_points_count),
        ctypes.c_longdouble(isoline_merge_segments_epsilon),
        ctypes.c_int(sharp),
        ctypes.c_longdouble(local_angle_staircase_threshold),
        ctypes.c_longdouble(total_angle_threshold),
        ctypes.c_longdouble(concentration_threshold),
        ctypes.c_longdouble(local_angle_sharp_threshold),
        ctypes.c_longdouble(det_threshold),
    )

    graphic_solutions = []

    print(f"Найдено корней: {roots.n_roots}")
    for i in range(roots.n_roots):
        print(f"  re={roots.roots[i].re:.6f}, im={roots.roots[i].im:.6f}, sin_r_s={roots.roots[i].sin_r_s:.6e}")
        graphic_solutions.extend([(roots.roots[i].re,
                                       roots.roots[i].im
                                      )])

    cu, cv, pl = plotDetContoursResult(contours)
    pl += list_plot(graphic_solutions, color="green", size=32)
    if make_plot:
        pl.show()
    lib.free_det_contours(contours)
    lib.free_characteristic_roots(roots)
    return graphic_solutions, (cu, cv, pl)

def prec_graphic_solutions(lib, solver, preciser, sol_num, graphic_solutions, precise,
                           verbose=True, return_all_stages=False):
    """
    Уточняет графическое решение методом Ньютона.

    Args:
        return_all_stages: если True, возвращает результаты всех стадий уточнения

    Returns:
        dict с информацией о корне или None при ошибке
    """
    from .mendrivelibc.common import det_eval_fast

    re0, im0 = graphic_solutions[sol_num]
    det0_re, det0_im = det_eval_fast(lib, re0, im0)
    det0_residual = float(np.sqrt((det0_re)**2 + (det0_im)**2))

    if verbose:
        print('sol_index =', sol_num - len(graphic_solutions)/2)
        print("Root Init:", re0, im0, f"|det|={det0_residual:.4e}")

    # Уточняем корень методом Ньютона (первая стадия)
    if solver is not None:
        re_root_num, im_root_num = solver.solve(re0=re0, im0=im0)
        det_re, det_im = det_eval_fast(lib, re_root_num, im_root_num)
        sHl, sEl, sHr, sEr = get_sign_K_fast(lib)
        det_residual_num = float(np.sqrt((det_re)**2 + (det_im)**2))

        if verbose:
            print("Root Sol:", re_root_num, im_root_num, f"|det|={det_residual_num:.4e}")
    else:
        re_root_num = re0
        im_root_num = im0

    if precise:
        # Уточняем корень с повышенной точностью (вторая стадия)
        re_prec, im_prec = preciser.solve(re0=re_root_num, im0=im_root_num)
        det_re_prec, det_im_prec = det_eval_fast(lib, re_prec, im_prec)
        sHl, sEl, sHr, sEr = get_sign_K_fast(lib)
        det_residual_prec = float(np.sqrt((det_re_prec)**2 + (det_im_prec)**2))

        if verbose:
            print("Root Prec:", re_prec, im_prec, f"|det|={det_residual_prec:.4e}")
            print("sign K:", sHl, sEl, sHr, sEr)

        if det_residual_prec > 1e-5:
            print(f"    ⚠️ Большой детерминант системы: {det_residual_prec:.2e}")
            return None

        # Формируем результат
        result = {
            're': re_prec,
            'im': im_prec,
            'det_residual': det_residual_prec,
            'sign_K': (sHl, sEl, sHr, sEr),
            'sol_num': sol_num,
            're_init': re0,
            'im_init': im0,
            'det_residual_init': det0_residual
        }
    else:
        # Формируем результат
        result = {
            're': re_root_num,
            'im': im_root_num,
            'det_residual': det_residual_num,
            'sign_K': (sHl, sEl, sHr, sEr),
            'sol_num': sol_num,
            're_init': re0,
            'im_init': im0,
            'det_residual_init': det0_residual
        }


    if return_all_stages:
        result['stages'] = {
            'init': {'re': re0, 'im': im0, 'det': det0_residual},
            'newton': {'re': re_root_num, 'im': im_root_num, 'det': det_residual_num},
            'prec': {'re': re_prec, 'im': im_prec, 'det': det_residual_prec}
        }

    return result

def find_resonance_roots_refine_and_cluster_them(
        lib, param_overrides, solver, preciser,
        re_range, im_range,
        n_re, n_im,
        sharp,
        eps_sin_r_s=1e-3,
        window_size=5,
        cluster_eps=1e-3,
        use_clustering=True,
        clustering_method='simple',  # 'simple' или 'hierarchical'
        verbose=True,
        make_plot=False):
    """
    Выполняет расчёт корней с кластеризацией .

    Args:
        lib
        solver
        preciser
        re_range, im_range: диапазон поиска корней
        n_re, n_im: количество точек сетки
        sharp: параметр резкости для поиска корней
        window_size: размер окна для сглаживания
        cluster_eps: параметр eps для кластеризации (расстояние между корнями)
        use_clustering: использовать ли кластеризацию корней
        clustering_method: 'simple' (быстрый) или 'hierarchical' (более точный)
        verbose: выводить ли отладочную информацию
        make_plot: строить ли графики

    Returns:
        clustered_roots или None при ошибке
    """

    # 3. Находим резонансные корни
    if False:
        graphic_solutions = find_resonance_roots_fast(lib, re_range, im_range)
    else:
        graphic_solutions, (cu, cv, pl) = find_resonance_roots_fast_с(
            lib, re_range, im_range,
            n_re=n_re, n_im=n_im,
            make_plot=make_plot,
            use_tracing=0,
            sharp=sharp,
            eps_sin_r_s = eps_sin_r_s,
            window_size=window_size
        )

    if len(graphic_solutions) == 0:
        print(f"    ⚠️ Корни не найдены для параметров {param_overrides}")
        return None

    # 4. Уточняем все корни и собираем базу
    if verbose:
        print(f"\n  🔍 Уточнение {len(graphic_solutions)} графических решений для параметров {param_overrides}...")

    refined_roots = []
    for sol_num in range(len(graphic_solutions)):
        result = prec_graphic_solutions(
            lib, solver, preciser,
            sol_num, graphic_solutions, precise=False,
            verbose=verbose
        )

        if result is not None:
            refined_roots.append(result)

    if len(refined_roots) == 0:
        print(f"    ⚠️ Не удалось уточнить ни один корень")
        return None

    if verbose:
        print(f"  ✓ Успешно уточнено {len(refined_roots)} корней")

    # 5. Кластеризация корней (опционально)
    if use_clustering:
        if verbose:
            print(f"\n  🔬 Кластеризация корней (eps={cluster_eps}, метод={clustering_method})...")

        if clustering_method == 'hierarchical':
            clustered_roots = cluster_roots_hierarchical(refined_roots, eps=cluster_eps, verbose=verbose)
        else:
            clustered_roots = cluster_roots(refined_roots, eps=cluster_eps, verbose=verbose)
    else:
        clustered_roots = refined_roots

    if verbose:
        print(f"  ✓ Будет рассчитано {len(clustered_roots)} уникальных корней для параметров {param_overrides}\n")

    refined_clustered_roots = []
    for sol_num in range(len(clustered_roots)):
        result = prec_graphic_solutions(
            lib, solver, preciser,
            sol_num, graphic_solutions, precise=True,
            verbose=verbose
        )

        if result is not None:
            refined_clustered_roots.append(result)

    return graphic_solutions, refined_roots, clustered_roots, refined_clustered_roots, (cu, cv, pl)


# ============================================================================
# ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ ДЛЯ КЛАСТЕРИЗАЦИИ КОРНЕЙ
# ============================================================================

def euclidean_distance(root1, root2):
    """Вычисляет евклидово расстояние между двумя корнями."""
    return np.sqrt((root1['re'] - root2['re'])**2 + (root1['im'] - root2['im'])**2)


def cluster_roots(roots_list, eps=1e-3, verbose=True):
    """
    Кластеризует близкие корни используя простой алгоритм на основе расстояний.

    Args:
        roots_list: список словарей с ключами 'kz', 'sz', 'det_residual', 'kz_sol_num'
        eps: максимальное расстояние между корнями в одном кластере
        verbose: выводить ли информацию о кластеризации

    Returns:
        список представителей кластеров (лучший корень из каждого кластера)
    """
    if len(roots_list) == 0:
        return []

    # Копируем список, чтобы не изменять оригинал
    remaining_roots = roots_list.copy()
    clusters = []

    while remaining_roots:
        # Берём первый оставшийся корень как центр нового кластера
        seed = remaining_roots.pop(0)
        current_cluster = [seed]

        # Ищем все корни, близкие к seed
        i = 0
        while i < len(remaining_roots):
            if euclidean_distance(seed, remaining_roots[i]) <= eps:
                current_cluster.append(remaining_roots.pop(i))
            else:
                i += 1

        clusters.append(current_cluster)

    if verbose:
        n_clusters = len(clusters)
        total_roots = sum(len(c) for c in clusters)
        print(f"  📊 Найдено кластеров: {n_clusters} из {total_roots} корней")

    # Выбираем лучший корень из каждого кластера (с минимальным det_residual)
    representative_roots = []
    for cluster_idx, cluster_roots in enumerate(clusters):
        # Выбираем корень с минимальной невязкой детерминанта
        best_root = min(cluster_roots, key=lambda r: r['det_residual'])

        if verbose:
            avg_kz = np.mean([r['kz'] for r in cluster_roots])
            avg_sz = np.mean([r['sz'] for r in cluster_roots])
            print(f"    Кластер {cluster_idx}: {len(cluster_roots)} корней, "
                  f"центр: kz≈{avg_kz:.6f}, sz≈{avg_sz:.6f}")
            print(f"      → выбран: kz={best_root['kz']:.6f}, sz={best_root['sz']:.6f}, "
                  f"|det|={best_root['det_residual']:.2e}")

        representative_roots.append(best_root)

    return representative_roots

def cluster_roots_hierarchical(roots_list, eps=1e-3, verbose=True):
    """
    Альтернативная версия: иерархическая кластеризация.
    Объединяет ближайшие пары корней, пока расстояние не превысит eps.

    Args:
        roots_list: список словарей с ключами 'kz', 'sz', 'det_residual', 'kz_sol_num'
        eps: максимальное расстояние между корнями в одном кластере
        verbose: выводить ли информацию о кластеризации

    Returns:
        список представителей кластеров
    """
    if len(roots_list) == 0:
        return []

    # Каждый корень начинает в своём кластере
    clusters = [[root] for root in roots_list]

    while True:
        # Находим два ближайших кластера
        min_dist = float('inf')
        merge_i, merge_j = -1, -1

        for i in range(len(clusters)):
            for j in range(i + 1, len(clusters)):
                # Расстояние между кластерами = минимальное расстояние между их элементами
                dist = min(euclidean_distance(r1, r2)
                          for r1 in clusters[i]
                          for r2 in clusters[j])

                if dist < min_dist:
                    min_dist = dist
                    merge_i, merge_j = i, j

        # Если минимальное расстояние больше eps, останавливаемся
        if min_dist > eps:
            break

        # Объединяем два ближайших кластера
        if merge_i != -1:
            clusters[merge_i].extend(clusters[merge_j])
            clusters.pop(merge_j)
        else:
            break

    if verbose:
        n_clusters = len(clusters)
        total_roots = sum(len(c) for c in clusters)
        print(f"  📊 Найдено кластеров: {n_clusters} из {total_roots} корней (иерархический метод)")

    # Выбираем лучший корень из каждого кластера
    representative_roots = []
    for cluster_idx, cluster_roots in enumerate(clusters):
        best_root = min(cluster_roots, key=lambda r: r['det_residual'])

        if verbose:
            avg_kz = np.mean([r['kz'] for r in cluster_roots])
            avg_sz = np.mean([r['sz'] for r in cluster_roots])
            print(f"    Кластер {cluster_idx}: {len(cluster_roots)} корней, "
                  f"центр: kz≈{avg_kz:.6f}, sz≈{avg_sz:.6f}")
            print(f"      → выбран: kz={best_root['kz']:.6f}, sz={best_root['sz']:.6f}, "
                  f"|det|={best_root['det_residual']:.2e}")

        representative_roots.append(best_root)

    return representative_roots

class NewtonPrecC:
    def __init__(self, lib, precision):
        self.precision = precision
        from .mendrivelibc.common import get_numeric_type
        self.num_type = get_numeric_type(self.precision)

        from ctypes import POINTER, c_int

        self.lib = lib
        self.lib.newton_adaptive_step.argtypes = [
            POINTER(self.num_type),  # kz
            POINTER(self.num_type),  # sz
            POINTER(self.num_type),  # step_re_re
            POINTER(self.num_type),  # step_im_re
            POINTER(self.num_type),  # step_re_im
            POINTER(self.num_type),  # step_im_im
            POINTER(self.num_type),  # f_abs_out
            self.num_type,           # step_decrease
            self.num_type,           # step_increase
            self.num_type,           # delta_eps
            self.num_type,           # f_abs_eps
            c_int                   # max_retries
        ]

        self.lib.newton_adaptive_step.restype = c_int

    def solve(self, re0, im0, nsteps=1000, logging=False):
        from ctypes import byref

        kz = self.num_type(re0)
        sz = self.num_type(im0)

        # 4 адаптивных шага (как в Python newton_prec)
        step_init = 0.95
        step_re_re = self.num_type(step_init)
        step_im_re = self.num_type(step_init)
        step_re_im = self.num_type(step_init)
        step_im_im = self.num_type(step_init)

        f_abs = self.num_type()

        # Параметры алгоритма
        step_decrease = self.num_type(0.999)
        step_increase = self.num_type(0.999999)
        delta_eps     = self.num_type(1e-64)
        f_abs_eps     = self.num_type(1e-64)
        max_retries   = c_int(2)

        for i in range(nsteps):
            f_abs_prev = f_abs

            ret = self.lib.newton_adaptive_step(
                byref(kz), byref(sz),
                byref(step_re_re), byref(step_im_re),
                byref(step_re_im), byref(step_im_im),
                byref(f_abs),
                step_decrease, step_increase,
                delta_eps, f_abs_eps,
                max_retries
            )

            if logging:
                print(f"  iter {i}: kz={kz.value:.10e}, sz={sz.value:.10e}, "
                      f"|f|={f_abs.value:.36e}, ret={ret}")

            delta_f = abs(f_abs.value - f_abs_prev.value)
            # Защита от зацикливания
            if (delta_f < 1e-128 and i > 5 and -2 != ret):
                if logging:
                    print(f"  ⚠️  Сходимость остановилась (|Δf| {delta_f}" )
                break

            if ret == 1:
                if logging:
                    print(f"✅ converged by delta at step {i}, |f|={f_abs.value:.36e}")
                break
            if ret == 2:
                if logging:
                    print(f"✅ converged by |f| at step {i}, |f|={f_abs.value:.36e}")
                break
            if -1 == ret < 0:
                raise RuntimeError(f"Derivative near zero at step {i}")
        else:
            if logging:
                print(f"⚠️  max iterations ({nsteps}) reached, |f|={f_abs.value:.36e}")

        return real128(kz.value), real128(sz.value)

class NewtonRootAdaptiveSolver:
    def __init__(self, lib, precision):
        self.lib = lib
        self.precision = precision

        from .mendrivelibc.common import get_numeric_type
        self.num_type = get_numeric_type(precision)

        from ctypes import POINTER, c_int, byref

        # solve_newton_root
        self.lib.newton_adaptive.argtypes = [
            POINTER(self.num_type), POINTER(self.num_type),
        ]
        self.lib.newton_adaptive.restype = c_int

    def solve(self, re0, im0):
        re_sol = self.num_type(re0)
        im_sol = self.num_type(im0)
        ret = self.lib.newton_adaptive(byref(re_sol), byref(im_sol))

        return real128(re_sol.value), real128(im_sol.value)

class NewtonRootClassicSolver:
    def __init__(self, lib, precision):
        self.lib = lib
        self.precision = precision
        from .mendrivelibc.common import get_numeric_type
        self.num_type = get_numeric_type(precision)

        from ctypes import POINTER, c_int, byref

        # solve_newton_root
        self.lib.solve_newton_root_classic.argtypes = [
            self.num_type, self.num_type,
            POINTER(self.num_type), POINTER(self.num_type),
            POINTER(self.num_type), c_int
        ]
        self.lib.solve_newton_root_classic.restype = c_int

    def solve(self, re0, im0, max_iter=100):
        re_sol = self.num_type()
        im_sol = self.num_type()
        f_abs_out = self.num_type()
        from ctypes import byref
        ret = self.lib.solve_newton_root_classic(self.num_type(re0), self.num_type(im0),
                                    byref(re_sol), byref(im_sol),
                                    byref(f_abs_out), max_iter)

        return real128(re_sol.value), real128(im_sol.value)

class NewtonRootSolver:
    def __init__(self, lib):
        self.lib = lib

        from ctypes import POINTER, c_double, c_int, byref

        # solve_newton_root
        self.lib.solve_newton_root.argtypes = [
            POINTER(c_double), POINTER(c_double),
            POINTER(c_double), POINTER(c_double),
            c_double, c_int
        ]
        self.lib.solve_newton_root.restype = c_int

    def solve(self, re0, im0, epsabs=1e-12, max_iter=100):
        from ctypes import c_double, c_int, byref
        re = c_double(re0)
        im = c_double(im0)
        omegare_sol = c_double()
        gamma_sol = c_double()
        ret = self.lib.solve_newton_root(byref(re), byref(im),
                                    byref(omegare_sol), byref(gamma_sol),
                                    c_double(epsabs), c_int(max_iter))

        return real128(omegare_sol.value), real128(gamma_sol.value)



def fix_K_sign(K):
    K_d = []
    for K_i in K:
        K_val = K_i.rhs()
        if K_val.imag().n(prec=128) < 0:
            K_d_i = K_i.lhs() == -K_val
        else:
            K_d_i = K_i.lhs() == K_val
        K_d += [K_d_i]
    return K_d

def svd_X4_old(M4_d, verbose = False, phase_ref_index=0, cond_threshold=1e12, residual_threshold=1e-8):
    """
    Надёжное вычисление приближённого ядра 4×4 комплексной матрицы через SVD.

    Параметры:
        M4_d: символьная или численная матрица 4×4
        verbose: вывод отладочной информации
        phase_ref_index: индекс компоненты, по которой нормализуется глобальная фаза
        cond_threshold: порог числа обусловленности (если выше — решение неустойчиво)
        residual_threshold: порог невязки ||M·X|| (если выше — режим не существует)

    Возвращает:
        X4_norm: нормализованный вектор ядра (длина 4, комплексный)
        result: dict с ключами:
            - 'residual': ||M·X||
            - 'cond': число обусловленности
            - 'svals': сингулярные значения
            - 'valid': bool — можно ли доверять решению
    """
    import numpy as np
    from scipy.linalg import svd

    # Преобразуем символьную матрицу M4 в numpy-массив комплексных чисел
    M4_np = np.array(M4_d, dtype=np.complex128)  # M4_d — ваша численная матрица 4×4

    # SVD: M = U · Σ · Vᴴ
    U, S, Vh = svd(M4_np)

    # Последняя строка Vh — собственный вектор для минимального сингулярного значения
    X4 = Vh[-1]  # комплексный вектор длины 4: [B1_zl, A1_z, A2_z, B1_zr]

    # Нормируем (например, по модулю A2_z)
    # X4 = X4 / X4[2]  # A2_z → 1

    residual = np.linalg.norm(M4_np @ X4)

    # Диагностика

    print("Сингулярные значения:", S)
    cond = S[0] / S[-1]  if S[-1] != 0 else np.inf # число обусловленности
    print("Число обусловленности:", cond)
    valid = (cond < cond_threshold) and (residual < residual_threshold)

    phase_ref = np.angle(X4[phase_ref_index])
    print("phase_ref:", phase_ref)

    # Нормировка по модулю
    norm = np.linalg.norm(X4)
    print("norm:", norm)

    if verbose:
        print("Коэффициенты:", X4)
        print("Невязка:", residual)  # должна быть ~1e-15

    return X4, residual

def svd_X4(M4_d, verbose=False, residual_threshold=1e-8):
    """
    Ищет решение M·X = 0 и автоматически исправляет фазовые/перестановочные ошибки.

    Поддерживает исправление:
      - глобального фазового сдвига (умножение на ±i)
      - перестановки A1 ↔ A2 (для симметричных мод)
    """
    import numpy as np
    from scipy.linalg import svd

    M4_np = np.array(M4_d, dtype=np.complex128)
    U, s, Vh = svd(M4_np)
    X4_raw = Vh[-1]
    residual_raw = np.linalg.norm(M4_np @ X4_raw)

    # Кандидаты на исправление
    candidates = [
        X4_raw,                     # оригинал
        1j * X4_raw,                # фаза +π/2
        -1j * X4_raw,               # фаза -π/2
        np.conj(X4_raw),            # комплексное сопряжение
    ]

    # Для 4-компонентного вектора [B_l, A1, A2, B_r] — пробуем A1↔A2
    if len(X4_raw) == 4:
        X_swapped = np.array([X4_raw[0], X4_raw[2], X4_raw[1], X4_raw[3]])
        candidates.extend([
            X_swapped,
            1j * X_swapped,
            -1j * X_swapped,
            np.conj(X_swapped),
        ])

    best_X = None
    best_residual = np.inf
    best_candidate = None

    for i, X in enumerate(candidates):
        # Нормализуем для устойчивости
        norm = np.linalg.norm(X)
        if norm < 1e-15:
            continue
        X_norm = X / norm
        res = np.linalg.norm(M4_np @ X_norm)
        if res < best_residual:
            best_residual = res
            best_X = X_norm
            best_candidate = i

    valid = best_residual < residual_threshold

    if verbose:
        print(f"Лучший кандидат: #{best_candidate}, остаток = {best_residual:.2e}")
        if best_candidate != 0:
            print("⚠️ Применена коррекция фазы/перестановки")

    return best_X, {'residual': best_residual, 'valid': valid, 'candidate_id': best_candidate}