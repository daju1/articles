def save_detplots(detplots, param_values, param_name='omega',
                  output_dir='detplots', file_format='png',
                  figsize=(8, 6), dpi=150, show=False):
    """
    Visualize and save determinant contour plots from a scan.

    Parameters
    ----------
    detplots : list of tuples (cu, cv, pl)
        Output from scan_parameter() or similar — each entry corresponds
        to one parameter value. cu/cv are Re=0/Im=0 contour point lists,
        pl is a SageMath Graphics object.
    param_values : list or array
        Parameter values corresponding to each entry in detplots.
        Must have the same length as detplots.
    param_name : str
        Name of the scanned parameter (used in filenames and titles).
    output_dir : str
        Directory to save plots into. Created if it doesn't exist.
    file_format : str
        Image format: 'png', 'pdf', 'svg', etc.
    figsize : tuple
        Figure size in inches (width, height).
    dpi : int
        Resolution for raster formats (png, jpg).
    show : bool
        If True, display each plot inline as well as saving it.
        Useful for spot-checking; set False for bulk saves.

    Returns
    -------
    saved_paths : list of str
        Paths to all saved files.
    """
    import os

    os.makedirs(output_dir, exist_ok=True)
    saved_paths = []

    if len(detplots) != len(param_values):
        raise ValueError(
            f"detplots has {len(detplots)} entries but param_values has "
            f"{len(param_values)} — they must match."
        )

    for idx, ((cu, cv, pl), pval) in enumerate(zip(detplots, param_values)):
        # Build a human-readable label for the parameter value
        pval_str = f"{float(pval):.4g}".replace('+', '').replace(' ', '')

        # Count contour intersections (approximate root count)
        n_re_contours = len(cu) if cu else 0
        n_im_contours = len(cv) if cv else 0

        # Compose title
        title = (
            f"{param_name} = {pval_str}  |  "
            f"Re=0 segments: {n_re_contours}, Im=0 segments: {n_im_contours}"
        )

        # SageMath Graphics objects accept a title kwarg in show()/save()
        filename = f"{param_name}_{idx:04d}_{pval_str}.{file_format}"
        filepath = os.path.join(output_dir, filename)

        try:
            pl.save(
                filepath,
                title=title,
                figsize=figsize,
                dpi=dpi
            )
            saved_paths.append(filepath)
            print(f"[{idx+1}/{len(detplots)}] Saved: {filepath}")
        except Exception as e:
            print(f"[{idx+1}/{len(detplots)}] ERROR saving {filepath}: {e}")

        if show:
            pl.show(title=title, figsize=figsize)

    print(f"\nDone. {len(saved_paths)}/{len(detplots)} plots saved to '{output_dir}/'")
    return saved_paths

def match_branches_between_intervals(results, param_name, distance_threshold=0.1,
                                     use_K_params=False):
    """
    Сопоставляет ветви между соседними интервалами по близости корней и K-параметров.

    Parameters:
    -----------
    results : list
        Список результатов из scan_parameter
    param_name : str
        Имя параметра сканирования
    distance_threshold : float
        Максимальное расстояние для считывания ветвей одинаковыми
    use_K_params : bool
        Использовать ли K_left, K_vacuum, K_right для сопоставления

    Returns:
    --------
    matched_branches : dict
        Словарь {global_branch_id: list of data points}
    """

    def extract_K_params(result):
        """Извлекает K-параметры из результата."""
        K_params = {}
        for key in ['K_left', 'K_vacuum', 'K_right']:
            if key in result:
                for eq in result[key]:
                    val = eq.rhs()
                    # Преобразуем в комплексное число
                    if hasattr(val, 'real') and hasattr(val, 'imag'):
                        K_params[key] = complex(float(val.real()), float(val.imag()))
                    else:
                        try:
                            K_params[key] = complex(val)
                        except:
                            K_params[key] = complex(float(val), 0)
            else:
                K_params[key] = complex(0, 0)
        return K_params

    def compute_distance(point1, point2):
        """Вычисляет расстояние между двумя точками."""
        import numpy as np
        # Расстояние по kz, sz
        kz_dist = (point1['kz'] - point2['kz'])**2
        sz_dist = (point1['sz'] - point2['sz'])**2
        root_dist = np.sqrt(kz_dist + sz_dist)
        return root_dist


    # Сначала собираем все данные по интервалам
    intervals_data = []

    for results_list in results:
        if not results_list:
            continue

        interval_branches = {}

        for result in results_list:
            if 'cluster_index' not in result:
                continue

            cluster_idx = result['cluster_index']

            if cluster_idx not in interval_branches:
                interval_branches[cluster_idx] = []

            K_params = extract_K_params(result)

            # Основная точка
            if 'params' in result and param_name in result['params']:
                point = {
                    'x': float(result['params'][param_name]),
                    'thrust': float(result['thrust_N_per_kW']),
                    'thrust_left': float(result['thrust_N_per_kW_left']),
                    'thrust_right': float(result['thrust_N_per_kW_right']),
                    'kz': float(result['kz']),
                    'sz': float(result['sz']),
                    'det_residual': float(result['det_residual']),
                    'K_params': K_params
                }
                interval_branches[cluster_idx].append(point)

            # Аппендикс
            if 'thrust_appendix' in result and len(result['thrust_appendix']) > 0:
                for app_point in result['thrust_appendix']:
                    if 'params' in app_point and param_name in app_point['params']:
                        K_params_app = extract_K_params(app_point)
                        point = {
                            'x': float(app_point['params'][param_name]),
                            'thrust': float(app_point['thrust_N_per_kW']),
                            'thrust_left': float(app_point['thrust_N_per_kW_left']),
                            'thrust_right': float(app_point['thrust_N_per_kW_right']),
                            'kz': float(app_point['kz']),
                            'sz': float(app_point['sz']),
                            'det_residual': float(app_point['det_residual']),
                            'K_params': K_params_app
                        }
                        interval_branches[cluster_idx].append(point)

        if interval_branches:
            intervals_data.append(interval_branches)

    if len(intervals_data) == 0:
        return {}

    # Теперь сопоставляем ветви между интервалами
    matched_branches = {}
    next_global_id = 0

    # Первый интервал - инициализация
    for local_idx, points in intervals_data[0].items():
        matched_branches[next_global_id] = points
        next_global_id += 1

    # Сопоставляем последующие интервалы
    for interval_idx in range(1, len(intervals_data)):
        current_interval = intervals_data[interval_idx]

        # Для каждой ветви в текущем интервале ищем ближайшую в предыдущих
        unmatched_locals = set(current_interval.keys())

        for global_id, global_points in matched_branches.items():
            if len(global_points) == 0:
                continue

            # Берём последнюю точку из глобальной ветви
            last_point = global_points[-1]

            # Ищем ближайшую локальную ветвь
            best_local_idx = None
            best_distance = float('inf')

            for local_idx in unmatched_locals:
                local_points = current_interval[local_idx]
                if len(local_points) == 0:
                    continue

                # Берём первую точку локальной ветви
                first_point = local_points[0]

                # Вычисляем расстояние
                distance = compute_distance(last_point, first_point)

                if distance < best_distance:
                    best_distance = distance
                    best_local_idx = local_idx

            # Если нашли достаточно близкую ветвь, добавляем к глобальной
            if best_local_idx is not None and best_distance < distance_threshold:
                matched_branches[global_id].extend(current_interval[best_local_idx])
                unmatched_locals.remove(best_local_idx)

        # Оставшиеся локальные ветви создают новые глобальные
        for local_idx in unmatched_locals:
            matched_branches[next_global_id] = current_interval[local_idx]
            next_global_id += 1

    return matched_branches


def get_branches_count(results, param_name, distance_threshold=0.1,
                       use_K_params=False):
    """
    Возвращает количество сопоставленных ветвей.

    Parameters:
    -----------
    results : list
        Список результатов из scan_parameter
    param_name : str
        Имя параметра сканирования
    distance_threshold : float
        Максимальное расстояние для считывания ветвей одинаковыми
    use_K_params : bool
        Использовать ли K-параметры для сопоставления

    Returns:
    --------
    count : int
        Количество ветвей
    """
    matched_branches = match_branches_between_intervals(
        results, param_name, distance_threshold, use_K_params,
    )
    return len(matched_branches)


def plot_matched_branch_results(results, param_name, figfilename, branch_id=0,
                                distance_threshold=0.1,
                                use_K_params=False,
                                language='ru'):
    """
    Построение графиков для конкретной сопоставленной ветви решения.
    Включает графики K-параметров и улучшенный график детерминанта с промежуточными точками.

    Parameters:
    -----------
    results : list
        Список результатов из scan_parameter
    param_name : str
        Имя параметра сканирования
    branch_id : int
        ID глобальной ветви для отображения
    distance_threshold : float
        Порог расстояния для сопоставления ветвей
    use_K_params : bool
        Использовать ли K-параметры для сопоставления
    language : str
        'ru' или 'en'
    """

    translations = {
        'ru': {
            'no_data': "Нет данных для построения графика",
            'branch': 'Ветвь',
            'specific_thrust': 'Удельная тяга (Н/кВт)',
            'thrust_vs': 'Зависимость тяги от',
            'kz_vs': 'Зависимость kz от',
            'sz_vs': 'Зависимость sz от',
            'det_vs': 'Невязка детерминанта от',
            'K_vs': 'Зависимость K-параметров от',
            'left': 'левая',
            'right': 'правая',
            'total': 'общая',
            'real_part': 'Re(kz)',
            'imag_part': 'Im(kz) = sz',
            'determinant': '|det|',
            'K_real': 'Re(K)',
            'K_imag': 'Im(K)',
            'K_abs': '|K|',
            'main_points': 'основные точки',
            'appendix_points': 'промежуточные точки'
        },
        'en': {
            'no_data': "No data to plot",
            'branch': 'Branch',
            'specific_thrust': 'Specific Thrust (N/kW)',
            'thrust_vs': 'Thrust vs',
            'kz_vs': 'kz vs',
            'sz_vs': 'sz vs',
            'det_vs': 'Determinant residual vs',
            'K_vs': 'K-parameters vs',
            'left': 'left',
            'right': 'right',
            'total': 'total',
            'real_part': 'Re(kz)',
            'imag_part': 'Im(kz) = sz',
            'determinant': '|det|',
            'K_real': 'Re(K)',
            'K_imag': 'Im(K)',
            'K_abs': '|K|',
            'main_points': 'main points',
            'appendix_points': 'intermediate points'
        }
    }

    lang = translations.get(language, translations['ru'])

    # Сопоставляем ветви
    matched_branches = match_branches_between_intervals(
        results, param_name, distance_threshold, use_K_params
    )

    if branch_id not in matched_branches:
        print(f"{lang['no_data']} для ветви {branch_id}")
        print(f"Доступные ветви: {list(matched_branches.keys())}")
        return

    branch_data = matched_branches[branch_id]

    if len(branch_data) == 0:
        print(f"{lang['no_data']} для ветви {branch_id}")
        return

    # Сортируем данные по x
    branch_data = sorted(branch_data, key=lambda p: p['x'])

    # Извлекаем данные
    x            = [p['x'] for p in branch_data]
    thrust       = [p['thrust'] for p in branch_data]
    thrust_left  = [p['thrust_left'] for p in branch_data]
    thrust_right = [p['thrust_right'] for p in branch_data]
    kz           = [p['kz'] for p in branch_data]
    sz           = [p['sz'] for p in branch_data]
    det_residual = [p['det_residual'] for p in branch_data]

    if use_K_params:
        # Извлекаем K-параметры
        K_left   = [p['K_params']['K_left'] for p in branch_data]
        K_vacuum = [p['K_params']['K_vacuum'] for p in branch_data]
        K_right  = [p['K_params']['K_right'] for p in branch_data]

        grid_rows = 4
        # Создаём графики - увеличенная сетка для K-параметров
        fig = plt.figure(figsize=(20, 14))
    else:
        grid_rows = 3
        # Создаём графики
        fig = plt.figure(figsize=(18, 10))

    gs = fig.add_gridspec(grid_rows, 3, hspace=0.35, wspace=0.3)

    grid_row = 0
    # График 1: Тяга (все компоненты)
    ax1 = fig.add_subplot(gs[grid_row, :2])
    ax1.semilogx(x, thrust, 'o-', linewidth=2, markersize=6, color='blue', label=lang['total'])
    ax1.semilogx(x, thrust_left, 's-', linewidth=1.5, markersize=5, color='red', alpha=0.7, label=lang['left'])
    ax1.semilogx(x, thrust_right, '^-', linewidth=1.5, markersize=5, color='green', alpha=0.7, label=lang['right'])
    ax1.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    ax1.set_xlabel(param_name, fontsize=12)
    ax1.set_ylabel(lang['specific_thrust'], fontsize=12)
    ax1.set_title(f"{lang['branch']} {branch_id}: {lang['thrust_vs']} {param_name}", fontsize=14)
    ax1.grid(True, which='both', linestyle=':', alpha=0.7)
    ax1.legend()

    # График 2: kz (действительная часть)
    ax2 = fig.add_subplot(gs[grid_row, 2])
    ax2.semilogx(x, kz, 'o-', linewidth=2, markersize=6, color='purple')
    ax2.set_xlabel(param_name, fontsize=12)
    ax2.set_ylabel(lang['real_part'], fontsize=12)
    ax2.set_title(f"{lang['kz_vs']} {param_name}", fontsize=14)
    ax2.grid(True, which='both', linestyle=':', alpha=0.7)
    grid_row += 1

    # График 3: sz (мнимая часть)
    ax3 = fig.add_subplot(gs[grid_row, 0])
    ax3.semilogx(x, sz, 'o-', linewidth=2, markersize=6, color='orange')
    ax3.set_xlabel(param_name, fontsize=12)
    ax3.set_ylabel(lang['imag_part'], fontsize=12)
    ax3.set_title(f"{lang['sz_vs']} {param_name}", fontsize=14)
    ax3.grid(True, which='both', linestyle=':', alpha=0.7)

    # График 4: Невязка детерминанта (с ВСЕМИ точками включая промежуточные)
    ax4 = fig.add_subplot(gs[grid_row, 1])
    # Разделяем на основные и промежуточные точки для разного отображения
    valid_data = [(x[i], det_residual[i], i) for i in range(len(det_residual))
                  if not np.isnan(det_residual[i]) and det_residual[i] > 0]

    if len(valid_data) > 0:
        # Определяем основные точки (первая точка каждого интервала сканирования)
        # Предполагаем, что основные точки имеют большие скачки в x
        x_diffs = [valid_data[i+1][0] / valid_data[i][0] if i < len(valid_data)-1 else 1
                   for i in range(len(valid_data))]
        threshold_ratio = 1.5  # Если следующая точка в 1.5+ раз больше, это скачок к новому интервалу

        main_x, main_det = [], []
        inter_x, inter_det = [], []

        for i, (xi, deti, idx) in enumerate(valid_data):
            if i == 0 or x_diffs[i-1] > threshold_ratio:
                # Основная точка
                main_x.append(xi)
                main_det.append(deti)
            else:
                # Промежуточная точка
                inter_x.append(xi)
                inter_det.append(deti)

        # Рисуем все точки линией
        all_x = [d[0] for d in valid_data]
        all_det = [d[1] for d in valid_data]
        ax4.loglog(all_x, all_det, '-', linewidth=1.5, color='brown', alpha=0.5)

        # Основные точки крупными маркерами
        if len(main_x) > 0:
            ax4.loglog(main_x, main_det, 'o', markersize=8, color='brown',
                      label=lang['main_points'], zorder=3)

        # Промежуточные точки мелкими маркерами
        if len(inter_x) > 0:
            ax4.loglog(inter_x, inter_det, 's', markersize=4, color='brown',
                      alpha=0.6, label=lang['appendix_points'], zorder=2)

        ax4.set_xlabel(param_name, fontsize=12)
        ax4.set_ylabel(lang['determinant'], fontsize=12)
        ax4.set_title(f"{lang['det_vs']} {param_name}", fontsize=14)
        ax4.grid(True, which='both', linestyle=':', alpha=0.7)
        ax4.legend()

    # График 5: Траектория в комплексной плоскости (kz, sz)
    ax5 = fig.add_subplot(gs[grid_row, 2])
    scatter = ax5.scatter(kz, sz, c=x, cmap='viridis', s=50)
    ax5.plot(kz, sz, '-', linewidth=1, color='gray', alpha=0.5)
    ax5.set_xlabel(lang['real_part'], fontsize=12)
    ax5.set_ylabel(lang['imag_part'], fontsize=12)
    ax5.set_title(f"Траектория в плоскости (kz, sz)", fontsize=14)
    ax5.grid(True, linestyle=':', alpha=0.7)
    cbar = plt.colorbar(scatter, ax=ax5)
    cbar.set_label(param_name, fontsize=10)

    if use_K_params:
        # Графики 6-8: K-параметры
        # График 6: K_left
        grid_row += 1
        ax6 = fig.add_subplot(gs[grid_row, 0])
        K_left_re  = [k.real for k in K_left]
        K_left_im  = [k.imag for k in K_left]
        K_left_abs = [abs(k) for k in K_left]
        #ax6.semilogx(x, K_left_re, 'o-', linewidth=2, markersize=5, color='blue', label=lang['K_real'])
        ax6.semilogx(x, K_left_im, 's-', linewidth=2, markersize=5, color='red', label=lang['K_imag'])
        #ax6.semilogx(x, K_left_abs, '^-', linewidth=2, markersize=5, color='green', label=lang['K_abs'])
        ax6.set_xlabel(param_name, fontsize=12)
        ax6.set_ylabel('K_left', fontsize=12)
        ax6.set_title(f"K_left vs {param_name}", fontsize=14)
        ax6.grid(True, which='both', linestyle=':', alpha=0.7)
        ax6.legend()

        # График 7: K_vacuum
        ax7 = fig.add_subplot(gs[grid_row, 1])
        K_vacuum_re  = [k.real for k in K_vacuum]
        K_vacuum_im  = [k.imag for k in K_vacuum]
        K_vacuum_abs = [abs(k) for k in K_vacuum]
        ax7.semilogx(x, K_vacuum_re, 'o-', linewidth=2, markersize=5, color='blue', label=lang['K_real'])
        ax7.semilogx(x, K_vacuum_im, 's-', linewidth=2, markersize=5, color='red', label=lang['K_imag'])
        #ax7.semilogx(x, K_vacuum_abs, '^-', linewidth=2, markersize=5, color='green', label=lang['K_abs'])
        ax7.set_xlabel(param_name, fontsize=12)
        ax7.set_ylabel('K_vacuum', fontsize=12)
        ax7.set_title(f"K_vacuum vs {param_name}", fontsize=14)
        ax7.grid(True, which='both', linestyle=':', alpha=0.7)
        ax7.legend()

        # График 8: K_right
        ax8 = fig.add_subplot(gs[grid_row, 2])
        K_right_re  = [k.real for k in K_right]
        K_right_im  = [k.imag for k in K_right]
        K_right_abs = [abs(k) for k in K_right]
        #ax8.semilogx(x, K_right_re, 'o-', linewidth=2, markersize=5, color='blue', label=lang['K_real'])
        ax8.semilogx(x, K_right_im, 's-', linewidth=2, markersize=5, color='red', label=lang['K_imag'])
        #ax8.semilogx(x, K_right_abs, '^-', linewidth=2, markersize=5, color='green', label=lang['K_abs'])
        ax8.set_xlabel(param_name, fontsize=12)
        ax8.set_ylabel('K_right', fontsize=12)
        ax8.set_title(f"K_right vs {param_name}", fontsize=14)
        ax8.grid(True, which='both', linestyle=':', alpha=0.7)
        ax8.legend()

    grid_row += 1
    # График 9: Отдельно общая тяга (крупно)
    ax9 = fig.add_subplot(gs[grid_row, :])
    ax9.semilogx(x, thrust, 'o-', linewidth=2.5, markersize=8, color='blue')
    ax9.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    ax9.fill_between(x, 0, thrust, alpha=0.3, color='blue')
    ax9.set_xlabel(param_name, fontsize=12)
    ax9.set_ylabel(lang['specific_thrust'], fontsize=12)
    ax9.set_title(f"{lang['branch']} {branch_id}: {lang['total']} {lang['thrust_vs']} {param_name}",
                  fontsize=14)
    ax9.grid(True, which='both', linestyle=':', alpha=0.7)

    plt.suptitle(f"{lang['branch']} {branch_id} - детальный анализ (сопоставленная)",
                 fontsize=16, fontweight='bold')
    # plt.show()
    plt.savefig(f"{figdir}/{figfilename}")
    plt.close()


def plot_all_branches_with_detplots(
    scan_results,
    detplots,
    param_values,
    param_name='omega',
    distance_threshold=0.1,
    output_dir="branch_detplots",
    figsize=(800, 600),
    marker_size=50,
    marker_color='red',
    thrust_key='thrust_N_per_kW',
):
    import os
    from sage.plot.point import point2d
    from sage.plot.line import line2d
    from sage.plot.text import text

    os.makedirs(output_dir, exist_ok=True)

    # --- строим глобальные ветви той же логикой что и plot_matched_branch_results ---
    matched_branches = match_branches_between_intervals(
        scan_results, param_name, distance_threshold
    )

    # --- индекс детплотов по значению omega ---
    detplot_by_omega = {float(pval): detplots[i] for i, pval in enumerate(param_values)}

    # --- для каждой глобальной ветви итерируем по точкам ---
    for global_id, branch_points in matched_branches.items():
        if not branch_points:
            continue

        for point in branch_points:
            omega      = point['x']
            kz_f       = point['kz']
            sz_f       = point['sz']
            thrust_val = point['thrust']

            detplot_tuple = detplot_by_omega.get(omega)
            if detplot_tuple is None:
                # ищем ближайший omega если float-ключи не совпадают точно
                closest_omega = min(detplot_by_omega.keys(), key=lambda o: abs(o - omega))
                if abs(closest_omega - omega) < 1e-6 * abs(omega + 1):
                    detplot_tuple = detplot_by_omega[closest_omega]
                else:
                    print(f"  [branch={global_id}, omega={omega:.4g}] No matching detplot, skipping.")
                    continue

            cu, cv, pl_base = detplot_tuple

            # --- границы из контурных сегментов ---
            all_pts = [pt for seg in cu for pt in seg] + [pt for seg in cv for pt in seg]
            if all_pts:
                xs = [float(pt[0]) for pt in all_pts]
                ys = [float(pt[1]) for pt in all_pts]
                kz_min, kz_max = min(xs), max(xs)
                sz_min, sz_max = min(ys), max(ys)
            else:
                margin = 0.5
                kz_min, kz_max = kz_f - margin, kz_f + margin
                sz_min, sz_max = sz_f - margin, sz_f + margin

            # --- маркер ---
            root_marker = point2d(
                (kz_f, sz_f),
                size=marker_size,
                color=marker_color,
                zorder=10,
            )

            # --- перекрестие ---
            crosshair = (
                line2d(
                    [(kz_min, sz_f), (kz_max, sz_f)],
                    color=marker_color, linestyle='--', thickness=0.8,
                )
                + line2d(
                    [(kz_f, sz_min), (kz_f, sz_max)],
                    color=marker_color, linestyle='--', thickness=0.8,
                )
            )

            # --- аннотация ---
            label_offset_kz = (kz_max - kz_min) * 0.02
            label_offset_sz = (sz_max - sz_min) * 0.02
            annotation = text(
                f"branch{global_id} ({kz_f:.3f}, {sz_f:.3f})\n{thrust_val:.2f} N/kW",
                (kz_f + label_offset_kz, sz_f + label_offset_sz),
                fontsize=8,
                color=marker_color,
                horizontal_alignment='left',
            )

            title_str = (
                f"{param_name}={omega:.4g}  branch={global_id}  "
                f"kz={kz_f:.5f}  sz={sz_f:.5f}  "
                f"thrust={thrust_val:.4g} N/kW"
            )

            combined = pl_base + root_marker + crosshair + annotation

            # --- имя файла: branch_id первым — удобно сортировать ---
            fname = (
                f"branch{global_id}_{param_name}_{omega:.6g}"
                f"_kz{kz_f:.5f}_sz{sz_f:.5f}.png"
            )
            fpath = os.path.join(output_dir, fname)

            combined.save(
                fpath,
                title=title_str,
                figsize=(figsize[0] / 100, figsize[1] / 100),
                dpi=100,
            )
            print(f"  Saved: {fname}  |  thrust={thrust_val:.4g} N/kW")

        print(f"[branch={global_id}] Done — {len(branch_points)} point(s) saved.")


# функция для покомпонентного анализа тяги каждой ветви:

def plot_branch_thrust_components(results, param_name, figfilename,
                                  branch_id=0,
                                  distance_threshold=0.1,
                                  use_K_params=False,
                                  language='ru'):
    """
    Построение покомпонентных графиков тяги для конкретной ветви.
    Показывает все компоненты силы, которые складываются в общую тягу.

    Parameters:
    -----------
    results : list
        Список результатов из scan_parameter
    param_name : str
        Имя параметра сканирования
    branch_id : int
        ID глобальной ветви для отображения
    distance_threshold : float
        Порог расстояния для сопоставления ветвей
    use_K_params : bool
        Использовать ли K-параметры для сопоставления
    language : str
        'ru' или 'en'
    """

    translations = {
        'ru': {
            'no_data': "Нет данных для построения графика",
            'branch': 'Ветвь',
            'components': 'Компоненты тяги',
            'left_side': 'Левая сторона',
            'right_side': 'Правая сторона',
            'total': 'Общая',
            'magnetic_eddy': 'Магнитный ток Фуко (F_jmD)',
            'electric_eddy_B': 'Электрический ток Фуко в B (p_Bt_jeB)',
            'electric_eddy_H': 'Электрический ток Фуко в H (p_Ht_jeH)',
            'convective_I_tang': 'Конвективная сила I танг. (p_I_conv_Ht)',
            'convective_I_norm': 'Конвективная сила I норм. (t_I_conv_Hn)',
            'rotforce_I': 'Роторная сила намагниченности (F_I_rotB)',
            'rotforce_P': 'Роторная сила поляризации (F_P_rotE)',
            'dynamic_I': 'Динамика намагниченности (F_dyn_I_D)',
            'dynamic_P': 'Динамика поляризации (F_dyn_P_B)',
            'tension_D': 'Натяжение D (t_Dn)',
            'pressure_E': 'Давление E (p_Et)',
            'vs': 'от',
            'specific_thrust': 'Удельная тяга (Н/кВт)'
        },
        'en': {
            'no_data': "No data to plot",
            'branch': 'Branch',
            'components': 'Thrust Components',
            'left_side': 'Left Side',
            'right_side': 'Right Side',
            'total': 'Total',
            'magnetic_eddy': 'Magnetic eddy current (F_jmD)',
            'electric_eddy_B': 'Electric eddy current in B (p_Bt_jeB)',
            'electric_eddy_H': 'Electric eddy current in H (p_Ht_jeH)',
            'convective_I_tang': 'Convective force I tang. (p_I_conv_Ht)',
            'convective_I_norm': 'Convective force I norm. (t_I_conv_Hn)',
            'rotforce_I': 'Rotor force magnetization (F_I_rotB)',
            'rotforce_P': 'Rotor force polarization (F_P_rotE)',
            'dynamic_I': 'Dynamic magnetization (F_dyn_I_D)',
            'dynamic_P': 'Dynamic polarization (F_dyn_P_B)',
            'tension_D': 'Tension D (t_Dn)',
            'pressure_E': 'Pressure E (p_Et)',
            'vs': 'vs',
            'specific_thrust': 'Specific Thrust (N/kW)'
        }
    }

    lang = translations.get(language, translations['ru'])

    # Сопоставляем ветви
    matched_branches = match_branches_between_intervals(
        results, param_name, distance_threshold, use_K_params
    )

    if branch_id not in matched_branches:
        print(f"{lang['no_data']} для ветви {branch_id}")
        print(f"Доступные ветви: {list(matched_branches.keys())}")
        return

    # Получаем точки ветви для сопоставления
    branch_points = matched_branches[branch_id]

    # Создаем словарь для быстрого поиска точек по (x, kz)
    branch_points_map = {(round(p['x'], 10), round(p['kz'], 6)): p for p in branch_points}

    # Собираем все компоненты
    components_data = {
        'x': [],
        'total': [],
        'total_left': [],
        'total_right': [],
        # Компоненты левой стороны
        'F_jmD_left': [],
        'p_Bt_jeB_left': [],
        'p_Ht_jeH_left': [],
        'p_I_conv_Ht_left': [],
        't_I_conv_Hn_left': [],
        'F_I_rotB_left': [],
        'F_dyn_I_D_left': [],
        'F_P_rotE_left': [],
        'F_dyn_P_B_left': [],
        't_Dn_left': [],
        'p_Et_left': [],
        # Компоненты правой стороны
        'F_jmD_right': [],
        'p_Bt_jeB_right': [],
        'p_Ht_jeH_right': [],
        'p_I_conv_Ht_right': [],
        't_I_conv_Hn_right': [],
        'F_I_rotB_right': [],
        'F_dyn_I_D_right': [],
        'F_P_rotE_right': [],
        'F_dyn_P_B_right': [],
        't_Dn_right': [],
        'p_Et_right': []
    }

    # Проходим по всем результатам и собираем данные
    for results_list in results:
        if not results_list:
            continue

        for result in results_list:
            # Обрабатываем основную точку
            if 'params' in result and param_name in result['params'] and 'kz' in result:
                x_val = float(result['params'][param_name])
                kz_val = float(result['kz'])
                key = (round(x_val, 10), round(kz_val, 6))

                if key in branch_points_map:
                    _add_component_data(components_data, result, x_val)

            # Обрабатываем аппендикс
            if 'thrust_appendix' in result and len(result['thrust_appendix']) > 0:
                for app_point in result['thrust_appendix']:
                    if 'params' in app_point and param_name in app_point['params'] and 'kz' in app_point:
                        x_app = float(app_point['params'][param_name])
                        kz_app = float(app_point['kz'])
                        key_app = (round(x_app, 10), round(kz_app, 6))

                        if key_app in branch_points_map:
                            _add_component_data(components_data, app_point, x_app)

    if len(components_data['x']) == 0:
        print(f"{lang['no_data']} для ветви {branch_id}")
        print(f"Проверьте, что результаты содержат все необходимые компоненты тяги")
        return

    # Сортируем по x
    sorted_indices = np.argsort(components_data['x'])
    for key in components_data:
        components_data[key] = [components_data[key][i] for i in sorted_indices]

    x = components_data['x']

    print(f"Найдено {len(x)} точек для ветви {branch_id}")

    # Создаём графики
    fig = plt.figure(figsize=(20, 16))


    grid_rows = 5

    gs = fig.add_gridspec(grid_rows, 2, hspace=0.4, wspace=0.3)

    # График 1: Общая тяга (левая + правая)
    grid_row = 0
    ax1 = fig.add_subplot(gs[grid_row, :])
    ax1.semilogx(x, components_data['total'], 'o-', linewidth=2.5, markersize=8,
                 color='blue', label=lang['total'])
    ax1.semilogx(x, components_data['total_left'], 's-', linewidth=2, markersize=6,
                 color='red', alpha=0.7, label=lang['left_side'])
    ax1.semilogx(x, components_data['total_right'], '^-', linewidth=2, markersize=6,
                 color='green', alpha=0.7, label=lang['right_side'])
    ax1.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    ax1.set_xlabel(param_name, fontsize=12)
    ax1.set_ylabel(lang['specific_thrust'], fontsize=12)
    ax1.set_title(f"{lang['branch']} {branch_id}: {lang['total']} тяга {lang['vs']} {param_name}",
                  fontsize=14, fontweight='bold')
    ax1.grid(True, which='both', linestyle=':', alpha=0.7)
    ax1.legend()
    grid_row += 1

    # График 2: Компоненты левой стороны (токи и конвекция)
    ax2 = fig.add_subplot(gs[grid_row, 0])
    ax2.semilogx(x, components_data['F_jmD_left'], 'o-', linewidth=1.5, markersize=5,
                 label=lang['magnetic_eddy'])
    ax2.semilogx(x, components_data['p_Bt_jeB_left'], 's-', linewidth=1.5, markersize=5,
                 label=lang['electric_eddy_B'])
    ax2.semilogx(x, components_data['p_I_conv_Ht_left'], '^-', linewidth=1.5, markersize=5,
                 label=lang['convective_I_tang'])
    ax2.semilogx(x, components_data['t_I_conv_Hn_left'], 'v-', linewidth=1.5, markersize=5,
                 label=lang['convective_I_norm'])
    ax2.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    ax2.set_xlabel(param_name, fontsize=11)
    ax2.set_ylabel(lang['specific_thrust'], fontsize=11)
    ax2.set_title(f"{lang['left_side']} - Токи и конвективные силы", fontsize=13)
    ax2.grid(True, which='both', linestyle=':', alpha=0.7)
    ax2.legend(fontsize=9)

    # График 3: Компоненты правой стороны (токи и конвекция)
    ax3 = fig.add_subplot(gs[grid_row, 1])
    ax3.semilogx(x, components_data['F_jmD_right'], 'o-', linewidth=1.5, markersize=5,
                 label=lang['magnetic_eddy'])
    ax3.semilogx(x, components_data['p_Bt_jeB_right'], 's-', linewidth=1.5, markersize=5,
                 label=lang['electric_eddy_B'])
    ax3.semilogx(x, components_data['p_I_conv_Ht_right'], '^-', linewidth=1.5, markersize=5,
                 label=lang['convective_I_tang'])
    ax3.semilogx(x, components_data['t_I_conv_Hn_right'], 'v-', linewidth=1.5, markersize=5,
                 label=lang['convective_I_norm'])
    ax3.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    ax3.set_xlabel(param_name, fontsize=11)
    ax3.set_ylabel(lang['specific_thrust'], fontsize=11)
    ax3.set_title(f"{lang['right_side']} - Токи и конвективные силы", fontsize=13)
    ax3.grid(True, which='both', linestyle=':', alpha=0.7)
    ax3.legend(fontsize=9)
    grid_row += 1

    # График 4: Роторные и динамические силы, натяжение полей (левая)
    ax4 = fig.add_subplot(gs[grid_row, 0])
    ax4.semilogx(x, components_data['F_I_rotB_left'], 'o-', linewidth=1.5, markersize=5,
                 label=lang['rotforce_I'])
    ax4.semilogx(x, components_data['F_P_rotE_left'], 's-', linewidth=1.5, markersize=5,
                 label=lang['rotforce_P'])
    ax4.semilogx(x, components_data['F_dyn_I_D_left'], 'o-', linewidth=1.5, markersize=5,
                 label=lang['dynamic_I'])
    ax4.semilogx(x, components_data['F_dyn_P_B_left'], 's-', linewidth=1.5, markersize=5,
                 label=lang['dynamic_P'])
    ax4.semilogx(x, components_data['t_Dn_left'], '^-', linewidth=1.5, markersize=5,
                 label=lang['tension_D'])
    ax4.semilogx(x, components_data['p_Et_left'], 'v-', linewidth=1.5, markersize=5,
                 label=lang['pressure_E'])
    ax4.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    ax4.set_xlabel(param_name, fontsize=11)
    ax4.set_ylabel(lang['specific_thrust'], fontsize=11)
    ax4.set_title(f"{lang['left_side']} - Роторные и динамические силы, натяжение полей", fontsize=13)
    ax4.grid(True, which='both', linestyle=':', alpha=0.7)
    ax4.legend(fontsize=9)

    # График 5: Роторные и динамические силы, натяжение полей (правая)
    ax5 = fig.add_subplot(gs[grid_row, 1])
    ax5.semilogx(x, components_data['F_I_rotB_right'], 'o-', linewidth=1.5, markersize=5,
                 label=lang['rotforce_I'])
    ax5.semilogx(x, components_data['F_P_rotE_right'], 's-', linewidth=1.5, markersize=5,
                 label=lang['rotforce_P'])
    ax5.semilogx(x, components_data['F_dyn_I_D_right'], 'o-', linewidth=1.5, markersize=5,
                 label=lang['dynamic_I'])
    ax5.semilogx(x, components_data['F_dyn_P_B_right'], 's-', linewidth=1.5, markersize=5,
                 label=lang['dynamic_P'])
    ax5.semilogx(x, components_data['t_Dn_right'], '^-', linewidth=1.5, markersize=5,
                 label=lang['tension_D'])
    ax5.semilogx(x, components_data['p_Et_right'], 'v-', linewidth=1.5, markersize=5,
                 label=lang['pressure_E'])
    ax5.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    ax5.set_xlabel(param_name, fontsize=11)
    ax5.set_ylabel(lang['specific_thrust'], fontsize=11)
    ax5.set_title(f"{lang['right_side']} - Роторные и динамические силы, натяжение полей", fontsize=13)
    ax5.grid(True, which='both', linestyle=':', alpha=0.7)
    ax5.legend(fontsize=9)
    grid_row += 1

    # График 6: Все компоненты левой стороны
    ax6 = fig.add_subplot(gs[grid_row, 0])
    ax6.semilogx(x, components_data['F_jmD_left'], 'o-', linewidth=1.2, markersize=4, label='F_jmD')
    ax6.semilogx(x, components_data['p_Bt_jeB_left'], 's-', linewidth=1.2, markersize=4, label='p_Bt_jeB')
    ax6.semilogx(x, components_data['p_I_conv_Ht_left'], '^-', linewidth=1.2, markersize=4, label='p_I_conv_Ht')
    ax6.semilogx(x, components_data['t_I_conv_Hn_left'], 'v-', linewidth=1.2, markersize=4, label='t_I_conv_Hn')
    ax6.semilogx(x, components_data['F_I_rotB_left'], 'd-', linewidth=1.2, markersize=4, label='F_I_rotB')
    ax6.semilogx(x, components_data['F_P_rotE_left'], 'p-', linewidth=1.2, markersize=4, label='F_P_rotE')
    ax6.semilogx(x, components_data['F_dyn_I_D_left'], 'd-', linewidth=1.2, markersize=4, label='F_dyn_I_D')
    ax6.semilogx(x, components_data['F_dyn_P_B_left'], 'p-', linewidth=1.2, markersize=4, label='F_dyn_P_B')
    ax6.semilogx(x, components_data['t_Dn_left'], '*-', linewidth=1.2, markersize=4, label='t_Dn')
    ax6.semilogx(x, components_data['p_Et_left'], 'h-', linewidth=1.2, markersize=4, label='p_Et')
    ax6.semilogx(x, components_data['total_left'], 'o-', linewidth=2.5, markersize=6,
                 color='black', label='TOTAL', zorder=10)
    ax6.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    ax6.set_xlabel(param_name, fontsize=11)
    ax6.set_ylabel(lang['specific_thrust'], fontsize=11)
    ax6.set_title(f"{lang['left_side']} - Все компоненты", fontsize=13)
    ax6.grid(True, which='both', linestyle=':', alpha=0.7)
    ax6.legend(fontsize=8, ncol=2)

    # График 7: Все компоненты правой стороны
    ax7 = fig.add_subplot(gs[grid_row, 1])
    ax7.semilogx(x, components_data['F_jmD_right'], 'o-', linewidth=1.2, markersize=4, label='F_jmD')
    ax7.semilogx(x, components_data['p_Bt_jeB_right'], 's-', linewidth=1.2, markersize=4, label='p_Bt_jeB')
    ax7.semilogx(x, components_data['p_I_conv_Ht_right'], '^-', linewidth=1.2, markersize=4, label='p_I_conv_Ht')
    ax7.semilogx(x, components_data['t_I_conv_Hn_right'], 'v-', linewidth=1.2, markersize=4, label='t_I_conv_Hn')
    ax7.semilogx(x, components_data['F_I_rotB_right'], 'd-', linewidth=1.2, markersize=4, label='F_I_rotB')
    ax7.semilogx(x, components_data['F_P_rotE_right'], 'p-', linewidth=1.2, markersize=4, label='F_P_rotE')
    ax7.semilogx(x, components_data['F_dyn_I_D_right'], 'd-', linewidth=1.2, markersize=4, label='F_dyn_I_D')
    ax7.semilogx(x, components_data['F_dyn_P_B_right'], 'p-', linewidth=1.2, markersize=4, label='F_dyn_P_B')
    ax7.semilogx(x, components_data['t_Dn_right'], '*-', linewidth=1.2, markersize=4, label='t_Dn')
    ax7.semilogx(x, components_data['p_Et_right'], 'h-', linewidth=1.2, markersize=4, label='p_Et')
    ax7.semilogx(x, components_data['total_right'], 'o-', linewidth=2.5, markersize=6,
                 color='black', label='TOTAL', zorder=10)
    ax7.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    ax7.set_xlabel(param_name, fontsize=11)
    ax7.set_ylabel(lang['specific_thrust'], fontsize=11)
    ax7.set_title(f"{lang['right_side']} - Все компоненты", fontsize=13)
    ax7.grid(True, which='both', linestyle=':', alpha=0.7)
    ax7.legend(fontsize=8, ncol=2)
    grid_row += 1

    # График 8: Относительный вклад компонент (левая)
    ax8 = fig.add_subplot(gs[grid_row, 0])
    components_left = [
        components_data['F_jmD_left'],
        components_data['p_Bt_jeB_left'],
        components_data['p_I_conv_Ht_left'],
        components_data['t_I_conv_Hn_left'],
        components_data['F_I_rotB_left'],
        components_data['F_P_rotE_left'],
        components_data['F_dyn_I_D_left'],
        components_data['F_dyn_P_B_left'],
        components_data['t_Dn_left'],
        components_data['p_Et_left']
    ]
    labels_short = ['F_jmD', 'p_Bt_jeB', 'p_I_conv_Ht', 't_I_conv_Hn',
                    'F_I_rotB', 'F_P_rotE', 'F_dyn_I_D', 'F_dyn_P_B', 't_Dn', 'p_Et']

    # Вычисляем относительный вклад в процентах
    for i, comp in enumerate(components_left):
        relative = [100 * c / t if abs(t) > 1e-10 else 0
                   for c, t in zip(comp, components_data['total_left'])]
        ax8.semilogx(x, relative, 'o-', linewidth=1.5, markersize=4, label=labels_short[i])

    ax8.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    ax8.set_xlabel(param_name, fontsize=11)
    ax8.set_ylabel('Вклад (%)', fontsize=11)
    ax8.set_title(f"{lang['left_side']} - Относительный вклад компонент", fontsize=13)
    ax8.grid(True, which='both', linestyle=':', alpha=0.7)
    ax8.legend(fontsize=8, ncol=2)

    # График 9: Относительный вклад компонент (правая)
    ax9 = fig.add_subplot(gs[4, 1])
    components_right = [
        components_data['F_jmD_right'],
        components_data['p_Bt_jeB_right'],
        components_data['p_I_conv_Ht_right'],
        components_data['t_I_conv_Hn_right'],
        components_data['F_I_rotB_right'],
        components_data['F_P_rotE_right'],
        components_data['F_dyn_I_D_right'],
        components_data['F_dyn_P_B_right'],
        components_data['t_Dn_right'],
        components_data['p_Et_right']
    ]

    for i, comp in enumerate(components_right):
        relative = [100 * c / t if abs(t) > 1e-10 else 0
                   for c, t in zip(comp, components_data['total_right'])]
        ax9.semilogx(x, relative, 'o-', linewidth=1.5, markersize=4, label=labels_short[i])

    ax9.axhline(y=0, color='gray', linestyle='--', linewidth=0.5)
    ax9.set_xlabel(param_name, fontsize=11)
    ax9.set_ylabel('Вклад (%)', fontsize=11)
    ax9.set_title(f"{lang['right_side']} - Относительный вклад компонент", fontsize=13)
    ax9.grid(True, which='both', linestyle=':', alpha=0.7)
    ax9.legend(fontsize=8, ncol=2)

    plt.suptitle(f"{lang['branch']} {branch_id} - {lang['components']}",
                 fontsize=16, fontweight='bold')
    # plt.show()
    plt.savefig(f"{figdir}/{figfilename}")
    plt.close()

def _add_component_data(components_data, result, x_val):
    """Вспомогательная функция для добавления данных компонент."""
    components_data['x'].append(x_val)
    components_data['total'].append(float(result['thrust_N_per_kW']))
    components_data['total_left'].append(float(result['thrust_N_per_kW_left']))
    components_data['total_right'].append(float(result['thrust_N_per_kW_right']))

    components = ['F_jmD',
                  'p_Bt_jeB',
                  'p_Ht_jeH',
                  'p_I_conv_Ht',
                  't_I_conv_Hn',
                  'F_I_rotB',
                  'F_P_rotE',
                  'F_dyn_I_D',
                  'F_dyn_P_B',
                  't_Dn',
                  'p_Et',
                 ]

    for component in components:
        components_data[component + '_left' ].append(float(result['N_per_kW_' + component + '_conj_left']))
        components_data[component + '_right'].append(float(result['N_per_kW_' + component + '_conj_right']))



def plot_field_report_for_branch(results, param_name, branch_id, base_digit_values,
        figfilename,
        distance_threshold=0.1,
        use_K_params=False,
        language='ru',
        show_appendix=False,
        appendix_step=5,
        max_points_to_show=5):
    """
    Построение отчёта по полям для конкретной ветви решения.

    Parameters:
    -----------
    results : list
        Список результатов из scan_parameter
    param_name : str
        Имя параметра сканирования
    branch_id : int
        ID ветви для построения
    base_digit_values : list
        Базовые параметры для расчёта полей
    distance_threshold : float
        Порог расстояния для сопоставления ветвей
    use_K_params : bool
        Использовать ли K-параметры для сопоставления
    language : str
        'ru' или 'en'
    show_appendix : bool
        Показывать ли промежуточные точки на одном графике
    appendix_step : int
        Показывать каждую N-ю точку из ветви
    max_points_to_show : int
        Максимальное количество точек ветви для отображения
    """

    translations = {
        'ru': {
            'field_report': 'Отчёт по полям для ветви',
            'wave_decay': 'Затухание E(z)',
            'normal_B': 'Нормальная B',
            'tangential_E': 'Тангенциальная E',
            'tangential_H': 'Тангенциальная H',
            'real_part': 'Re',
            'imag_part': 'Im',
            'amplitude': 'Abs',
            'branch': 'Ветвь',
            'point': 'Точка',
            'thrust': 'Тяга',
            'position': 'x (см)',
            'decay_length': 'длина затухания',
            'param_value': 'Значение параметра'
        },
        'en': {
            'field_report': 'Field Report for Branch',
            'wave_decay': 'E(z) decay',
            'normal_B': 'Normal B',
            'tangential_E': 'Tangential E',
            'tangential_H': 'Tangential H',
            'real_part': 'Re',
            'imag_part': 'Im',
            'amplitude': 'Abs',
            'branch': 'Branch',
            'point': 'Point',
            'thrust': 'Thrust',
            'position': 'x (cm)',
            'decay_length': 'decay length',
            'param_value': 'Parameter value'
        }
    }

    lang = translations.get(language, translations['ru'])

    # Получаем точки ветви
    matched_branches = match_branches_between_intervals(
        results, param_name, distance_threshold, use_K_params
    )

    if branch_id not in matched_branches:
        print(f"❌ Ветвь {branch_id} не найдена")
        return

    branch_points = matched_branches[branch_id]

    if len(branch_points) == 0:
        print(f"❌ Ветвь {branch_id} не содержит точек")
        return

    print(f"\n{'='*80}")
    print(f"{lang['field_report']} {branch_id}")
    print(f"Параметр: {param_name}")
    print(f"Количество точек в ветви: {len(branch_points)}")
    print(f"{'='*80}\n")

    import matplotlib.pyplot as plt
    import numpy as np

    # Определяем тип волны (TM или TE)
    try:
        wave_type = 'TM' if is_TM else 'TE'
    except:
        wave_type = 'TM'

    # Функция для вычисления полей в точке
    def compute_fields_for_point(kz_val, sz_val, params_dict, dv_local,
                                 sign_K_H_l_d,
                                 sign_K_E_l_d,
                                 sign_K_H_r_d,
                                 sign_K_E_r_d,
                                ):
        """Вычисляет поля для заданной точки"""
        k_z_cur = [kz == kz_val, sz == sz_val, k_z == kz_val + I * sz_val]

        dv = override_params(dv_local, params_dict)

        K_v = [K_v_s.subs(k_subs).subs(kappa_vacuum_sol).subs(k_z_cur).subs(dv)
               for K_v_s in K_vacuum_subs]
        K_l = [K_l_s.subs(k_subs).subs(kappa_vacuum_sol).subs(k_z_cur).subs(dv)
               for K_l_s in K_left_conductor_subs]
        K_r = [K_r_s.subs(k_subs).subs(kappa_vacuum_sol).subs(k_z_cur).subs(dv)
               for K_r_s in K_right_conductor_subs]

        K_l = fix_K_sign(K_l)
        K_r = fix_K_sign(K_r)

        M4_cur = M4_K \
            .subs(K_v).subs(K_l).subs(K_r) \
            .subs(kappa_vacuum_sol).subs(epsilon_mu_subs) \
            .subs(k_subs).subs(dv).subs(k_z_cur) \
            .subs(sign_K_1_subs)

        # Mb_cur = M_boundaries_K \
        #     .subs(K_v).subs(K_l).subs(K_r) \
        #     .subs(kappa_vacuum_sol).subs(epsilon_mu_subs) \
        #     .subs(k_subs).subs(dv).subs(k_z_cur) \
        #     .subs(sign_K_1_subs)

        X4, residual = svd_X4(M4_cur)
        coeffs = [vars4[i] == X4[i] for i in range(len(vars4))]

        # boundaries_res = (Mb_cur*(Matrix(vars4).transpose())).subs(coeffs).n()
        # print("boundaries_res", boundaries_res)

        fields = subs_fields(k_z_cur, K_l, K_v, K_r, coeffs, dv)
        return fields

    # Функция для вычисления значений поля на массиве x
    def eval_field_array_real(field_expr, x_array):
        """Вычисляет значения поля для массива x"""
        values = []
        for x_val in x_array:
            val = field_expr.subs(x==x_val)
            values.append(val.real())
        return values

    def eval_field_array_imag(field_expr, x_array):
        """Вычисляет значения поля для массива x"""
        values = []
        for x_val in x_array:
            val = field_expr.subs(x==x_val)
            values.append(val.imag())
        return values

    # Находим соответствующие результаты для точек ветви
    points_with_results = []
    for bp in branch_points:
        # Ищем result для этой точки
        for results_list in results:
            if results_list:
                for result in results_list:
                    if (abs(result['kz'] - bp['kz']) < 1e-6 and
                        abs(result['sz'] - bp['sz']) < 1e-6):
                        points_with_results.append({
                            'branch_point': bp,
                            'result': result
                        })
                        break

    if len(points_with_results) == 0:
        print(f"❌ Не найдено результатов для точек ветви {branch_id}")
        return

    # Сортируем точки по значению параметра
    points_with_results.sort(key=lambda p: p['branch_point']['x'])

    # Выбираем точки для отображения (равномерно распределённые)
    if len(points_with_results) > max_points_to_show:
        indices = np.linspace(0, len(points_with_results)-1, max_points_to_show, dtype=int)
        selected_points = [points_with_results[i] for i in indices]
    else:
        selected_points = points_with_results[::max(1, len(points_with_results)//max_points_to_show)]

    # Получаем A и толщины из первой точки
    first_result = selected_points[0]['result']
    param_overrides = first_result['params']
    dv = override_params(base_digit_values, param_overrides)

    A_val = None
    h_l_val = h_conductor_l
    h_r_val = h_conductor_r
    for eq in dv:
        if str(eq.lhs()) == 'a':
            A_val = float(eq.rhs())

    if A_val is None:
        A_val = A

    # Создаём массивы x для каждой области
    n_points = 100
    x_left   = np.linspace(-A_val-h_l_val, -A_val, n_points)
    x_vacuum = np.linspace(-A_val, A_val, n_points)
    x_right  = np.linspace( A_val, A_val+h_r_val, n_points)

    # Цвета для разных точек
    colors = plt.cm.viridis(np.linspace(0, 1, len(selected_points)))

    # Подготовка точек для построения
    points_to_plot = []
    for i, pw in enumerate(selected_points):
        bp = pw['branch_point']
        result = pw['result']

        points_to_plot.append({
            'kz': float(bp['kz']),
            'sz': float(bp['sz']),
            'params': result['params'],
            'dv': override_params(base_digit_values, result['params']),
            'label': f"{param_name}={bp['x']:.2e}",
            'color': colors[i],
            'param_value': bp['x'],
            'thrust': bp['thrust']
        })

    # Вычисляем поля для всех точек
    print(f"Вычисление полей для {len(points_to_plot)} точек...")
    all_fields = []
    for point in points_to_plot:
        fields = compute_fields_for_point(point['kz'], point['sz'],
                                          point['params'], point['dv'],
                                          sign_K_H_l_d = 1,
                                          sign_K_E_l_d = 1,
                                          sign_K_H_r_d = 1,
                                          sign_K_E_r_d = 1,
                                         )
        all_fields.append(fields)

    # Создаём сетку графиков: 3 строки × 3 колонки
    fig = plt.figure(figsize=(18, 12))
    gs = fig.add_gridspec(3, 3, hspace=0.35, wspace=0.3)

    # Заголовок
    fig.suptitle(f"{lang['field_report']} {branch_id}: {param_name}",
                 fontsize=14, fontweight='bold')

    # Строка 0: График затухания E(z) и нормальная компонента Bx
    # График 0,0: Затухание вдоль z для первой точки
    ax = fig.add_subplot(gs[0, 0])
    kz_sol = points_to_plot[0]['kz']
    sz_sol = points_to_plot[0]['sz']
    kz_re = float(abs(kz_sol))
    kz_im = float(abs(sz_sol))

    if kz_im > 0:
        decay_length = 5.0 / kz_im
    else:
        decay_length = 10.0

    z_array = np.linspace(0, decay_length, 500).astype(float)
    field_z = np.exp(-kz_im * z_array) * np.cos(kz_re * z_array)
    ax.plot(z_array, field_z, 'b-', linewidth=2)
    ax.set_xlabel('z (см)')
    ax.set_ylabel('E(z)')
    ax.set_title(f"{lang['wave_decay']}\n{lang['decay_length']}={decay_length:.2f} см")
    ax.grid(True, alpha=0.3)

    # График 0,1: Нормальная компонента Bx - Re
    ax = fig.add_subplot(gs[0, 1])
    for i, (point, fields) in enumerate(zip(points_to_plot, all_fields)):
        Bx_l_d = fields[27]  # Bx_l
        Bx_r_d = fields[28]  # Bx_r
        Bx_v_d = fields[29]  # Bx_v

        # print("Bx_l_d", Bx_l_d)
        # print("Bx_v_d", Bx_v_d)
        # print("Bx_r_d", Bx_r_d)

        Bx_l_vals_real = eval_field_array_real(Bx_l_d, x_left)
        Bx_v_vals_real = eval_field_array_real(Bx_v_d, x_vacuum)
        Bx_r_vals_real = eval_field_array_real(Bx_r_d, x_right)

        # print("Bx_l_vals", Bx_l_vals)

        alpha = 1.0 if i == 0 else 0.5
        lw = 2 if i == 0 else 1

        ax.plot(x_left,   Bx_l_vals_real, color=point['color'], alpha=alpha, linewidth=lw)
        ax.plot(x_vacuum, Bx_v_vals_real, color=point['color'], alpha=alpha, linewidth=lw)
        ax.plot(x_right,  Bx_r_vals_real, color=point['color'], alpha=alpha, linewidth=lw,
                label=point['label'])

    ax.axvline(x=-A_val, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x= A_val, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel(lang['position'])
    ax.set_ylabel(f"{lang['normal_B']} - {lang['real_part']}")
    ax.set_title(f"{lang['normal_B']} - {lang['real_part']}")
    ax.grid(True, alpha=0.3)
    ax.legend(fontsize=7, loc='best')

    # График 0,2: Нормальная компонента Bx - Im
    ax = fig.add_subplot(gs[0, 2])
    for i, (point, fields) in enumerate(zip(points_to_plot, all_fields)):
        Bx_l_d = fields[27]
        Bx_r_d = fields[28]
        Bx_v_d = fields[29]

        Bx_l_vals_imag = eval_field_array_imag(Bx_l_d, x_left)
        Bx_v_vals_imag = eval_field_array_imag(Bx_v_d, x_vacuum)
        Bx_r_vals_imag = eval_field_array_imag(Bx_r_d, x_right)

        alpha = 1.0 if i == 0 else 0.5
        lw = 2 if i == 0 else 1

        ax.plot(x_left,   Bx_l_vals_imag, color=point['color'], alpha=alpha, linewidth=lw)
        ax.plot(x_vacuum, Bx_v_vals_imag, color=point['color'], alpha=alpha, linewidth=lw)
        ax.plot(x_right,  Bx_r_vals_imag, color=point['color'], alpha=alpha, linewidth=lw)

    ax.axvline(x=-A_val, color='gray', linestyle='--', alpha=0.5)
    ax.axvline(x=A_val, color='gray', linestyle='--', alpha=0.5)
    ax.set_xlabel(lang['position'])
    ax.set_ylabel(f"{lang['normal_B']} - {lang['imag_part']}")
    ax.set_title(f"{lang['normal_B']} - {lang['imag_part']}")
    ax.grid(True, alpha=0.3)

    # Строка 1: Тангенциальная компонента E
    if wave_type == 'TM':
        E_indices = (6, 7, 8)  # Ez_l, Ez_r, Ez_v
    else:
        E_indices = (3, 4, 5)  # Ey_l, Ey_r, Ey_v

    for col, (component_type, title) in enumerate([('real', lang['real_part']),
                                                     ('imag', lang['imag_part']),
                                                     ('abs', lang['amplitude'])]):
        ax = fig.add_subplot(gs[1, col])

        for i, (point, fields) in enumerate(zip(points_to_plot, all_fields)):
            E_l_d = fields[E_indices[0]]
            E_r_d = fields[E_indices[1]]
            E_v_d = fields[E_indices[2]]

            E_l_vals_real = eval_field_array_real(E_l_d, x_left)
            E_v_vals_real = eval_field_array_real(E_v_d, x_vacuum)
            E_r_vals_real = eval_field_array_real(E_r_d, x_right)

            E_l_vals_imag = eval_field_array_imag(E_l_d, x_left)
            E_v_vals_imag = eval_field_array_imag(E_v_d, x_vacuum)
            E_r_vals_imag = eval_field_array_imag(E_r_d, x_right)

            alpha = 1.0 if i == 0 else 0.5
            lw = 2 if i == 0 else 1

            if component_type == 'real':
                vals_l, vals_v, vals_r = E_l_vals_real, E_v_vals_real, E_r_vals_real
            elif component_type == 'imag':
                vals_l, vals_v, vals_r = E_l_vals_imag, E_v_vals_imag, E_r_vals_imag
#             else:  # abs
#                 vals_l, vals_v, vals_r = np.abs(E_l_vals), np.abs(E_v_vals), np.abs(E_r_vals)

            ax.plot(x_left, vals_l, color=point['color'], alpha=alpha, linewidth=lw)
            ax.plot(x_vacuum, vals_v, color=point['color'], alpha=alpha, linewidth=lw)
            ax.plot(x_right, vals_r, color=point['color'], alpha=alpha, linewidth=lw)

        ax.axvline(x=-A_val, color='gray', linestyle='--', alpha=0.5)
        ax.axvline(x=A_val, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlabel(lang['position'])
        ax.set_ylabel(f"{lang['tangential_E']} - {title}")
        ax.set_title(f"{lang['tangential_E']} - {title}")
        ax.grid(True, alpha=0.3)

    # Строка 2: Тангенциальная компонента H
    if wave_type == 'TM':
        H_indices = (21, 22, 23)  # Hy_l, Hy_r, Hy_v
    else:
        H_indices = (24, 25, 26)  # Hz_l, Hz_r, Hz_v

    for col, (component_type, title) in enumerate([('real', lang['real_part']),
                                                     ('imag', lang['imag_part']),
                                                     ('abs', lang['amplitude'])]):
        ax = fig.add_subplot(gs[2, col])

        for i, (point, fields) in enumerate(zip(points_to_plot, all_fields)):
            H_l_d = fields[H_indices[0]]
            H_r_d = fields[H_indices[1]]
            H_v_d = fields[H_indices[2]]

            H_l_vals_real = eval_field_array_real(H_l_d, x_left)
            H_v_vals_real = eval_field_array_real(H_v_d, x_vacuum)
            H_r_vals_real = eval_field_array_real(H_r_d, x_right)

            H_l_vals_imag = eval_field_array_imag(H_l_d, x_left)
            H_v_vals_imag = eval_field_array_imag(H_v_d, x_vacuum)
            H_r_vals_imag = eval_field_array_imag(H_r_d, x_right)

            alpha = 1.0 if i == 0 else 0.5
            lw = 2 if i == 0 else 1

            if component_type == 'real':
                vals_l, vals_v, vals_r = H_l_vals_real, H_v_vals_real, H_r_vals_real
            elif component_type == 'imag':
                vals_l, vals_v, vals_r = H_l_vals_imag, H_v_vals_imag, H_r_vals_imag
#             else:  # abs
#                 vals_l, vals_v, vals_r = np.abs(H_l_vals), np.abs(H_v_vals), np.abs(H_r_vals)

            ax.plot(x_left, vals_l, color=point['color'], alpha=alpha, linewidth=lw)
            ax.plot(x_vacuum, vals_v, color=point['color'], alpha=alpha, linewidth=lw)
            ax.plot(x_right, vals_r, color=point['color'], alpha=alpha, linewidth=lw)

        ax.axvline(x=-A_val, color='gray', linestyle='--', alpha=0.5)
        ax.axvline(x=A_val, color='gray', linestyle='--', alpha=0.5)
        ax.set_xlabel(lang['position'])
        ax.set_ylabel(f"{lang['tangential_H']} - {title}")
        ax.set_title(f"{lang['tangential_H']} - {title}")
        ax.grid(True, alpha=0.3)

    # plt.show()
    plt.savefig(f"{figdir}/{figfilename}")
    plt.close()


def plot_tensor_report_for_branch(results, param_name, branch_id, base_digit_values,
        figfilename, stress_tensor_name = 'maxwell_stress_tensor',
        distance_threshold=0.1,
        use_K_params=False,
        language='ru',
        show_appendix=False,
        appendix_step=5,
        max_points_to_show=5):
    """
    Построение отчёта по тензору Максвелла для конкретной ветви решения.
    Визуализирует дивергенцию тензора натяжений div T_x в трёх областях.

    Parameters:
    -----------
    results : list
        Список результатов из scan_parameter
    param_name : str
        Имя параметра сканирования
    branch_id : int
        ID ветви для построения
    base_digit_values : list
        Базовые параметры для расчёта полей
    distance_threshold : float
        Порог расстояния для сопоставления ветвей
    use_K_params : bool
        Использовать ли K-параметры для сопоставления
    language : str
        'ru' или 'en'
    show_appendix : bool
        Показывать ли промежуточные точки на одном графике
    appendix_step : int
        Показывать каждую N-ю точку из ветви
    max_points_to_show : int
        Максимальное количество точек ветви для отображения
    """

    translations = {
        'ru': {
            'tensor_report': f"Отчёт по тензору {stress_tensor_name} для ветви",
            'divergence': 'Дивергенция div T_x',
            'left_conductor': 'Левый проводник (феррит)',
            'vacuum': 'Вакуумный зазор',
            'right_conductor': 'Правый проводник (металл)',
            'branch': 'Ветвь',
            'point': 'Точка',
            'thrust': 'Тяга',
            'position': 'x (см)',
            'force_density': 'Плотность силы (дин/см³)',
            'param_value': 'Значение параметра',
            'total_force': 'Интегральная сила',
            'left_force': 'Сила слева',
            'right_force': 'Сила справа',
            'real_part': 'Действительная часть',
            'imag_part': 'Действительная часть',
        },
        'en': {
            'tensor_report': f"{stress_tensor_name} Stress Tensor Report for Branch",
            'divergence': 'Divergence div T_x',
            'left_conductor': 'Left Conductor (Ferrite)',
            'vacuum': 'Vacuum Gap',
            'right_conductor': 'Right Conductor (Metal)',
            'branch': 'Branch',
            'point': 'Point',
            'thrust': 'Thrust',
            'position': 'x (cm)',
            'force_density': 'Force Density (dyn/cm³)',
            'param_value': 'Parameter value',
            'total_force': 'Total Force',
            'left_force': 'Left Force',
            'right_force': 'Right Force',
            'real_part': 'Real part',
            'imag_part': 'Image part',
        }
    }

    lang = translations.get(language, translations['ru'])

    # Получаем точки ветви
    matched_branches = match_branches_between_intervals(
        results, param_name, distance_threshold, use_K_params
    )

    if branch_id not in matched_branches:
        print(f"❌ Ветвь {branch_id} не найдена")
        return

    branch_points = matched_branches[branch_id]

    if len(branch_points) == 0:
        print(f"❌ Ветвь {branch_id} не содержит точек")
        return

    print(f"\n{'='*80}")
    print(f"{lang['tensor_report']} {branch_id}")
    print(f"Параметр: {param_name}")
    print(f"Количество точек в ветви: {len(branch_points)}")
    print(f"{'='*80}\n")

    import matplotlib.pyplot as plt
    import numpy as np

    # Функция для вычисления значений дивергенции на массиве x
    def eval_divergence_array(div_expr, x_array):
        """Вычисляет значения дивергенции для массива x"""
        values = []
        for x_val in x_array:
            # Подставляем x, kz, sz и параметры
            val_complex = div_expr.subs(x == x_val)
            values.append(val_complex)
        return values

    # Находим соответствующие результаты для точек ветви
    points_with_results = []
    for bp in branch_points:
        # Ищем result для этой точки
        for results_list in results:
            if results_list:
                for result in results_list:
                    if (abs(result['kz'] - bp['kz']) < 1e-6 and
                        abs(result['sz'] - bp['sz']) < 1e-6):
                        points_with_results.append({
                            'branch_point': bp,
                            'result': result
                        })
                        break

    if len(points_with_results) == 0:
        print(f"❌ Не найдено результатов для точек ветви {branch_id}")
        return

    # Сортируем точки по значению параметра
    points_with_results.sort(key=lambda p: p['branch_point']['x'])

    # Выбираем точки для отображения (равномерно распределённые)
    if len(points_with_results) > max_points_to_show:
        indices = np.linspace(0, len(points_with_results)-1, max_points_to_show, dtype=int)
        selected_points = [points_with_results[i] for i in indices]
    else:
        selected_points = points_with_results[::max(1, len(points_with_results)//max_points_to_show)]

    # Получаем A и толщины из первой точки
    first_result = selected_points[0]['result']
    param_overrides = first_result['params']
    dv = override_params(base_digit_values, param_overrides)

    A_val = None
    h_l_val = h_conductor_l
    h_r_val = h_conductor_r
    for eq in dv:
        if str(eq.lhs()) == 'a':
            A_val = float(eq.rhs())

    if A_val is None:
        A_val = A

    # Создаём массивы x для каждой области
    n_points = 100
    x_left = np.linspace(-A_val-h_l_val, -A_val, n_points)
    x_vacuum = np.linspace(-A_val, A_val, n_points)
    x_right = np.linspace(A_val, A_val+h_r_val, n_points)

    # Цвета для разных точек
    colors = plt.cm.viridis(np.linspace(0, 1, len(selected_points)))

    # Подготовка точек для построения
    points_to_plot = []
    for i, pw in enumerate(selected_points):
        bp = pw['branch_point']
        result = pw['result']

        points_to_plot.append({
            'kz': float(bp['kz']),
            'sz': float(bp['sz']),
            'params': result['params'],
            'dv': override_params(base_digit_values, result['params']),
            'label': f"{param_name}={bp['x']:.2e}",
            'color': colors[i],
            'param_value': bp['x'],
            'thrust': bp['thrust'],
            stress_tensor_name: result [stress_tensor_name]
        })

    # Вычисляем дивергенцию для всех точек
    print(f"Вычисление дивергенции тензора для {len(points_to_plot)} точек...")
    all_divergences = []

    for point in points_to_plot:
        kz_val = point['kz']
        sz_val = point['sz']
        dv_point = point['dv']

        div_l = point[stress_tensor_name]['div_T_x_l']
        div_v = point[stress_tensor_name]['div_T_x_v']
        div_r = point[stress_tensor_name]['div_T_x_r']

        if div_l is None or div_v is None or div_r is None:
            print(f"  ⚠️  Предупреждение: дивергенция не вычислена для {point['label']}")
            continue

        # Вычисляем значения на сетках
        div_l_vals = eval_divergence_array(div_l, x_left)
        div_v_vals = eval_divergence_array(div_v, x_vacuum)
        div_r_vals = eval_divergence_array(div_r, x_right)

        all_divergences.append({
            'x_left': x_left,
            'x_vacuum': x_vacuum,
            'x_right': x_right,
            'div_l': div_l_vals,
            'div_v': div_v_vals,
            'div_r': div_r_vals,
            'color': point['color'],
            'label': point['label'],
            'param_value': point['param_value'],
            'thrust': point['thrust']
        })

    if len(all_divergences) == 0:
        print("❌ Не удалось вычислить ни одной точки дивергенции")
        return

    # Создаём сетку графиков
    fig = plt.figure(figsize=(18, 14))
    gs = fig.add_gridspec(3, 2, hspace=0.4, wspace=0.3)

    # Заголовок
    fig.suptitle(f"{lang['tensor_report']} {branch_id}: {param_name}",
                 fontsize=14, fontweight='bold')

    # График 1: Дивергенция - действительная часть (все области)
    ax1 = fig.add_subplot(gs[0, 0])
    for idx, div_data in enumerate(all_divergences):
        is_first = (idx == 0)
        alpha = 1.0 if is_first else 0.5
        lw = 2 if is_first else 1

        # Объединяем все области для непрерывного графика
        x_all = np.concatenate([div_data['x_left'], div_data['x_vacuum'], div_data['x_right']])

        # Для SageMath объектов нужно явно извлекать действительную часть
        div_l_real = [val.n().real() for val in div_data['div_l']]
        div_v_real = [val.n().real() for val in div_data['div_v']]
        div_r_real = [val.n().real() for val in div_data['div_r']]

        div_all = div_l_real + div_v_real + div_r_real

        # Преобразуем в Python float для matplotlib
        x_all_float = [float(x) for x in x_all]
        div_all_float = [float(d) for d in div_all]

        ax1.plot(x_all_float, div_all_float, color=div_data['color'], alpha=alpha, linewidth=lw,
                 label=div_data['label'])

        # Границы
        ax1.axvline(x=-A_val, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
        ax1.axvline(x=A_val, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)

    ax1.set_xlabel(lang['position'], fontsize=12)
    ax1.set_ylabel(f"{lang['force_density']} - Re", fontsize=12)
    ax1.set_title(f"{lang['divergence']} - {lang['real_part']}", fontsize=13)
    ax1.grid(True, which='both', linestyle=':', alpha=0.5)
    ax1.legend(fontsize=8, loc='best')

    # График 2: Дивергенция - мнимая часть (все области)
    ax2 = fig.add_subplot(gs[0, 1])
    for idx, div_data in enumerate(all_divergences):
        is_first = (idx == 0)
        alpha = 1.0 if is_first else 0.5
        lw = 2 if is_first else 1

        # Объединяем все области для непрерывного графика
        x_all = np.concatenate([div_data['x_left'], div_data['x_vacuum'], div_data['x_right']])

        # Для SageMath объектов нужно явно извлекать действительную часть
        div_l_imag = [val.n().imag() for val in div_data['div_l']]
        div_v_imag = [val.n().imag() for val in div_data['div_v']]
        div_r_imag = [val.n().imag() for val in div_data['div_r']]

        div_all = div_l_imag + div_v_imag + div_r_imag

        # Преобразуем в Python float для matplotlib
        x_all_float = [float(x) for x in x_all]
        div_all_float = [float(d) for d in div_all]

        ax2.plot(x_all_float, div_all_float, color=div_data['color'], alpha=alpha, linewidth=lw,
                 label=div_data['label'])

        ax2.axvline(x=-A_val, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
        ax2.axvline(x=A_val, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)

    ax2.set_xlabel(lang['position'], fontsize=12)
    ax2.set_ylabel(f"{lang['force_density']} - Im", fontsize=12)
    ax2.set_title(f"{lang['divergence']} - {lang['imag_part']}", fontsize=13)
    ax2.grid(True, which='both', linestyle=':', alpha=0.5)
    ax2.legend(fontsize=8, loc='best')

    # График 3: Модуль дивергенции (логарифмическая шкала)
    ax3 = fig.add_subplot(gs[1, :])
    for idx, div_data in enumerate(all_divergences):
        is_first = (idx == 0)
        alpha = 1.0 if is_first else 0.5
        lw = 2 if is_first else 1

        # Объединяем все области для непрерывного графика
        x_all = np.concatenate([div_data['x_left'], div_data['x_vacuum'], div_data['x_right']])

        # Для SageMath объектов нужно явно извлекать действительную часть
        div_l_real = [val.n().abs() for val in div_data['div_l']]
        div_v_real = [val.n().abs() for val in div_data['div_v']]
        div_r_real = [val.n().abs() for val in div_data['div_r']]

        div_all = div_l_real + div_v_real + div_r_real

        # Преобразуем в Python float для matplotlib
        x_all_float = [float(x) for x in x_all]
        div_all_float = [float(d) for d in div_all]

        ax3.plot(x_all_float, div_all_float, color=div_data['color'], alpha=alpha, linewidth=lw,
                 label=div_data['label'])

        ax3.axvline(x=-A_val, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
        ax3.axvline(x=A_val, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)

    ax3.set_xlabel(lang['position'], fontsize=12)
    ax3.set_ylabel(f"|{lang['force_density']}| (log scale)", fontsize=12)
    ax3.set_title(f"{lang['divergence']} - Absolute Value (Log Scale)", fontsize=13)
    ax3.grid(True, which='both', linestyle=':', alpha=0.5)
    ax3.legend(fontsize=8, loc='best')

    # Вспомогательная функция: правило трапеций для списка значений Sage
    def trapz_sage(values, dx):
        """
        Численное интегрирование по правилу трапеций для списка комплексных чисел Sage.

        Параметры:
        ----------
        values : list
            Список комплексных чисел Sage (результат вычисления дивергенции)
        dx : float
            Шаг сетки по координате

        Возвращает:
        -----------
        float : интеграл от действительной части значений
        """
        if len(values) < 2:
            return 0.0

        # Извлекаем действительные части и конвертируем в обычные float
        vals_float = [float(val.real()) for val in values]

        # Правило трапеций для равномерной сетки:
        # ∫f(x)dx ≈ dx * (0.5*f₀ + f₁ + f₂ + ... + fₙ₋₂ + 0.5*fₙ₋₁)
        result = dx * (0.5*vals_float[0] + sum(vals_float[1:-1]) + 0.5*vals_float[-1])

        return result

    # График 4: Интегральная сила по областям
    ax4 = fig.add_subplot(gs[2, 0])

    x_params = []
    thrust_values = []
    force_l_values = []
    force_v_values = []
    force_r_values = []

    for div_data in all_divergences:
        x_params.append(div_data['param_value'])
        thrust_values.append(div_data['thrust'])

        # Численное интегрирование дивергенции (сила = интеграл от плотности силы)
        dx_l = div_data['x_left'][1] - div_data['x_left'][0]
        dx_v = div_data['x_vacuum'][1] - div_data['x_vacuum'][0]
        dx_r = div_data['x_right'][1] - div_data['x_right'][0]

        # Интеграл от действительной части дивергенции в каждой области
        force_l = trapz_sage(div_data['div_l'], dx_l)
        force_v = trapz_sage(div_data['div_v'], dx_v)
        force_r = trapz_sage(div_data['div_r'], dx_r)

        force_l_values.append(force_l)
        force_v_values.append(force_v)
        force_r_values.append(force_r)

    ax4.semilogx(x_params, force_l_values, 'o-', linewidth=2, markersize=8,
                 color='red', alpha=0.7, label=lang['left_force'])
    ax4.semilogx(x_params, force_v_values, 's-', linewidth=2, markersize=8,
                 color='blue', alpha=0.7, label=lang['vacuum'])
    ax4.semilogx(x_params, force_r_values, '^-', linewidth=2, markersize=8,
                 color='green', alpha=0.7, label=lang['right_force'])
    ax4.semilogx(x_params, thrust_values, 'd-', linewidth=2, markersize=8,
                 color='black', alpha=0.9, label=lang['total_force'])

    ax4.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)
    ax4.set_xlabel(param_name, fontsize=12)
    ax4.set_ylabel(lang['force_density'], fontsize=12)
    ax4.set_title(f"{lang['total_force']} by Region", fontsize=13)
    ax4.grid(True, which='both', linestyle=':', alpha=0.5)
    ax4.legend(fontsize=10, loc='best')

    # График 5: Сравнение тяги из тензора и полной тяги
    ax5 = fig.add_subplot(gs[2, 1])

    # Суммарная сила из дивергенции
    total_force_tensor = [fl + fv + fr for fl, fv, fr in zip(force_l_values,
                                                               force_v_values,
                                                               force_r_values)]

    ax5.semilogx(x_params, thrust_values, 'o-', linewidth=2, markersize=8,
                 color='blue', alpha=0.8, label=lang['total_force'] + ' (calc_thrust)')
    ax5.semilogx(x_params, total_force_tensor, 's-', linewidth=2, markersize=8,
                 color='red', alpha=0.8, label=lang['total_force'] + ' (∫div T dV)')

    ax5.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)
    ax5.set_xlabel(param_name, fontsize=12)
    ax5.set_ylabel(lang['force_density'], fontsize=12)
    ax5.set_title("Force Comparison: Direct vs Tensor Method", fontsize=13)
    ax5.grid(True, which='both', linestyle=':', alpha=0.5)
    ax5.legend(fontsize=10, loc='best')

    # plt.show()
    plt.savefig(f"{figdir}/{figfilename}")
    plt.close()


def plot_tensor_report_for_branch_2(results, param_name, branch_id, base_digit_values,
        figfilename, stress_tensor_name = 'maxwell_stress_tensor',
        distance_threshold=0.1, use_K_params=False,
        language='ru', max_points_to_show=5):
    """
    Построение отчёта по тензору Максвелла для конкретной ветви решения.
    Визуализирует компоненты тензора T_xx, T_xy, T_xz и их дивергенцию.
    """

    translations = {
        'ru': {
            'tensor_report': f"Отчёт по тензору {stress_tensor_name} для ветви",
            'divergence': 'Дивергенция div T_x',
            'T_xx': 'Нормальное напряжение T_xx',
            'T_xy': 'Сдвиговое напряжение T_xy',
            'T_xz': 'Сдвиговое напряжение T_xz',
            'left_conductor': 'Левый проводник (феррит)',
            'vacuum': 'Вакуумный зазор',
            'right_conductor': 'Правый проводник (металл)',
            'branch': 'Ветвь',
            'thrust': 'Тяга',
            'position': 'x (см)',
            'stress': 'Напряжение (дин/см²)',
            'force_density': 'Плотность силы (дин/см³)',
            'param_value': 'Значение параметра',
            'total_force': 'Интегральная сила',
            'surface_force': 'Поверхностная сила',
            'volume_force': 'Объёмная сила',
            'real_part': 'Re',
            'imag_part': 'Im'
        },
        'en': {
            'tensor_report': f"{stress_tensor_name} Stress Tensor Report for Branch",
            'divergence': 'Divergence div T_x',
            'T_xx': 'Normal Stress T_xx',
            'T_xy': 'Shear Stress T_xy',
            'T_xz': 'Shear Stress T_xz',
            'left_conductor': 'Left Conductor (Ferrite)',
            'vacuum': 'Vacuum Gap',
            'right_conductor': 'Right Conductor (Metal)',
            'branch': 'Branch',
            'thrust': 'Thrust',
            'position': 'x (cm)',
            'stress': 'Stress (dyn/cm²)',
            'force_density': 'Force Density (dyn/cm³)',
            'param_value': 'Parameter value',
            'total_force': 'Total Force',
            'surface_force': 'Surface Force',
            'volume_force': 'Volume Force',
            'real_part': 'Re',
            'imag_part': 'Im'
        }
    }

    lang = translations.get(language, translations['ru'])

    # Получаем точки ветви
    matched_branches = match_branches_between_intervals(
        results, param_name, distance_threshold, use_K_params
    )

    if branch_id not in matched_branches:
        print(f"❌ Ветвь {branch_id} не найдена")
        return

    branch_points = matched_branches[branch_id]
    if len(branch_points) == 0:
        print(f"❌ Ветвь {branch_id} не содержит точек")
        return

    print(f"\n{'='*80}")
    print(f"{lang['tensor_report']} {branch_id}")
    print(f"Параметр: {param_name}")
    print(f"Количество точек в ветви: {len(branch_points)}")
    print(f"{'='*80}\n")

    import matplotlib.pyplot as plt

    # Функция для вычисления значений тензора на массиве x
    def eval_tensor_component(tensor_expr, x_array, dv_local):
        """Вычисляет значения компоненты тензора для массива x"""
        values = []

        for x_val in x_array:
            val = tensor_expr.subs(x == x_val).subs(dv_local)
            # Конвертируем в комплексное число Python
            val_num = val.n().real()
            values.append(val_num)

        return values

    # Находим результаты для точек ветви
    points_with_results = []
    for bp in branch_points:
        for results_list in results:
            if results_list:
                for result in results_list:
                    if (abs(result['kz'] - bp['kz']) < 1e-6 and
                        abs(result['sz'] - bp['sz']) < 1e-6):
                        points_with_results.append({
                            'branch_point': bp,
                            'result': result
                        })
                        break

    if len(points_with_results) == 0:
        print(f"❌ Не найдено результатов для точек ветви {branch_id}")
        return

    # Сортируем и выбираем точки
    points_with_results.sort(key=lambda p: p['branch_point']['x'])

    if len(points_with_results) > max_points_to_show:
        indices = np.linspace(0, len(points_with_results)-1, max_points_to_show, dtype=int)
        selected_points = [points_with_results[i] for i in indices]
    else:
        selected_points = points_with_results

    # Получаем геометрию из первой точки
    first_result = selected_points[0]['result']
    param_overrides = first_result['params']
    dv = override_params(base_digit_values, param_overrides)

    A_val = None
    h_l_val = h_conductor_l
    h_r_val = h_conductor_r
    for eq in dv:
        if str(eq.lhs()) == 'a':
            A_val = float(eq.rhs())

    if A_val is None:
        A_val = A

    # Создаём массивы x
    n_points = 100
    x_left = np.linspace(-A_val-h_l_val, -A_val, n_points)
    x_vacuum = np.linspace(-A_val, A_val, n_points)
    x_right = np.linspace(A_val, A_val+h_r_val, n_points)

    # Цвета для точек
    colors = plt.cm.viridis(np.linspace(0, 1, len(selected_points)))

    # Подготовка данных
    points_to_plot = []
    for i, pw in enumerate(selected_points):
        bp = pw['branch_point']
        result = pw['result']

        points_to_plot.append({
            'kz': float(bp['kz']),
            'sz': float(bp['sz']),
            'params': result['params'],
            'dv': override_params(base_digit_values, result['params']),
            'label': f"{param_name}={bp['x']:.2e}",
            'color': colors[i],
            'param_value': bp['x'],
            'thrust': bp['thrust'],
            stress_tensor_name: result[stress_tensor_name],
            # Тензорные компоненты (символьные выражения)
        })

    # Вычисляем значения для всех точек
    print(f"Вычисление компонент тензора для {len(points_to_plot)} точек...")
    all_tensor_data = []

    for point in points_to_plot:
        dv_point = point['dv']

        # Вычисляем компоненты тензора
        T_xx_l = eval_tensor_component(point[stress_tensor_name]['T_xx_l'], x_left, dv_point)
        T_xx_v = eval_tensor_component(point[stress_tensor_name]['T_xx_v'], x_vacuum, dv_point)
        T_xx_r = eval_tensor_component(point[stress_tensor_name]['T_xx_r'], x_right, dv_point)

        T_xy_l = eval_tensor_component(point[stress_tensor_name]['T_xy_l'], x_left, dv_point)
        T_xy_v = eval_tensor_component(point[stress_tensor_name]['T_xy_v'], x_vacuum, dv_point)
        T_xy_r = eval_tensor_component(point[stress_tensor_name]['T_xy_r'], x_right, dv_point)

        T_xz_l = eval_tensor_component(point[stress_tensor_name]['T_xz_l'], x_left, dv_point)
        T_xz_v = eval_tensor_component(point[stress_tensor_name]['T_xz_v'], x_vacuum, dv_point)
        T_xz_r = eval_tensor_component(point[stress_tensor_name]['T_xz_r'], x_right, dv_point)

        # Дивергенция
        div_l = eval_tensor_component(point[stress_tensor_name]['div_T_x_l'], x_left, dv_point)
        div_v = eval_tensor_component(point[stress_tensor_name]['div_T_x_v'], x_vacuum, dv_point)
        div_r = eval_tensor_component(point[stress_tensor_name]['div_T_x_r'], x_right, dv_point)

        all_tensor_data.append({
            'x_left': x_left,
            'x_vacuum': x_vacuum,
            'x_right': x_right,
            'T_xx_l': T_xx_l,
            'T_xx_v': T_xx_v,
            'T_xx_r': T_xx_r,
            'T_xy_l': T_xy_l,
            'T_xy_v': T_xy_v,
            'T_xy_r': T_xy_r,
            'T_xz_l': T_xz_l,
            'T_xz_v': T_xz_v,
            'T_xz_r': T_xz_r,
            'div_l': div_l,
            'div_v': div_v,
            'div_r': div_r,
            'color': point['color'],
            'label': point['label'],
            'param_value': point['param_value'],
            'thrust': point['thrust']
        })

    # Создаём графики
    fig = plt.figure(figsize=(18, 14))
    gs = fig.add_gridspec(3, 3, hspace=0.4, wspace=0.3)
    fig.suptitle(f"{lang['tensor_report']} {branch_id}: {param_name}",
                 fontsize=14, fontweight='bold')

    # График 1: T_xx (действительная часть)
    ax1 = fig.add_subplot(gs[0, 0])
    for idx, tensor_data in enumerate(all_tensor_data):
        is_first = (idx == 0)
        alpha = 1.0 if is_first else 0.5
        lw = 2 if is_first else 1

        x_all = np.concatenate([tensor_data['x_left'], tensor_data['x_vacuum'], tensor_data['x_right']])
        # === FIX: используем .real() как метод и конвертируем в float ===
        T_xx_all = np.concatenate([
            [val.n().real() for val in tensor_data['T_xx_l']],
            [val.n().real() for val in tensor_data['T_xx_v']],
            [val.n().real() for val in tensor_data['T_xx_r']]
        ])

        ax1.plot(x_all, T_xx_all, color=tensor_data['color'],
                 alpha=alpha, linewidth=lw, label=tensor_data['label'])
        ax1.axvline(x=-A_val, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
        ax1.axvline(x=A_val, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)

    ax1.set_xlabel(lang['position'], fontsize=12)
    ax1.set_ylabel(lang['stress'], fontsize=12)
    ax1.set_title(f"{lang['T_xx']} - {lang['real_part']}", fontsize=13)
    ax1.grid(True, which='both', linestyle=':', alpha=0.5)
    ax1.legend(fontsize=8, loc='best')

    # График 2: T_xy (действительная часть)
    ax2 = fig.add_subplot(gs[0, 1])
    for idx, tensor_data in enumerate(all_tensor_data):
        is_first = (idx == 0)
        alpha = 1.0 if is_first else 0.5
        lw = 2 if is_first else 1

        x_all = np.concatenate([tensor_data['x_left'], tensor_data['x_vacuum'], tensor_data['x_right']])
        T_xz_all = np.concatenate([
            [val.n().real() for val in tensor_data['T_xy_l']],
            [val.n().real() for val in tensor_data['T_xy_v']],
            [val.n().real() for val in tensor_data['T_xy_r']],
        ])

        ax2.plot(x_all, T_xz_all, color=tensor_data['color'],
                 alpha=alpha, linewidth=lw, label=tensor_data['label'])
        ax2.axvline(x=-A_val, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
        ax2.axvline(x=A_val, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)

    ax2.set_xlabel(lang['position'], fontsize=12)
    ax2.set_ylabel(lang['stress'], fontsize=12)
    ax2.set_title(f"{lang['T_xy']} - {lang['real_part']}", fontsize=13)
    ax2.grid(True, which='both', linestyle=':', alpha=0.5)

    # График 3: T_xz (действительная часть)
    ax3 = fig.add_subplot(gs[0, 2])
    for idx, tensor_data in enumerate(all_tensor_data):
        is_first = (idx == 0)
        alpha = 1.0 if is_first else 0.5
        lw = 2 if is_first else 1

        x_all = np.concatenate([tensor_data['x_left'], tensor_data['x_vacuum'], tensor_data['x_right']])
        T_xz_all = np.concatenate([
            [val.n().real() for val in tensor_data['T_xz_l']],
            [val.n().real() for val in tensor_data['T_xz_v']],
            [val.n().real() for val in tensor_data['T_xz_r']],
        ])

        ax3.plot(x_all, T_xz_all, color=tensor_data['color'],
                 alpha=alpha, linewidth=lw, label=tensor_data['label'])
        ax3.axvline(x=-A_val, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
        ax3.axvline(x=A_val, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)

    ax3.set_xlabel(lang['position'], fontsize=12)
    ax3.set_ylabel(lang['stress'], fontsize=12)
    ax3.set_title(f"{lang['T_xz']} - {lang['real_part']}", fontsize=13)
    ax3.grid(True, which='both', linestyle=':', alpha=0.5)

    # График 4: div T_x (действительная часть)
    ax4 = fig.add_subplot(gs[1, :])
    for idx, tensor_data in enumerate(all_tensor_data):
        is_first = (idx == 0)
        alpha = 1.0 if is_first else 0.5
        lw = 2 if is_first else 1

        x_all = np.concatenate([tensor_data['x_left'], tensor_data['x_vacuum'], tensor_data['x_right']])
        div_all = np.concatenate([
            [val.n().real() for val in tensor_data['div_l']],
            [val.n().real() for val in tensor_data['div_v']],
            [val.n().real() for val in tensor_data['div_r']],
        ])

        ax4.plot(x_all, div_all, color=tensor_data['color'],
                 alpha=alpha, linewidth=lw, label=tensor_data['label'])
        ax4.axvline(x=-A_val, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)
        ax4.axvline(x=A_val, color='gray', linestyle='--', alpha=0.5, linewidth=0.8)

    ax4.set_xlabel(lang['position'], fontsize=12)
    ax4.set_ylabel(lang['force_density'], fontsize=12)
    ax4.set_title(f"{lang['divergence']} - {lang['real_part']}", fontsize=13)
    ax4.grid(True, which='both', linestyle=':', alpha=0.5)

    # График 5: Интегральные силы
    ax5 = fig.add_subplot(gs[2, 0])

    x_params = []
    thrust_values = []
    surface_force_left = []
    surface_force_right = []
    volume_force_total = []

    for tensor_data in all_tensor_data:
        x_params.append(tensor_data['param_value'])
        thrust_values.append(tensor_data['thrust'])

        # Поверхностная сила = скачок T_xx на границе
        T_xx_at_left_boundary_in = tensor_data['T_xx_v'][0].n().real()  # x = -A+0
        T_xx_at_left_boundary_out = tensor_data['T_xx_l'][-1].n().real()  # x = -A-0
        F_surface_left = T_xx_at_left_boundary_in - T_xx_at_left_boundary_out

        T_xx_at_right_boundary_in = tensor_data['T_xx_v'][-1].n().real()  # x = +A-0
        T_xx_at_right_boundary_out = tensor_data['T_xx_r'][0].n().real()  # x = +A+0
        F_surface_right = T_xx_at_right_boundary_out - T_xx_at_right_boundary_in

        # Объёмная сила = интеграл от div T_x
        dx_l = tensor_data['x_left'][1] - tensor_data['x_left'][0]
        dx_v = tensor_data['x_vacuum'][1] - tensor_data['x_vacuum'][0]
        dx_r = tensor_data['x_right'][1] - tensor_data['x_right'][0]

        F_vol_l = np.trapz([val.n().real() for val in tensor_data['div_l']], dx=dx_l)
        F_vol_v = np.trapz([val.n().real() for val in tensor_data['div_v']], dx=dx_v)
        F_vol_r = np.trapz([val.n().real() for val in tensor_data['div_r']], dx=dx_r)

        F_volume = F_vol_l + F_vol_v + F_vol_r

        surface_force_left.append(F_surface_left)
        surface_force_right.append(F_surface_right)
        volume_force_total.append(F_volume)

    ax5.semilogx(x_params, surface_force_left, 'o-', linewidth=2, markersize=8,
                 color='red', alpha=0.7, label=f"{lang['surface_force']} (x=-A)")
    ax5.semilogx(x_params, surface_force_right, 's-', linewidth=2, markersize=8,
                 color='blue', alpha=0.7, label=f"{lang['surface_force']} (x=+A)")
    ax5.semilogx(x_params, volume_force_total, '^-', linewidth=2, markersize=8,
                 color='green', alpha=0.7, label=f"{lang['volume_force']} (∫div T dV)")
    ax5.semilogx(x_params, thrust_values, 'd-', linewidth=2, markersize=8,
                 color='black', alpha=0.9, label=f"{lang['total_force']} (calc_thrust)")
    ax5.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)
    ax5.set_xlabel(param_name, fontsize=12)
    ax5.set_ylabel(lang['force_density'], fontsize=12)
    ax5.set_title(f"{lang['total_force']} by Region", fontsize=13)
    ax5.grid(True, which='both', linestyle=':', alpha=0.5)
    ax5.legend(fontsize=10, loc='best')

    # График 6: Сравнение поверхностной силы и тяги
    ax6 = fig.add_subplot(gs[2, 1])

    total_surface_force = [sl + sr for sl, sr in zip(surface_force_left, surface_force_right)]

    ax6.semilogx(x_params, thrust_values, 'o-', linewidth=2, markersize=8,
                 color='blue', alpha=0.8, label=f"{lang['total_force']} (calc_thrust)")
    ax6.semilogx(x_params, total_surface_force, 's-', linewidth=2, markersize=8,
                 color='red', alpha=0.8, label=f"{lang['surface_force']} (∮T_xx dS)")
    ax6.semilogx(x_params, volume_force_total, '^-', linewidth=2, markersize=8,
                 color='green', alpha=0.8, label=f"{lang['volume_force']} (∫div T dV)")

    ax6.axhline(y=0, color='gray', linestyle='--', linewidth=0.8)
    ax6.set_xlabel(param_name, fontsize=12)
    ax6.set_ylabel(lang['force_density'], fontsize=12)
    ax6.set_title("Force Comparison: Surface vs Volume vs Direct", fontsize=13)
    ax6.grid(True, which='both', linestyle=':', alpha=0.5)
    ax6.legend(fontsize=10, loc='best')

    # plt.show()
    plt.savefig(f"{figdir}/{figfilename}")
    plt.close()

def compute_surface_integral_force(results, param_name, branch_id, base_digit_values,
                                   stress_tensor_name = 'maxwell_stress_tensor',
                                   distance_threshold=0.1,
                                   use_K_params=False,
                                   language='ru'):
    """
    Вычисление поверхностных интегралов тензора Максвелла для трёх областей:
    левый металл, вакуумный зазор, правый феррит.

    Применяет формулу Тамма (105.9) + теорему Гаусса (32.5) отдельно к каждому объёму:

        ∫_V ∂T_ik/∂x_k dV = ∮_S T_xn dS = dP_x/dt + F_x^mech

    Возвращает силы в динах и в Н/кВт для каждой области.

    Parameters:
    -----------
    results : list
        Список результатов из scan_parameter
    param_name : str
        Имя параметра сканирования (обычно 'omega')
    branch_id : int
        ID ветви для анализа
    base_digit_values : list
        Базовые параметры для подстановки
    distance_threshold: параметры сопоставления ветвей
    language : 'ru' или 'en'

    Returns:
    --------
    dict : словарь с результатами для каждой области
    """

    translations = {
        'ru': {
            'surface_integral': f"Поверхностный интеграл тензора {stress_tensor_name}",
            'left_metal': 'Левый металл',
            'vacuum_gap': 'Вакуумный зазор',
            'right_ferrite': 'Правый феррит',
            'branch': 'Ветвь',
            'param': 'Параметр',
            'force_dyn': 'Сила [дин/см²]',
            'force_NkW': 'Удельная сила [Н/кВт]',
            'power_W': 'Мощность [Вт/см²]',
            'momentum_rate': 'dP/dt [дин/см²]',
            'total': 'Сумма',
            'balance_error': 'Погрешность баланса'
        },
        'en': {
            'surface_integral': f"{stress_tensor_name} Stress Tensor Surface Integral",
            'left_metal': 'Left Metal',
            'vacuum_gap': 'Vacuum Gap',
            'right_ferrite': 'Right Ferrite',
            'branch': 'Branch',
            'param': 'Parameter',
            'force_dyn': 'Force [dyn/cm²]',
            'force_NkW': 'Specific Force [N/kW]',
            'power_W': 'Power [W/cm²]',
            'momentum_rate': 'dP/dt [dyn/cm²]',
            'total': 'Total',
            'balance_error': 'Balance Error'
        }
    }

    lang = translations.get(language, translations['ru'])

    # Получаем точки ветви
    matched_branches = match_branches_between_intervals(
        results, param_name, distance_threshold, use_K_params
    )

    if branch_id not in matched_branches:
        print(f"❌ Ветвь {branch_id} не найдена")
        return None

    branch_points = matched_branches[branch_id]
    if len(branch_points) == 0:
        print(f"❌ Ветвь {branch_id} не содержит точек")
        return None

    print(f"\n{'='*80}")
    print(f"{lang['surface_integral']} для ветви {branch_id}")
    print(f"{lang['param']}: {param_name}")
    print(f"{'='*80}\n")

    import numpy as np

    # Результаты для вывода
    output = []

    for bp in branch_points:
        # Ищем соответствующий result
        result = None
        for results_list in results:
            if results_list:
                for r in results_list:
                    if (abs(r['kz'] - bp['kz']) < 1e-6 and
                        abs(r['sz'] - bp['sz']) < 1e-6):
                        result = r
                        break
            if result:
                break

        if not result:
            continue

        # Извлекаем данные
        dv = override_params(base_digit_values, result['params'])
        kz_val = result['kz']
        sz_val = result['sz']
        W_total_W   = result['W_total_W']  # [Вт/см²]
        thrust_N_per_kW = result['thrust_N_per_kW']    # [Н/кВт] из calc_thrust

        # Геометрия
        A_val = None
        for eq in dv:
            if str(eq.lhs()) == 'a':
                A_val = float(eq.rhs())
                break
        if A_val is None:
            A_val = 0.5

        # Эффективная длина по z (для нормировки интеграла по z)
        k_z_Im = abs(float(sz_val))
        L_z_eff = 1.0 / (2.0 * k_z_Im) if k_z_Im > 1e-10 else 1.0  # [см]

        # === ВЫЧИСЛЕНИЕ T_xx НА ГРАНИЦАХ ===
        # T_xx = (1/8π) * (μ·H_x² - ε·E_y² - μ·H_z²) для TE-моды

        def eval_T_xx(T_xx_expr, x_val, dv_local):
            """Вычисляет T_xx в точке"""

            val = T_xx_expr.subs(x == x_val).subs(dv_local)
            # Усреднение по периоду: Re[...]/2 для комплексных амплитуд
            return val.n().real()  # усреднённое значение

        # Граничные значения T_xx (усреднённые по периоду)
        # x = -A: граница металл/вакуум
        T_xx_metal_side    = eval_T_xx (result[stress_tensor_name]['T_xx_l'], -A_val, dv)  # со стороны металла
        T_xx_vac_side_left = eval_T_xx (result[stress_tensor_name]['T_xx_v'], -A_val, dv)  # со стороны вакуума
        print("T_xx_metal_side", T_xx_metal_side)
        print("T_xx_vac_side_left", T_xx_vac_side_left)

        # x = +A: граница вакуум/феррит
        T_xx_vac_side_right = eval_T_xx(result[stress_tensor_name]['T_xx_v'], A_val, dv)  # со стороны вакуума
        T_xx_ferrite_side   = eval_T_xx(result[stress_tensor_name]['T_xx_r'], A_val, dv)  # со стороны феррита
        print("T_xx_vac_side_right", T_xx_vac_side_right)
        print("T_xx_ferrite_side", T_xx_ferrite_side)

        # === ПОВЕРХНОСТНЫЕ ИНТЕГРАЛЫ ДЛЯ КАЖДОЙ ОБЛАСТИ ===
        # Формула: F_x = ∮ T_xn dS = [T_xx]_boundary * L_y * L_z_eff
        # Нормировка на L_y = 1 см, результат на единицу площади по y

        # 1. ЛЕВЫЙ МЕТАЛЛ (объём: x ∈ (-∞, -A])
        # Поверхности: x → -∞ (T=0) и x = -A (нормаль направлена ВЛЕВО, т.е. -x^)
        # F_metal = T_xx|_{x→-∞} - T_xx|_{x=-A^-} = 0 - T_xx_metal_side
        F_metal_dyn = -T_xx_metal_side * L_z_eff  # [дин/см²]
        print("F_metal_dyn", F_metal_dyn)

        # 2. ВАКУУМНЫЙ ЗАЗОР (объём: x ∈ [-A, +A])
        # Поверхности: x = -A (нормаль +x^) и x = +A (нормаль -x^)
        # F_vacuum = T_xx|_{x=-A^+} - T_xx|_{x=+A^-}
        F_vacuum_dyn = (T_xx_vac_side_left - T_xx_vac_side_right) * L_z_eff  # [дин/см²]
        print("F_vacuum_dyn", F_vacuum_dyn)

        # 3. ПРАВЫЙ ФЕРРИТ (объём: x ∈ [+A, +∞))
        # Поверхности: x = +A (нормаль +x^) и x → +∞ (T=0)
        # F_ferrite = T_xx|_{x=+A^+} - T_xx|_{x→+∞} = T_xx_ferrite_side - 0
        F_ferrite_dyn = T_xx_ferrite_side * L_z_eff  # [дин/см²]
        print("L_z_eff", L_z_eff)
        print("F_ferrite_dyn", F_ferrite_dyn)

        # === ПЕРЕВОД В Н/кВт ===
        # 1 Н = 10^5 дин, 1 кВт = 1000 Вт
        # Удельная сила: [Н/кВт] = [дин/см²] * (10^-5 Н/дин) / [Вт/см²] * 1000
        factor_dyn_to_NkW = 1e-5 * 1000 / W_total_W if W_total_W > 1e-20 else 0
        print("W_total_W", W_total_W)
        print("thrust_N_per_kW", thrust_N_per_kW)
        print("factor_dyn_to_NkW", factor_dyn_to_NkW)

        F_metal_NkW   = F_metal_dyn   * factor_dyn_to_NkW
        F_vacuum_NkW  = F_vacuum_dyn  * factor_dyn_to_NkW
        F_ferrite_NkW = F_ferrite_dyn * factor_dyn_to_NkW
        print("F_metal_NkW", F_metal_NkW)
        print("F_vacuum_NkW", F_vacuum_NkW)
        print("F_ferrite_NkW", F_ferrite_NkW)

        # === БАЛАНС ИМПУЛЬСА ===
        # По Тамму §105: сумма поверхностных интегралов = -d/dt ∫g dV
        # Для стационарного режима усреднённое dP/dt = 0, поэтому сумма должна быть ~0
        F_total_dyn = F_metal_dyn + F_vacuum_dyn + F_ferrite_dyn
        F_total_NkW = F_metal_NkW + F_vacuum_NkW + F_ferrite_NkW

        # Погрешность баланса
        balance_error = F_total_NkW - thrust_N_per_kW

        # === СОХРАНЕНИЕ РЕЗУЛЬТАТОВ ===
        entry = {
            'param_value': bp['x'],
            'kz': kz_val,
            'sz': sz_val,
            'L_z_eff': L_z_eff,
            'W_total_W': W_total_W,
            'thrust_N_per_kW': thrust_N_per_kW,
            # Силы в динах
            'F_metal_dyn': F_metal_dyn,
            'F_vacuum_dyn': F_vacuum_dyn,
            'F_ferrite_dyn': F_ferrite_dyn,
            'F_total_dyn': F_total_dyn,
            # Силы в Н/кВт
            'F_metal_NkW': F_metal_NkW,
            'F_vacuum_NkW': F_vacuum_NkW,
            'F_ferrite_NkW': F_ferrite_NkW,
            'F_total_NkW': F_total_NkW,
            # Баланс
            'balance_error': balance_error,
            # Значения T_xx на границах (для отладки)
            'T_xx_metal': T_xx_metal_side,
            'T_xx_vac_left': T_xx_vac_side_left,
            'T_xx_vac_right': T_xx_vac_side_right,
            'T_xx_ferrite': T_xx_ferrite_side,
        }
        output.append(entry)


    return output

def plot_surface_integral_results(surface_results,
                                  param_name,
                                  stress_tensor_name,
                                  figfilename,
                                  language='ru'):
    """
    Визуализация результатов поверхностных интегралов.

    Parameters:
    -----------
    surface_results : list
        Результаты из compute_surface_integral_force
    param_name : str
        Имя параметра для оси X
    language : 'ru' или 'en'
    """

    translations = {
        'ru': {
            'title': f"Поверхностные интегралы тензора {stress_tensor_name}",
            'force_NkW': 'Удельная сила [Н/кВт]',
            'param': 'Параметр',
            'metal': 'Металл',
            'vacuum': 'Вакуум',
            'ferrite': 'Феррит',
            'total': 'Сумма',
            'thrust_N_per_kW': 'thrust_N_per_kW',
            'balance': 'Баланс импульса'
        },
        'en': {
            'title': f"{stress_tensor_name} Stress Tensor Surface Integrals",
            'force_NkW': 'Specific Force [N/kW]',
            'param': 'Parameter',
            'metal': 'Metal',
            'vacuum': 'Vacuum',
            'ferrite': 'Ferrite',
            'total': 'Total',
            'thrust_N_per_kW': 'thrust_N_per_kW',
            'balance': 'Momentum Balance'
        }
    }

    lang = translations.get(language, translations['ru'])

    if not surface_results or len(surface_results) == 0:
        print("❌ Нет данных для построения графиков")
        return

    import matplotlib.pyplot as plt
    import numpy as np

    # Извлекаем данные
    x_vals      = [r['param_value']        for r in surface_results]
    F_metal     = [r['F_metal_NkW']        for r in surface_results]
    F_vacuum    = [r['F_vacuum_NkW']       for r in surface_results]
    F_ferrite   = [r['F_ferrite_NkW']      for r in surface_results]
    F_total     = [r['F_total_NkW']        for r in surface_results]
    thrust_N_per_kW = [r['thrust_N_per_kW'] for r in surface_results]
    balance_err = [r['balance_error']      for r in surface_results]

    # Графики
    fig, axes = plt.subplots(2, 1, figsize=(10, 8))

    # График 1: Силы по областям
    ax1 = axes[0]
    ax1.plot(x_vals, F_metal,     'o-',  label=lang['metal'], linewidth=2)
    ax1.plot(x_vals, F_vacuum,    's-',  label=lang['vacuum'], linewidth=2)
    ax1.plot(x_vals, F_ferrite,   '^-',  label=lang['ferrite'], linewidth=2)
    ax1.plot(x_vals, F_total,     'd--', label=lang['total'], color='black', linewidth=2)
    ax1.plot(x_vals, thrust_N_per_kW, 'x:',  label=lang['thrust_N_per_kW'], color='red', linewidth=2)

    ax1.axhline(y=0, color='gray', linestyle='--', alpha=0.5)
    ax1.set_xlabel(lang['param'])
    ax1.set_ylabel(lang['force_NkW'])
    ax1.set_title(lang['title'])
    ax1.grid(True, which='both', linestyle=':', alpha=0.5)
    ax1.legend(fontsize=9)

    # График 2: Погрешность баланса
    ax2 = axes[1]

    ax2.plot(x_vals, balance_err, 'm-', linewidth=2)
    ax2.set_xlabel(lang['param'])
    ax2.set_ylabel('Относительная погрешность баланса')
    ax2.set_title(lang['balance'])
    ax2.grid(True, which='both', linestyle=':', alpha=0.5)
    ax2.set_yscale('log')

    plt.tight_layout()
    # plt.show()
    plt.savefig(f"{figdir}/{figfilename}")
    plt.close()

# функция для получения лучшей ветви решения:

def get_best_solution_branch(results, param_name,
                             distance_threshold=0.1,
                             use_K_params=False,
                             criterion='max_abs_thrust',
                             thrust_direction='auto',  # 'auto', 'positive', 'negative'
                             verbose=True):
    """
    Находит лучшую ветвь решения по заданному критерию.

    Parameters:
    -----------
    results : list
        Список результатов из scan_parameter
    param_name : str
        Имя параметра сканирования
    distance_threshold : float
        Порог расстояния для сопоставления ветвей
    use_K_params : bool
        Использовать ли K-параметры для сопоставления
    criterion : str
        Критерий выбора лучшей ветви:
        - 'max_thrust': максимальная тяга (с учётом знака)
        - 'min_thrust': минимальная тяга (для отрицательной тяги)
        - 'avg_thrust': средняя тяга
        - 'max_abs_thrust': максимальная абсолютная тяга
        - 'avg_abs_thrust': средняя абсолютная тяга
        - 'integral_thrust': интегральная тяга по параметру
        - 'integral_abs_thrust': интегральная абсолютная тяга
    thrust_direction : str
        Направление тяги:
        - 'auto': автоматически определить по знаку средней тяги
        - 'positive': тяга положительная (вправо)
        - 'negative': тяга отрицательная (влево)
    verbose : bool
        Выводить ли информацию о найденной ветви

    Returns:
    --------
    best_branch_id : int
        ID лучшей ветви
    best_value : float
        Значение критерия для лучшей ветви
    branch_info : dict
        Дополнительная информация о ветви
    sorted_branches : list
        Отсортированный список ветвей
    """

    # Сопоставляем ветви
    matched_branches = match_branches_between_intervals(
        results, param_name, distance_threshold, use_K_params
    )

    if len(matched_branches) == 0:
        print("Нет доступных ветвей")
        return None, None, None, None

    # Автоматически определяем направление тяги
    if thrust_direction == 'auto':
        all_thrust = []
        for branch_id, points in matched_branches.items():
            all_thrust.extend([p['thrust'] for p in points])
        avg_all_thrust = np.mean(all_thrust)
        thrust_direction = 'negative' if avg_all_thrust < 0 else 'positive'
        if verbose:
            print(f"Автоопределение направления тяги: {thrust_direction} (средняя тяга = {avg_all_thrust:.6f})")

    best_branch_id = None
    best_value = -float('inf') if thrust_direction == 'positive' else float('inf')
    branches_stats = {}

    # Анализируем каждую ветвь
    for branch_id, points in matched_branches.items():
        if len(points) == 0:
            continue

        # Извлекаем тягу
        thrust_values = [p['thrust'] for p in points]
        x_values = [p['x'] for p in points]
        k_z_values = [{'param':p['x'], 'k_z' : (p['kz'], p['sz'])} for p in points]

        # Вычисляем критерий
        if criterion == 'max_thrust':
            value = max(thrust_values)
        elif criterion == 'min_thrust':
            value = min(thrust_values)
        elif criterion == 'avg_thrust':
            value = np.mean(thrust_values)
        elif criterion == 'max_abs_thrust':
            value = max(abs(t) for t in thrust_values)
        elif criterion == 'avg_abs_thrust':
            value = np.mean([abs(t) for t in thrust_values])
        elif criterion == 'integral_thrust':
            # Интегрируем по логарифмической шкале
            if len(x_values) > 1:
                log_x = np.log10(x_values)
                value = np.trapz(thrust_values, log_x)
            else:
                value = thrust_values[0] if thrust_values else 0
        elif criterion == 'integral_abs_thrust':
            # Интегрируем абсолютную тягу
            if len(x_values) > 1:
                log_x = np.log10(x_values)
                value = np.trapz([abs(t) for t in thrust_values], log_x)
            else:
                value = abs(thrust_values[0]) if thrust_values else 0
        else:
            raise ValueError(f"Неизвестный критерий: {criterion}")

        # Собираем статистику
        branches_stats[branch_id] = {
            'criterion_value': value,
            'max_thrust': max(thrust_values),
            'min_thrust': min(thrust_values),
            'max_abs_thrust': max(abs(t) for t in thrust_values),
            'avg_thrust': np.mean(thrust_values),
            'avg_abs_thrust': np.mean([abs(t) for t in thrust_values]),
            'std_thrust': np.std(thrust_values),
            'n_points': len(points),
            'x_range': (min(x_values), max(x_values)),
            'k_z_values': k_z_values,
        }

        # Обновляем лучшую ветвь с учётом направления тяги
        if thrust_direction == 'negative':
            # Для отрицательной тяги: лучше = более отрицательная (меньше по значению)
            if criterion in ['min_thrust', 'avg_thrust', 'integral_thrust']:
                if value < best_value:
                    best_value = value
                    best_branch_id = branch_id
            else:  # для abs критериев - больше лучше
                if value > best_value:
                    best_value = value
                    best_branch_id = branch_id
        else:  # positive
            # Для положительной тяги: больше = лучше
            if value > best_value:
                best_value = value
                best_branch_id = branch_id

    if best_branch_id is None:
        print("Не удалось найти лучшую ветвь")
        return None, None, None, None

    best_info = branches_stats[best_branch_id]

    # Сортируем ветви по значению критерия
    if thrust_direction == 'negative' and criterion in ['min_thrust', 'avg_thrust', 'integral_thrust']:
        # Для отрицательной тяги сортируем по возрастанию (более отрицательные первые)
        sorted_branches = sorted(branches_stats.items(),
                                key=lambda x: x[1]['criterion_value'],
                                reverse=False)
    else:
        # Для остальных случаев - по убыванию
        sorted_branches = sorted(branches_stats.items(),
                                key=lambda x: x[1]['criterion_value'],
                                reverse=True)

    if verbose:
        print(f"\n{'='*60}")
        print(f"Лучшая ветвь по критерию '{criterion}': Branch {best_branch_id}")
        print(f"Направление тяги: {thrust_direction}")
        print(f"{'='*60}")
        print(f"  Значение критерия:    {best_value:.6f}")
        print(f"  Максимальная тяга:    {best_info['max_thrust']:.6f} Н/кВт")
        print(f"  Минимальная тяга:     {best_info['min_thrust']:.6f} Н/кВт")
        print(f"  Макс. абс. тяга:      {best_info['max_abs_thrust']:.6f} Н/кВт")
        print(f"  Средняя тяга:         {best_info['avg_thrust']:.6f} Н/кВт")
        print(f"  Средняя абс. тяга:    {best_info['avg_abs_thrust']:.6f} Н/кВт")
        print(f"  Ст. отклонение:       {best_info['std_thrust']:.6f} Н/кВт")
        print(f"  Количество точек:     {best_info['n_points']}")
        print(f"  Диапазон {param_name}: {best_info['x_range'][0]:.2e} — {best_info['x_range'][1]:.2e}")
        print(f"\nСтатистика всех ветвей:")
        print(f"{'Branch':<10} {'Критерий':<12} {'Макс.|тяга|':<12} {'Средн.|тяга|':<12} {'Точек':<8}")
        print("-" * 60)

        for bid, stats in sorted_branches:
            marker = "★" if bid == best_branch_id else " "
            print(f"{marker} {bid:<8} {stats['criterion_value']:<12.6f} "
                  f"{stats['max_abs_thrust']:<12.6f} {stats['avg_abs_thrust']:<12.6f} {stats['n_points']:<8}")
        print(f"{'='*60}\n")

    return best_branch_id, best_value, best_info, sorted_branches

def compare_branches_by_criteria(results, param_name,
                                 distance_threshold=0.1,
                                 use_K_params=False,
                                 thrust_direction='auto'):
    """
    Сравнивает ветви по всем доступным критериям.

    Parameters:
    -----------
    results : list
        Список результатов из scan_parameter
    param_name : str
        Имя параметра сканирования
    distance_threshold : float
        Порог расстояния для сопоставления ветвей
    use_K_params : bool
        Использовать ли K-параметры для сопоставления
    thrust_direction : str
        'auto', 'positive', 'negative'

    Returns:
    --------
    comparison : dict
        Словарь {критерий: (branch_id, value, info)}
    """

    criteria = ['max_abs_thrust', 'avg_abs_thrust', 'integral_abs_thrust',
                'max_thrust', 'min_thrust', 'avg_thrust', 'integral_thrust']

    comparison = {}

    print(f"\n{'='*70}")
    print(f"СРАВНЕНИЕ ВЕТВЕЙ ПО РАЗЛИЧНЫМ КРИТЕРИЯМ")
    print(f"{'='*70}")

    for criterion in criteria:
        branch_id, value, info, _ = get_best_solution_branch(
            results, param_name, distance_threshold,
            use_K_params, criterion, thrust_direction,
            verbose=False
        )

        if branch_id is not None:
            comparison[criterion] = (branch_id, value, info)
            print(f"\n{criterion:20s}: Branch {branch_id}, значение = {value:.6f}")
            print(f"  {'':20s}  макс. |тяга| = {info['max_abs_thrust']:.6f} Н/кВт")
            print(f"  {'':20s}  средн. |тяга| = {info['avg_abs_thrust']:.6f} Н/кВт")
            print(f"  {'':20s}  макс. тяга = {info['max_thrust']:.6f} Н/кВт")
            print(f"  {'':20s}  мин. тяга = {info['min_thrust']:.6f} Н/кВт")

    print(f"\n{'='*70}\n")

    return comparison


def get_solution_branch(results, param_name, branch_id,
                        distance_threshold=0.1,
                        use_K_params=False):

    # Получаем точки ветви
    matched_branches = match_branches_between_intervals(
        results, param_name, distance_threshold, use_K_params
    )

    if branch_id not in matched_branches:
        print(f"❌ Ветвь {branch_id} не найдена")
        return None

    branch_points = matched_branches[branch_id]
    if len(branch_points) == 0:
        print(f"❌ Ветвь {branch_id} не содержит точек")
        return None

    return branch_points

# Перебор всех ветвей
def report_branch_results(branch_id):
    plot_matched_branch_results(results_omega, 'omega',
        figfilename=f"branch{branch_id}_details.png", branch_id=branch_id,
        distance_threshold=0.1,
        use_K_params=False,
        language=language)

    plot_branch_thrust_components(results_omega, 'omega',
        figfilename=f"branch{branch_id}_thrust_components.png",
        branch_id=branch_id, language=language)

    # Вычислить поверхностные интегралы для ветви
    maxwell_surface_results = compute_surface_integral_force(
        results_omega, 'omega', branch_id=branch_id,
        base_digit_values=digit_values_init.copy() + omega_value.copy(),
        stress_tensor_name = 'stress_tensor', language=language)

    # Построить графики поверхностных интегралов
    plot_surface_integral_results(maxwell_surface_results,
        'omega', figfilename=f"branch{branch_id}_stress_tensor.png",
        stress_tensor_name = 'stress_tensor', language=language)

    # Вычислить поверхностные интегралы для ветви
    maxwell_surface_results = compute_surface_integral_force(
        results_omega, 'omega', branch_id=branch_id,
        base_digit_values=digit_values_init.copy() + omega_value.copy(),
        stress_tensor_name = 'maxwell_stress_tensor', language=language)

    # Построить графики поверхностных интегралов
    plot_surface_integral_results(maxwell_surface_results,
        'omega', figfilename=f"branch{branch_id}_maxwell_tensor.png",
        stress_tensor_name = 'maxwell_stress_tensor', language=language)

    # Вычислить поверхностные интегралы для ветви
    convective_surface_results_PE = compute_surface_integral_force(
        results_omega, 'omega', branch_id=branch_id,
        base_digit_values=digit_values_init.copy() + omega_value.copy(),
        stress_tensor_name = 'convective_stress_tensor_PE', language=language)

    # Построить графики поверхностных интегралов
    plot_surface_integral_results(convective_surface_results_PE,
        'omega', figfilename=f"branch{branch_id}_conv_tensor_PE.png",
        stress_tensor_name = 'convective_stress_tensor_PE', language=language)

    # Вычислить поверхностные интегралы для ветви
    convective_surface_results_IH = compute_surface_integral_force(
        results_omega, 'omega', branch_id=branch_id,
        base_digit_values=digit_values_init.copy() + omega_value.copy(),
        stress_tensor_name = 'convective_stress_tensor_IH', language=language)

    # Построить графики поверхностных интегралов
    plot_surface_integral_results(convective_surface_results_IH,
        'omega', figfilename=f"branch{branch_id}_conv_tensor_IH.png",
        stress_tensor_name = 'convective_stress_tensor_IH', language=language)

    # Построить отчёт по полям для конкретной ветви
    plot_field_report_for_branch(results_omega, 'omega',
        base_digit_values=digit_values_init.copy() + omega_value.copy(),
        figfilename=f"branch{branch_id}_fields.png", branch_id=branch_id,
        language=language, max_points_to_show=5)

    # Построить отчёт по тензору для конкретной ветви
    plot_tensor_report_for_branch_2(results_omega, 'omega',
        branch_id=branch_id, base_digit_values=digit_values_init.copy() + omega_value.copy(),
        figfilename=f"branch{branch_id}_stress_tensor_components.png",
        stress_tensor_name = 'stress_tensor', language=language, max_points_to_show=5)

    plot_tensor_report_for_branch(results_omega, 'omega',
        base_digit_values=digit_values_init.copy() + omega_value.copy(),
        figfilename=f"branch{branch_id}_stress_tensor_divergence.png", branch_id=branch_id,
        stress_tensor_name = 'stress_tensor', language=language, max_points_to_show=5)

#     plot_tensor_report_for_branch_2(results_omega, 'omega',
#         base_digit_values=digit_values_init.copy() + omega_value.copy(),
#         figfilename=f"branch{branch_id}_maxwell_tensor_components.png", branch_id=branch_id,
#         stress_tensor_name = 'maxwell_stress_tensor',
#         language=language, max_points_to_show=5)

#     plot_tensor_report_for_branch(results_omega, 'omega',
#         figfilename=f"branch{branch_id}_maxwell_tensor_divergence.png", branch_id=branch_id,
#         base_digit_values=digit_values_init.copy() + omega_value.copy(),
#         stress_tensor_name = 'maxwell_stress_tensor',
#         language=language, max_points_to_show=5)

#     # Построить отчёт по тензору для конкретной ветви
#     plot_tensor_report_for_branch_2(results_omega, 'omega',
#         figfilename=f"branch{branch_id}_conv_tensor_PE_components.png", branch_id=branch_id,
#         base_digit_values=digit_values_init.copy() + omega_value.copy(),
#         stress_tensor_name = 'convective_stress_tensor_PE',
#         language=language, max_points_to_show=5)

#     plot_tensor_report_for_branch(results_omega, 'omega',
#         base_digit_values=digit_values_init.copy() + omega_value.copy(),
#         figfilename=f"branch{branch_id}_conv_tensor_PE_divergence.png", branch_id=branch_id,
#         stress_tensor_name = 'convective_stress_tensor_PE',
#         language=language, max_points_to_show=5)

#     # Построить отчёт по тензору для конкретной ветви
#     plot_tensor_report_for_branch_2(results_omega, 'omega',
#         base_digit_values=digit_values_init.copy() + omega_value.copy(),
#         figfilename=f"branch{branch_id}_conv_tensor_IH_components.png", branch_id=branch_id,
#         stress_tensor_name = 'convective_stress_tensor_IH',
#         language=language, max_points_to_show=5)

#     plot_tensor_report_for_branch(results_omega, 'omega',
#         base_digit_values=digit_values_init.copy() + omega_value.copy(),
#         figfilename=f"branch{branch_id}_conv_tensor_IH_divergence.png", branch_id=branch_id,
#         stress_tensor_name = 'convective_stress_tensor_IH',
#         language=language, max_points_to_show=5)
