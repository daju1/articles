

def find_contour_intersections_c(lib, seg_u, seg_v,
                                 eps_det=1e-12,
                                 eps_sin_r_s=1e-3,
                                 sharp=False,
                                 extrap_len = 1e-3,
                                 cos_max_angle = -0.94,
                                 sin_min_angle = 0.3,
                                 sin_max_angle = 0.8,
                                 window_size = 5,
                                 local_angle_staircase_threshold = 0.3,
                                 total_angle_threshold = 0.3,
                                 concentration_threshold = 0.4,
                                 local_angle_sharp_threshold = 0.6,
                                 det_threshold = 1.0
                                ):
    """
    Находит пересечения двух ломаных через C-библиотеку.

    Args:
        seg_u, seg_v: списки [(x,y), ...]
        eps_det: порог для определения параллельности

    Returns:
        Список (x, y) — найденных пересечений
    """
    import numpy as np
    from ctypes import Structure, c_longdouble, POINTER, c_int

    if len(seg_u) < 2 or len(seg_v) < 2:
        return []

    # === Определение структур ===
    class Point2D(ctypes.Structure):
        _fields_ = [("kz", ctypes.c_longdouble),
                    ("sz", ctypes.c_longdouble),
                    ("sin_r_s", c_longdouble)]

    # Преобразуем в массивы long double
    cu_x = (ctypes.c_longdouble * len(seg_u))(*[p[0] for p in seg_u])
    cu_y = (ctypes.c_longdouble * len(seg_u))(*[p[1] for p in seg_u])
    cv_x = (ctypes.c_longdouble * len(seg_v))(*[p[0] for p in seg_v])
    cv_y = (ctypes.c_longdouble * len(seg_v))(*[p[1] for p in seg_v])

    # Максимальное число пересечений
    max_n = 100
    intersections = (Point2D * int(max_n))()

    if sharp:
        # Объявляем
        lib.find_contour_intersections_with_corners.argtypes = [
            POINTER(c_longdouble), POINTER(c_longdouble), c_int,
            POINTER(c_longdouble), POINTER(c_longdouble), c_int,
            POINTER(Point2D), c_int,
            c_longdouble, # eps_det
            c_longdouble, # eps_sin_r_s
            c_longdouble, # extrap_len
            c_longdouble, # cos_max_angle
            c_longdouble, # sin_min_angle
            c_longdouble, # sin_max_angle
            c_int,        # window_size
            c_longdouble, # local_angle_staircase_threshold,  // 0.3L
            c_longdouble, # total_angle_threshold,            // 0.3L
            c_longdouble, # concentration_threshold,          // 0.4L
            c_longdouble, # local_angle_sharp_threshold,      // 0.6L
            c_longdouble, # det_threshold
        ]
        lib.find_contour_intersections_with_corners.restype = c_int

        n_found = lib.find_contour_intersections_with_corners(
            cu_x, cu_y, len(seg_u),
            cv_x, cv_y, len(seg_v),
            intersections, c_int(max_n),
            c_longdouble(eps_det),
            c_longdouble(eps_sin_r_s),
            c_longdouble(extrap_len),
            c_longdouble(cos_max_angle),
            c_longdouble(sin_min_angle),
            c_longdouble(sin_max_angle),
            c_int(window_size),
            c_longdouble(local_angle_staircase_threshold),
            c_longdouble(total_angle_threshold),
            c_longdouble(concentration_threshold),
            c_longdouble(local_angle_sharp_threshold),
            c_longdouble(det_threshold)
        )
    else:
        # Настройка вызова
        lib.find_contour_intersections.argtypes = [
            POINTER(c_longdouble), POINTER(c_longdouble), c_int,
            POINTER(c_longdouble), POINTER(c_longdouble), c_int,
            POINTER(Point2D),
            c_int,
            c_longdouble,
            c_longdouble
        ]
        lib.find_contour_intersections.restype = ctypes.c_int

        n_found = lib.find_contour_intersections(
            cu_x, cu_y, len(seg_u),
            cv_x, cv_y, len(seg_v),
            intersections, c_int(max_n),
            c_longdouble(eps_det),
            c_longdouble(eps_sin_r_s)
        )

    result = []
    for i in range(n_found):
        result.append((intersections[i].kz,
                       intersections[i].sz,
                       intersections[i].sin_r_s))
    return result

def test_sharp_corners(lib, contours, cu, cv,
                       cos_max_angle = -0.94,
                       sin_min_angle = 0.6,
                       sin_max_angle = 0.8,
                       window_size = 5,
                       local_angle_staircase_threshold = 0.3,
                       total_angle_threshold           = 0.3,
                       concentration_threshold         = 0.4,
                       local_angle_sharp_threshold     = 0.6,
                       det_threshold = 1.0,
                       debug_plot = True):

    class Point2D(ctypes.Structure):
        _fields_ = [("kz", c_longdouble),
                    ("sz", c_longdouble),
                    ("sin_r_s", c_longdouble)]

    class Corner2D(ctypes.Structure):
        _fields_ = [("kz", c_longdouble),
                    ("sz", c_longdouble),
                    ("cosine", c_longdouble),
                    ("sine", c_longdouble),
                    ("type", c_int),
                    ("i", c_int),
                    ("refined_i", c_int),
                    ("refined_pt", Point2D),
                   ]

    class ContourLine(ctypes.Structure):
        _fields_ = [("points", POINTER(Point2D)),
                    ("n_points", c_int)]

    # Тестовая функция для острых углов
    lib.test_sharp_corners.argtypes = [ POINTER(ContourLine),
                                       c_longdouble, # cos_max_angle - порог остроты угла, например -0.94L
                                       c_longdouble, # sin_min_angle
                                       c_longdouble, # sin_max_angle
                                       c_int,        #  window_size
                                       c_longdouble,  # local_angle_staircase_threshold,
                                       c_longdouble,  # total_angle_threshold
                                       c_longdouble,  # concentration_threshold
                                       c_longdouble,  # local_angle_sharp_threshold
                                       c_longdouble,  # det_threshold
                                       c_char_p,
                                       POINTER(Corner2D), c_int]

    lib.test_sharp_corners.ret_type = c_int

    # Максимальное число острых углов
    max_sharp_corners = 100
    sharp_corners = (Corner2D * int(max_sharp_corners))()

    cu_sharp_corners = []
    cv_sharp_corners = []

    # Тестируем линии Re=0
    for i in range (0, contours.n_re_contours):
        seg_u = cu[i]
        name = f"Re=0 line={i} points: {len(seg_u)}"
        # Получаем указатель на i-ю линию
        line_ptr = ctypes.cast(
            ctypes.addressof(contours.re_zero.contents) + i * ctypes.sizeof(ContourLine),
            ctypes.POINTER(ContourLine)
        )
        #print("line_ptr.n_points",  line_ptr.n_points)
        n_found = lib.test_sharp_corners(line_ptr,
                               c_longdouble(cos_max_angle),
                               c_longdouble(sin_min_angle),
                               c_longdouble(sin_max_angle),
                               c_int(window_size),
                               c_longdouble(local_angle_staircase_threshold),
                               c_longdouble(total_angle_threshold),
                               c_longdouble(concentration_threshold),
                               c_longdouble(local_angle_sharp_threshold),
                               c_longdouble(det_threshold),
                               c_char_p(name.encode('utf-8')),
                               sharp_corners, c_int(max_sharp_corners));
        print(name, "Найдено острых углов:", n_found)
        corners_u = []
        if debug_plot:  # или условие по длине
            pl_debug = Graphics()
            pl_debug += list_plot(seg_u, size=8)
        for j in range(0, n_found):
            print(sharp_corners[j].kz, sharp_corners[j].sz,
                  sharp_corners[j].cosine,
                  sharp_corners[j].sine,
                  sharp_corners[j].type,
                  sharp_corners[j].i,
                  sharp_corners[j].refined_i,
                  sharp_corners[j].refined_pt.kz,
                  sharp_corners[j].refined_pt.sz,
                 )
            corners_u.append((sharp_corners[j].kz, sharp_corners[j].sz))
            if debug_plot:  # или условие по длине
                if 1 == sharp_corners[j].type:
                    pl_debug += point((sharp_corners[j].kz, sharp_corners[j].sz), size = 36, color = "brown")
                elif 2 == sharp_corners[j].type:
                    pl_debug += point((sharp_corners[j].kz, sharp_corners[j].sz), size = 36, color = "green")
                elif 3 == sharp_corners[j].type:
                    pl_debug += point((sharp_corners[j].kz, sharp_corners[j].sz), size = 36, color = "yellow")
                pl_debug += point((sharp_corners[j].refined_pt.kz, sharp_corners[j].refined_pt.sz), size = 36, color = "pink")
        cu_sharp_corners.append(corners_u)
        if debug_plot:  # или условие по длине
            pl_debug.show()

    # Тестируем линии Im=0
    for i in range (0, contours.n_im_contours):
        seg_v = cv[i]
        name = f"Im=0 line={i} points: {len(seg_u)}"
        line_ptr = ctypes.cast(
            ctypes.addressof(contours.im_zero.contents) + i * ctypes.sizeof(ContourLine),
            ctypes.POINTER(ContourLine)
        )
        #print("line_ptr.n_points",  line_ptr.n_points)
        n_found = lib.test_sharp_corners(line_ptr,
                               c_longdouble(cos_max_angle),
                               c_longdouble(sin_min_angle),
                               c_longdouble(sin_max_angle),
                               c_int(window_size),
                               c_longdouble(local_angle_staircase_threshold),
                               c_longdouble(total_angle_threshold),
                               c_longdouble(concentration_threshold),
                               c_longdouble(local_angle_sharp_threshold),
                               c_longdouble(det_threshold),
                               c_char_p(name.encode('utf-8')),
                               sharp_corners, c_int(max_sharp_corners));
        print(name, "Найдено острых углов:", n_found)
        corners_v = []
        if debug_plot:  # или условие по длине
            pl_debug = Graphics()
            pl_debug += list_plot(seg_v, color="red", size=4)
        for j in range(0, n_found):
            print(sharp_corners[j].kz, sharp_corners[j].sz,
                  sharp_corners[j].cosine,
                  sharp_corners[j].sine,
                  sharp_corners[j].type,
                  sharp_corners[j].i,
                  sharp_corners[j].refined_i,
                  sharp_corners[j].refined_pt.kz,
                  sharp_corners[j].refined_pt.sz,
                 )
            corners_v.append((sharp_corners[j].kz, sharp_corners[j].sz))
            if debug_plot:  # или условие по длине
                if 1 == sharp_corners[j].type:
                    pl_debug += point((sharp_corners[j].kz, sharp_corners[j].sz), size = 36, color = "brown")
                elif 2 == sharp_corners[j].type:
                    pl_debug += point((sharp_corners[j].kz, sharp_corners[j].sz), size = 36, color = "green")
                elif 3 == sharp_corners[j].type:
                    pl_debug += point((sharp_corners[j].kz, sharp_corners[j].sz), size = 36, color = "yellow")
                pl_debug += point((sharp_corners[j].refined_pt.kz, sharp_corners[j].refined_pt.sz), size = 36, color = "pink")
        cv_sharp_corners.append(corners_v)
        if debug_plot:  # или условие по длине
            pl_debug.show()
    return cu_sharp_corners, cv_sharp_corners

def characteristic_plot_c(lib, kz_range, sz_range,
                          nk=100, ns=100,
                          use_traced=False,
                          make_plot=True,
                          debug_plot=False,
                          sharp=False,
                          cos_max_angle = -0.94,
                          sin_min_angle = 0.6,
                          sin_max_angle = 0.8,
                          window_size = 5,
                          local_angle_staircase_threshold = 0.3,
                          total_angle_threshold           = 0.3,
                          concentration_threshold         = 0.4,
                          local_angle_sharp_threshold     = 0.6,
                          det_threshold = 1.0,
                          min_isoline_points_count = 0,
                          isoline_merge_segments_epsilon = 1e-10):
    """
    Быстрый поиск резонансных корней с использованием C-библиотеки.

    Returns:
        Список кортежей (kz, sz) - приближённые положения корней
    """

    import numpy as np
    from ctypes import Structure, c_longdouble, POINTER, c_int, byref, c_char_p

    # === Определение структур ===
    class Point2D(ctypes.Structure):
        _fields_ = [("kz", c_longdouble),
                    ("sz", c_longdouble),
                    ("sin_r_s", c_longdouble)]

    class ContourLine(ctypes.Structure):
        _fields_ = [("points", POINTER(Point2D)),
                    ("n_points", c_int)]

    class DetContoursResult(Structure):
        _fields_ = [
            ("re_zero", POINTER(ContourLine)),
            ("n_re_contours", c_int),
            ("im_zero", POINTER(ContourLine)),
            ("n_im_contours", c_int)
        ]

    # === Настройка аргументов ===
    if use_traced:
        lib.compute_det_contours_traced.argtypes = [
            c_longdouble, c_longdouble, c_int,
            c_longdouble, c_longdouble, c_int,
            POINTER(DetContoursResult),
            c_longdouble
        ]
    else:
        lib.compute_det_contours.argtypes = [
            c_longdouble, c_longdouble, c_int,
            c_longdouble, c_longdouble, c_int,
            POINTER(DetContoursResult),
            c_longdouble, # eps_nan,      // порог для NaN/Inf: если |val| > eps_nan считаем invalid
            c_int,        # min_points_count
            c_longdouble, # СБОРКА ЛОМАНЫХ С ДОПУСКОМ eps_merge = 1e-10L; подберите под масштаб ваших данных
        ]

    lib.free_det_contours.argtypes = [
        POINTER(DetContoursResult)
    ]

    # === Вызов C-функции ===
    contours = DetContoursResult()
    if use_traced:
        ret = lib.compute_det_contours_traced(
            c_longdouble(kz_range[0]), c_longdouble(kz_range[1]), c_int(nk),
            c_longdouble(sz_range[0]), c_longdouble(sz_range[1]), c_int(ns),
            byref(contours),
            c_longdouble(1e300)
        )
    else:
        ret = lib.compute_det_contours(
            c_longdouble(kz_range[0]), c_longdouble(kz_range[1]), c_int(nk),
            c_longdouble(sz_range[0]), c_longdouble(sz_range[1]), c_int(ns),
            byref(contours),
            c_longdouble(1e300),
            c_int(min_isoline_points_count),
            c_longdouble(isoline_merge_segments_epsilon),
        )

    if ret != 0:
        print("Ошибка при построении изолиний")
        return []

    print(f"Re=0: {contours.n_re_contours} линий")
    print(f"Im=0: {contours.n_im_contours} линий")


    cu, cv, pl = plotDetContoursResult(contours)
    if sharp:
        # Тестовая функция для острых углов
        cu_sharp_corners, cv_sharp_corners = test_sharp_corners(lib, contours, cu, cv,
            cos_max_angle = cos_max_angle,
            sin_min_angle = sin_min_angle,
            sin_max_angle = sin_max_angle,
            window_size = window_size,
            local_angle_staircase_threshold = local_angle_staircase_threshold,
            total_angle_threshold           = total_angle_threshold,
            concentration_threshold         = concentration_threshold,
            local_angle_sharp_threshold     = local_angle_sharp_threshold,
            det_threshold = det_threshold,
        );

    # === Поиск пересечений ===
    k_z_graphic_solutions_ = []

    for i_u in range(len(cu)):
        seg_u = cu[i_u]
        if len(seg_u) < nk*10:
            continue
        print("len(seg_u)", len(seg_u))
        if sharp:
            corners_u = cu_sharp_corners[i_u]
        for i_v in range(len(cv)):
            seg_v = cv[i_v]
            if len(seg_v) < nk*10:
                continue
            print("len(seg_v)", len(seg_v))
            if sharp:
                corners_v = cv_sharp_corners[i_v]
            # Интерактивная отладка: показываем пару линий
            if debug_plot:  # или условие по длине
                pl_debug = Graphics()
                pl_debug += list_plot(seg_u, size=8)
                if sharp:
                    pl_debug += list_plot(corners_u, size=36)
                pl_debug += list_plot(seg_v, color="red", size=4)
                if sharp:
                    pl_debug += list_plot(corners_v, color="red", size=18)
                pl_debug.show()

            # Находим пересечения
            intersections = find_contour_intersections_c(lib, seg_u, seg_v, sharp=sharp)
            # print(intersections)
            k_z_graphic_solutions_.extend(intersections)

    pl += list_plot(k_z_graphic_solutions_, color="green", size=16)
    if make_plot:
        pl.show()

    lib.free_det_contours(contours)
    if sharp:
        return k_z_graphic_solutions_, cu_sharp_corners, cv_sharp_corners
    return k_z_graphic_solutions_
