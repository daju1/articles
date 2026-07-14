from sage.all import RealField, ComplexField
from ctypes import *

# Create a SageMath RealField with high precision
complex128 = ComplexField(128)
real128    = RealField(128)

import numpy as np

N_re = 1000
N_im = 1000

MENDRIVE_LIB_PRECISION='long_double'
# MENDRIVE_LIB_PRECISION='float128'

import struct
from ctypes import Structure, c_longdouble, c_void_p, c_byte, c_ubyte, c_double


def load_lib(name, precision=MENDRIVE_LIB_PRECISION):

    from ctypes import CDLL, POINTER, c_int

    num_type = get_numeric_type(precision)

    # Загрузка библиотеки
    lib = CDLL("./{name}.so".format(name=name))

    lib.det_eval.argtypes = [num_type,
                             num_type,
                             POINTER(num_type),
                             POINTER(num_type)]
    lib.det_eval.restype  = None

    lib.get_sign_K.argtypes = [POINTER(c_int),
                               POINTER(c_int),
                               POINTER(c_int),
                               POINTER(c_int)]
    lib.get_sign_K.restype  = None

    return lib

# Функция вызова:
def det_eval_fast(lib, kz, sz, precision=MENDRIVE_LIB_PRECISION):
    from ctypes import byref
    num_type = get_numeric_type(precision)
    re = num_type()
    im = num_type()
    lib.det_eval(num_type(kz), num_type(sz), byref(re), byref(im))
    return real128(re.value), real128(im.value)

# Функция вызова:
def get_sign_K_fast(lib):
    from ctypes import byref, c_int
    sign_K_H_l_d = c_int()
    sign_K_E_l_d = c_int()
    sign_K_H_r_d = c_int()
    sign_K_E_r_d = c_int()
    lib.get_sign_K(
        byref(sign_K_H_l_d), byref(sign_K_E_l_d),
        byref(sign_K_H_r_d), byref(sign_K_E_r_d))
    return \
        real128(sign_K_H_l_d.value), real128(sign_K_E_l_d.value), \
        real128(sign_K_H_r_d.value), real128(sign_K_E_r_d.value),

# === Глобальное определение Float128 ===
class Float128(Structure):
    """Структура для IEEE 754 binary128 (quad precision)"""
    # Create a byte array type of size 16
    # This ensures no extra padding is added, making the total size exactly 16 bytes
    from ctypes import c_ubyte
    _pack_ = 1
    _fields_ = [
        ("field_0", c_ubyte), # 1 byte
        ("field_1", c_ubyte), # 1 byte
        ("field_2", c_ubyte), # 1 byte
        ("field_3", c_ubyte), # 1 byte
        ("field_4", c_ubyte), # 1 byte
        ("field_5", c_ubyte), # 1 byte
        ("field_6", c_ubyte), # 1 byte
        ("field_7", c_ubyte), # 1 byte
        ("field_8", c_ubyte), # 1 byte
        ("field_9", c_ubyte), # 1 byte
        ("field_a", c_ubyte), # 1 byte
        ("field_b", c_ubyte), # 1 byte
        ("field_c", c_ubyte), # 1 byte
        ("field_d", c_ubyte), # 1 byte
        ("field_e", c_ubyte), # 1 byte
        ("field_f", c_ubyte), # 1 byte
    ]

    def __repr__(self):
        return f"Float128({bytes(self.bytes).hex()})"

def get_numeric_type(precision):
    """
    Возвращает ctypes-тип для заданной точности.

    Args:
        precision: 'long_double' | 'float128' | 'mpfr_512'

    Returns:
        ctypes type для использования в библиотеке
    """
    from ctypes import c_void_p, c_longdouble
    if precision == 'float128':
        return Float128
    elif precision == 'mpfr_512':
        # Для MPFR используем указатель на структуру mpfr_t
        return c_void_p
    else:  # long_double (80-bit x87)
        return c_longdouble


def to_numeric(val, precision):
    """
    Конвертирует значение SageMath/Python в ctypes-тип.

    Args:
        val: числовое значение (real128, float, Decimal, str)
        precision: 'long_double' | 'float128' | 'mpfr_512'

    Returns:
        ctypes-объект готовый для передачи в C
    """
    if precision == 'mpfr_512':
        # MPFR требует инициализации через библиотеку
        # Используем строковое представление для максимальной точности
        raise NotImplementedError(
            "MPFR requires library initialization. "
            "Use mpfr_init2() + mpfr_set_str() in C code."
        )

    elif precision == 'float128':
        return _to_float128(val)

    else:  # long_double
        return _to_long_double(val)


def _to_long_double(val):
    """
    Конвертирует в long double с сохранением точности.
    """
    # Если val — объект SageMath RealField, используем строку
    val_str = str(val)

    # c_longdouble может принимать строку через промежуточный float
    # Для большей точности используем numpy
    try:
        import numpy as np
        return c_longdouble(np.longdouble(val_str))
    except (ImportError, ValueError):
        # Fallback: теряем часть точности
        return c_longdouble(float(val))


def _to_float128(val):
    """
    Конвертирует значение в IEEE 754 binary128 (quad precision).

    Использует libquadmath через ctypes или mpmath для эмуляции.
    """
    f128 = Float128()

    # Метод 1: Через libquadmath (если доступна)
    try:
        from ctypes import CDLL, c_char_p
        libquadmath = CDLL("libquadmath.so.0")

        # strtoflt128 конвертирует строку в __float128
        libquadmath.strtoflt128.argtypes = [c_char_p, c_void_p]
        libquadmath.strtoflt128.restype = Float128

        val_str = str(val).encode('ascii')
        f128 = libquadmath.strtoflt128(val_str, None)
        return f128
    except OSError:
        pass

    # Метод 2: Через mpmath для получения байтов
    try:
        from mpmath import mpf, mp
        mp.prec = 113  # binary128 имеет 113 бит мантиссы

        mpval = mpf(str(val))

        # Получаем компоненты IEEE 754
        sign, man, exp, bc = mpval._mpf_

        # Собираем binary128: 1 бит знак + 15 бит экспонента + 112 бит мантисса
        if man == 0:
            # Ноль
            f128.bytes[:] = [0] * 16
        else:
            # Нормализованное число
            # bias = 16383 для binary128
            biased_exp = exp + bc - 1 + 16383

            if biased_exp <= 0:
                # Денормализованное или underflow
                biased_exp = 0
                man = man >> (1 - (exp + bc - 1 + 16383))
            elif biased_exp >= 0x7FFF:
                # Overflow -> Infinity
                biased_exp = 0x7FFF
                man = 0

            # Убираем скрытый бит
            man_bits = man & ((1 << 112) - 1)

            # Собираем 128 бит: sign(1) | exp(15) | mantissa(112)
            bits = (sign << 127) | (biased_exp << 112) | man_bits

            # Конвертируем в little-endian bytes
            for i in range(16):
                f128.bytes[i] = (bits >> (8 * i)) & 0xFF

        return f128

    except ImportError:
        pass

    # Метод 3: Fallback через double (потеря точности!)
    import warnings
    warnings.warn(
        f"Converting {val} to float128 via double - precision loss!",
        RuntimeWarning
    )

    dval = float(val)
    # Записываем как double в первые 8 байт (это НЕПРАВИЛЬНО для реальных расчётов!)
    dbytes = struct.pack('<d', dval)
    f128.bytes[0:8] = list(dbytes)
    f128.bytes[8:16] = [0] * 8

    return f128


# === Вспомогательные функции для отладки ===

def float128_to_str(f128):
    """Конвертирует Float128 обратно в строку (для проверки)."""
    try:
        from ctypes import CDLL, c_char_p, create_string_buffer
        libquadmath = CDLL("libquadmath.so.0")

        buf = create_string_buffer(50)
        libquadmath.quadmath_snprintf(buf, 50, b"%.36Qe", f128)
        return buf.value.decode('ascii')
    except:
        return f"<Float128: {bytes(f128.bytes).hex()}>"

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
            seg_u.append((float(pt.kz), float(pt.sz)))
        cu.append(seg_u)
        pl += list_plot(seg_u, size=5)

    for i_contour in range(contours.n_im_contours):
        line = contours.im_zero[i_contour]
        seg_v = []
        for i in range(line.n_points):
            pt = line.points[i]
            seg_v.append((float(pt.kz), float(pt.sz)))
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
        _fields_ = [("kz", c_longdouble),
                    ("sz", c_longdouble),
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
        result.append((float(pt.kz), float(pt.sz)))
    return result

def find_resonance_roots_fast(lib, omega_re_range, omega_im_range, n_re=100, n_im=100, make_plot=False, debug_plot=False):
    """
    Быстрый поиск резонансных корней с использованием C-библиотеки.

    Returns:
        Список кортежей (kz, sz) - приближённые положения корней
    """
    import numpy as np

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
    from ctypes import Structure, c_longdouble, POINTER, c_int
    class Point2D(ctypes.Structure):
        _fields_ = [("kz", c_longdouble),
                    ("sz", c_longdouble),
                    ("sin_r_s", c_longdouble)
                   ]

    # Определите структуру
    class CharacteristicRoots(ctypes.Structure):
        _fields_ = [("roots", ctypes.POINTER(Point2D)),
                    ("n_roots", ctypes.c_int)]

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

    omega_graphic_solutions = []

    print(f"Найдено корней: {roots.n_roots}")
    for i in range(roots.n_roots):
        print(f"  re={roots.roots[i].re:.6f}, im={roots.roots[i].im:.6f}, sin_r_s={roots.roots[i].sin_r_s:.6e}")
        omega_graphic_solutions.extend([(roots.roots[i].kz,
                                       roots.roots[i].sz
                                      )])

    cu, cv, pl = plotDetContoursResult(contours)
    pl += list_plot(omega_graphic_solutions, color="green", size=32)
    if make_plot:
        pl.show()
    lib.free_det_contours(contours)
    lib.free_characteristic_roots(roots)
    return omega_graphic_solutions, (cu, cv, pl)

def prec_graphic_solutions(lib, solver, preciser, omega_sol_num, omega_graphic_solutions, precise,
                           verbose=True, return_all_stages=False):
    """
    Уточняет графическое решение методом Ньютона.

    Args:
        return_all_stages: если True, возвращает результаты всех стадий уточнения

    Returns:
        dict с информацией о корне или None при ошибке
    """
    re0, im0 = omega_graphic_solutions[omega_sol_num]
    det0_re, det0_im = det_eval_fast(lib, re0, im0)
    det0_residual = float(np.sqrt((det0_re)**2 + (det0_im)**2))

    if verbose:
        print('omega_sol_index =', omega_sol_num - len(omega_graphic_solutions)/2)
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
            'omega_sol_num': omega_sol_num,
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
            'omega_sol_num': omega_sol_num,
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
    if True:
        omega_graphic_solutions = find_resonance_roots_fast(lib, re_range, im_range)
    else:
        omega_graphic_solutions, (cu, cv, pl) = find_resonance_roots_fast_с(
            lib, re_range, im_range,
            n_re=n_re, n_im=n_im,
            make_plot=make_plot,
            use_tracing=0,
            sharp=sharp,
            eps_sin_r_s = eps_sin_r_s,
            window_size=window_size
        )

    if len(omega_graphic_solutions) == 0:
        print(f"    ⚠️ Корни не найдены для параметров {param_overrides}")
        return None

    # 4. Уточняем все корни и собираем базу
    if verbose:
        print(f"\n  🔍 Уточнение {len(omega_graphic_solutions)} графических решений для параметров {param_overrides}...")

    refined_roots = []
    for omega_sol_num in range(len(omega_graphic_solutions)):
        result = prec_graphic_solutions(
            lib, solver, preciser,
            omega_sol_num, omega_graphic_solutions, precise=False,
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
    for omega_sol_num in range(len(clustered_roots)):
        result = prec_graphic_solutions(
            lib, solver, preciser,
            omega_sol_num, omega_graphic_solutions, precise=True,
            verbose=verbose
        )

        if result is not None:
            refined_clustered_roots.append(result)

    return omega_graphic_solutions, refined_roots, clustered_roots, refined_clustered_roots#, (cu, cv, pl)

class NewtonPrecC:
    def __init__(self, lib, precision):
        self.precision = precision
        self.num_type = get_numeric_type(self.precision)

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

        self.num_type = get_numeric_type(precision)

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
        self.num_type = get_numeric_type(precision)

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
        ret = self.lib.solve_newton_root_classic(self.num_type(re0), self.num_type(im0),
                                    byref(re_sol), byref(im_sol),
                                    byref(f_abs_out), max_iter)

        return real128(re_sol.value), real128(im_sol.value)

class NewtonRootSolver:
    def __init__(self, lib):
        self.lib = lib

        # solve_newton_root
        self.lib.solve_newton_root.argtypes = [
            POINTER(c_double), POINTER(c_double),
            POINTER(c_double), POINTER(c_double),
            c_double, c_int
        ]
        self.lib.solve_newton_root.restype = c_int

    def solve(self, re0, im0, epsabs=1e-12, max_iter=100):
        re = c_double(re0)
        im = c_double(im0)
        omegare_sol = c_double()
        gamma_sol = c_double()
        ret = self.lib.solve_newton_root(byref(re), byref(im),
                                    byref(omegare_sol), byref(gamma_sol),
                                    epsabs, max_iter)

        return real128(omegare_sol.value), real128(gamma_sol.value)

def compute_maxwell_stress_tensor_symbolic_ED(
    Ex, Ey, Ez,
    Dx, Dy, Dz,
    Hx, Hy, Hz,
    Bx, By, Bz,
):
    """
    Символьное вычисление компонент тензора Максвелла и их дивергенции

    Параметры:
    ----------
    Ex, Ey, Ez : символьные выражения компонент E поля
    Dx, Dy, Dz : символьные выражения компонент D поля
    Hx, Hy, Hz : символьные выражения компонент H поля
    Bx, By, Bz : символьные выражения компонент B поля

    Возвращает:
    -----------
    dict : словарь с символьными выражениями:
        - 'T_xx', 'T_xy', 'T_xz' : компоненты тензора
        - 'div_T_x' : дивергенция sum_k dT_xk/dx_k
        - dT_xx_dx = diff(T_xx, x)
        - dT_xy_dy = diff(T_xy, y)
        - dT_xz_dz = diff(T_xz, z)
    """

    # Скалярные произведения E·D и H·B
    ED = Ex*conjugate(Dx) + Ey*conjugate(Dy) + Ez*conjugate(Dz)
    # HB = Hx*conjugate(Bx) + Hy*conjugate(By) + Hz*conjugate(Bz)

    # Компоненты тензора Максвелла (формула 105.1 Тамма)
    # T_{ik} = 1/(8π) { E_i D_k + E_k D_i + H_i B_k + H_k B_i - δ_{ik}(E·D + H·B) }

    # Стандартный тензор Тамма (формула 105.1)
    # T_xx = 1/(8π) * { 2ExDx + 2HxBx - (ED + HB) }

    # T_xx (i=k=x, δ_xx = 1)
    T_xx = (1/(8*pi)) * (
        avg_over_y(2*Ex*conjugate(Dx) - (ED))/2
    )

    # T_xy (i=x, k=y, δ_xy = 0)
    T_xy = (1/(8*pi)) * (
        avg_over_y(Ex*conjugate(Dy) + Ey*conjugate(Dx))/2
    )

    # T_xz (i=x, k=z, δ_xz = 0)
    T_xz = (1/(8*pi)) * (
        avg_over_y(Ex*conjugate(Dz) + Ez*conjugate(Dx))/2
    )

    # Дивергенция: div T_x = ∂T_xx/∂x + ∂T_xy/∂y + ∂T_xz/∂z
    dT_xx_dx = diff(T_xx, x)
    dT_xy_dy = diff(T_xy, y)
    dT_xz_dz = diff(T_xz, z)

    div_T_x = dT_xx_dx + dT_xy_dy + dT_xz_dz

    return {
        'T_xx': T_xx,
        'T_xy': T_xy,
        'T_xz': T_xz,
        'div_T_x': div_T_x,
        'dT_xx_dx': dT_xx_dx,
        'dT_xy_dy': dT_xy_dy,
        'dT_xz_dz': dT_xz_dz,
    }

def compute_maxwell_stress_tensor_symbolic_HB(
    Ex, Ey, Ez,
    Dx, Dy, Dz,
    Hx, Hy, Hz,
    Bx, By, Bz,
):
    """
    Символьное вычисление компонент тензора Максвелла и их дивергенции

    Параметры:
    ----------
    Ex, Ey, Ez : символьные выражения компонент E поля
    Dx, Dy, Dz : символьные выражения компонент D поля
    Hx, Hy, Hz : символьные выражения компонент H поля
    Bx, By, Bz : символьные выражения компонент B поля

    Возвращает:
    -----------
    dict : словарь с символьными выражениями:
        - 'T_xx', 'T_xy', 'T_xz' : компоненты тензора
        - 'div_T_x' : дивергенция sum_k dT_xk/dx_k
        - dT_xx_dx = diff(T_xx, x)
        - dT_xy_dy = diff(T_xy, y)
        - dT_xz_dz = diff(T_xz, z)
    """

    # Скалярные произведения E·D и H·B
    # ED = Ex*conjugate(Dx) + Ey*conjugate(Dy) + Ez*conjugate(Dz)
    HB = Hx*conjugate(Bx) + Hy*conjugate(By) + Hz*conjugate(Bz)

    # Компоненты тензора Максвелла (формула 105.1 Тамма)
    # T_{ik} = 1/(8π) { E_i D_k + E_k D_i + H_i B_k + H_k B_i - δ_{ik}(E·D + H·B) }

    # Стандартный тензор Тамма (формула 105.1)
    # T_xx = 1/(8π) * { 2ExDx + 2HxBx - (ED + HB) }

    # T_xx (i=k=x, δ_xx = 1)
    T_xx = (1/(8*pi)) * (
        avg_over_y(2*Hx*conjugate(Bx) - (HB))/2
    )

    # T_xy (i=x, k=y, δ_xy = 0)
    T_xy = (1/(8*pi)) * (
        avg_over_y(Hx*conjugate(By) + Hy*conjugate(Bx))/2
    )

    # T_xz (i=x, k=z, δ_xz = 0)
    T_xz = (1/(8*pi)) * (
        avg_over_y(Hx*conjugate(Bz) + Hz*conjugate(Bx))/2
    )

    # Дивергенция: div T_x = ∂T_xx/∂x + ∂T_xy/∂y + ∂T_xz/∂z
    dT_xx_dx = diff(T_xx, x)
    dT_xy_dy = diff(T_xy, y)
    dT_xz_dz = diff(T_xz, z)

    div_T_x = dT_xx_dx + dT_xy_dy + dT_xz_dz

    return {
        'T_xx': T_xx,
        'T_xy': T_xy,
        'T_xz': T_xz,
        'div_T_x': div_T_x,
        'dT_xx_dx': dT_xx_dx,
        'dT_xy_dy': dT_xy_dy,
        'dT_xz_dz': dT_xz_dz,
    }

def compute_maxwell_stress_tensor_symbolic(
    Ex, Ey, Ez,
    Dx, Dy, Dz,
    Hx, Hy, Hz,
    Bx, By, Bz,
):
    """
    Символьное вычисление компонент тензора Максвелла и их дивергенции

    Параметры:
    ----------
    Ex, Ey, Ez : символьные выражения компонент E поля
    Dx, Dy, Dz : символьные выражения компонент D поля
    Hx, Hy, Hz : символьные выражения компонент H поля
    Bx, By, Bz : символьные выражения компонент B поля

    Возвращает:
    -----------
    dict : словарь с символьными выражениями:
        - 'T_xx', 'T_xy', 'T_xz' : компоненты тензора
        - 'div_T_x' : дивергенция sum_k dT_xk/dx_k
        - dT_xx_dx = diff(T_xx, x)
        - dT_xy_dy = diff(T_xy, y)
        - dT_xz_dz = diff(T_xz, z)
    """

    # Скалярные произведения E·D и H·B
    ED = Ex*conjugate(Dx) + Ey*conjugate(Dy) + Ez*conjugate(Dz)
    HB = Hx*conjugate(Bx) + Hy*conjugate(By) + Hz*conjugate(Bz)

    # Компоненты тензора Максвелла (формула 105.1 Тамма)
    # T_{ik} = 1/(8π) { E_i D_k + E_k D_i + H_i B_k + H_k B_i - δ_{ik}(E·D + H·B) }

    # Стандартный тензор Тамма (формула 105.1)
    # T_xx = 1/(8π) * { 2ExDx + 2HxBx - (ED + HB) }

    # T_xx (i=k=x, δ_xx = 1)
    T_xx = (1/(8*pi)) * (
        avg_over_y(2*Ex*conjugate(Dx) + 2*Hx*conjugate(Bx) - (ED + HB))/2
    )

    # T_xy (i=x, k=y, δ_xy = 0)
    T_xy = (1/(8*pi)) * (
        avg_over_y(Ex*conjugate(Dy) + Ey*conjugate(Dx) + Hx*conjugate(By) + Hy*conjugate(Bx))/2
    )

    # T_xz (i=x, k=z, δ_xz = 0)
    T_xz = (1/(8*pi)) * (
        avg_over_y(Ex*conjugate(Dz) + Ez*conjugate(Dx) + Hx*conjugate(Bz) + Hz*conjugate(Bx))/2
    )

    # Дивергенция: div T_x = ∂T_xx/∂x + ∂T_xy/∂y + ∂T_xz/∂z
    dT_xx_dx = diff(T_xx, x)
    dT_xy_dy = diff(T_xy, y)
    dT_xz_dz = diff(T_xz, z)

    div_T_x = dT_xx_dx + dT_xy_dy + dT_xz_dz

    return {
        'T_xx': T_xx,
        'T_xy': T_xy,
        'T_xz': T_xz,
        'div_T_x': div_T_x,
        'dT_xx_dx': dT_xx_dx,
        'dT_xy_dy': dT_xy_dy,
        'dT_xz_dz': dT_xz_dz,
    }
def compute_convective_stress_tensor_symbolic_PE(
    Ex, Ey, Ez,
    Dx, Dy, Dz,
    Hx, Hy, Hz,
    Bx, By, Bz,
):
    """
    Символьное вычисление компонент тензора Максвелла и их дивергенции

    Параметры:
    ----------
    Ex, Ey, Ez : символьные выражения компонент E поля
    Dx, Dy, Dz : символьные выражения компонент D поля
    Hx, Hy, Hz : символьные выражения компонент H поля
    Bx, By, Bz : символьные выражения компонент B поля

    Возвращает:
    -----------
    dict : словарь с символьными выражениями:
        - 'T_xx', 'T_xy', 'T_xz' : компоненты тензора
        - 'div_T_x' : дивергенция sum_k dT_xk/dx_k
        - dT_xx_dx = diff(T_xx, x)
        - dT_xy_dy = diff(T_xy, y)
        - dT_xz_dz = diff(T_xz, z)
    """

    # Скалярные произведения E·D и H·B
    EE = avg_over_y(Ex*conjugate(Ex) + Ey*conjugate(Ey) + Ez*conjugate(Ez))/2
    ED = avg_over_y(Ex*conjugate(Dx) + Ey*conjugate(Dy) + Ez*conjugate(Dz))/2
    HH = avg_over_y(Hx*conjugate(Hx) + Hy*conjugate(Hy) + Hz*conjugate(Hz))/2
    HB = avg_over_y(Hx*conjugate(Bx) + Hy*conjugate(By) + Hz*conjugate(Bz))/2

    # T_xx (i=k=x, δ_xx = 1)
    T_xx = (1/(8*pi)) * (
        ED - EE
    )

    # T_xy (i=x, k=y, δ_xy = 0)
    T_xy = (1/(8*pi)) * (
        0
    )

    # T_xz (i=x, k=z, δ_xz = 0)
    T_xz = (1/(8*pi)) * (
        0
    )

    # Дивергенция: div T_x = ∂T_xx/∂x + ∂T_xy/∂y + ∂T_xz/∂z
    dT_xx_dx = diff(T_xx, x)
    dT_xy_dy = diff(T_xy, y)
    dT_xz_dz = diff(T_xz, z)

    div_T_x = dT_xx_dx + dT_xy_dy + dT_xz_dz

    return {
        'T_xx': T_xx,
        'T_xy': T_xy,
        'T_xz': T_xz,
        'div_T_x': div_T_x,
        'dT_xx_dx': dT_xx_dx,
        'dT_xy_dy': dT_xy_dy,
        'dT_xz_dz': dT_xz_dz,
    }
def compute_convective_stress_tensor_symbolic_IH(
    Ex, Ey, Ez,
    Dx, Dy, Dz,
    Hx, Hy, Hz,
    Bx, By, Bz,
):
    """
    Символьное вычисление компонент тензора Максвелла и их дивергенции

    Параметры:
    ----------
    Ex, Ey, Ez : символьные выражения компонент E поля
    Dx, Dy, Dz : символьные выражения компонент D поля
    Hx, Hy, Hz : символьные выражения компонент H поля
    Bx, By, Bz : символьные выражения компонент B поля

    Возвращает:
    -----------
    dict : словарь с символьными выражениями:
        - 'T_xx', 'T_xy', 'T_xz' : компоненты тензора
        - 'div_T_x' : дивергенция sum_k dT_xk/dx_k
        - dT_xx_dx = diff(T_xx, x)
        - dT_xy_dy = diff(T_xy, y)
        - dT_xz_dz = diff(T_xz, z)
    """

    # Скалярные произведения E·D и H·B
    EE = avg_over_y(Ex*conjugate(Ex) + Ey*conjugate(Ey) + Ez*conjugate(Ez))/2
    ED = avg_over_y(Ex*conjugate(Dx) + Ey*conjugate(Dy) + Ez*conjugate(Dz))/2
    HH = avg_over_y(Hx*conjugate(Hx) + Hy*conjugate(Hy) + Hz*conjugate(Hz))/2
    HB = avg_over_y(Hx*conjugate(Bx) + Hy*conjugate(By) + Hz*conjugate(Bz))/2

    # T_xx (i=k=x, δ_xx = 1)
    T_xx = (1/(8*pi)) * (
        HB - HH
    )

    # T_xy (i=x, k=y, δ_xy = 0)
    T_xy = (1/(8*pi)) * (
        0
    )

    # T_xz (i=x, k=z, δ_xz = 0)
    T_xz = (1/(8*pi)) * (
        0
    )

    # Дивергенция: div T_x = ∂T_xx/∂x + ∂T_xy/∂y + ∂T_xz/∂z
    dT_xx_dx = diff(T_xx, x)
    dT_xy_dy = diff(T_xy, y)
    dT_xz_dz = diff(T_xz, z)

    div_T_x = dT_xx_dx + dT_xy_dy + dT_xz_dz

    return {
        'T_xx': T_xx,
        'T_xy': T_xy,
        'T_xz': T_xz,
        'div_T_x': div_T_x,
        'dT_xx_dx': dT_xx_dx,
        'dT_xy_dy': dT_xy_dy,
        'dT_xz_dz': dT_xz_dz,
    }