

MENDRIVE_LIB_PRECISION='long_double'
# MENDRIVE_LIB_PRECISION='float128'

import struct


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

def unload_lib(lib):
    """Выгрузить загруженную CDLL-библиотеку из памяти."""
    from ctypes.util import find_library
    from ctypes import CDLL, c_void_p, c_int

    if not hasattr(lib, "_handle"):
        raise ValueError("Объект не является загруженной CDLL-библиотекой.")

    # Ищем libdl

    libdl_path = find_library("dl")
    if not libdl_path:
        raise RuntimeError("libdl не найдена — выгрузка .so невозможна.")

    libdl = CDLL(libdl_path)

    dlclose = libdl.dlclose
    dlclose.argtypes = [c_void_p]
    dlclose.restype = c_int

    handle = lib._handle
    ret = dlclose(handle)

    if ret == 0:
        print(f"✅ Библиотека {lib._name} выгружена из памяти.")
        return True
    else:
        # можно вызвать ctypes.get_errno(), но dlclose не устанавливает errno на всех системах
        print(f"❌ dlclose вернул код {ret}. Библиотека, возможно, не выгружена.")
        return False


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
