def sage_expr_to_c_with_params(expr, precision, assign_to="det_v"):
    import sympy as sp

    # 1) sage → sympy
    sym_expr = expr._sympy_()

    # 2) подготавливаем замену
    symbol_to_cfield = {
        c:             "params.c",
        omega:         "params.omega",
        a:             "params.a",
        #ifdef KY
        b:             "params.b",
        m:             "params.m",
        #endif

        # левый проводник
        epsilon_l_xx:  "params.eps_l_xx",
        epsilon_l_yy:  "params.eps_l_yy",
        epsilon_l_zz:  "params.eps_l_zz",
        mu_l_xx:       "params.mu_l_xx",
        mu_l_yy:       "params.mu_l_yy",
        mu_l_zz:       "params.mu_l_zz",
        mu_l_yz_c:     "params.mu_l_yz",
        mu_l_zy_c:     "params.mu_l_zy",
        sigma_e_l_xx:  "params.sigma_e_l_xx",
        sigma_e_l_yy:  "params.sigma_e_l_yy",
        sigma_e_l_zz:  "params.sigma_e_l_zz",
        sigma_m_l_xx:  "params.sigma_m_l_xx",
        sigma_m_l_yy:  "params.sigma_m_l_yy",
        sigma_m_l_zz:  "params.sigma_m_l_zz",

        # правый проводник
        epsilon_r_xx:  "params.eps_r_xx",
        epsilon_r_yy:  "params.eps_r_yy",
        epsilon_r_zz:  "params.eps_r_zz",
        mu_r_xx:       "params.mu_r_xx",
        mu_r_yy:       "params.mu_r_yy",
        mu_r_zz:       "params.mu_r_zz",
        mu_r_yz_c:     "params.mu_r_yz",
        mu_r_zy_c:     "params.mu_r_zy",
        sigma_e_r_xx:  "params.sigma_e_r_xx",
        sigma_e_r_yy:  "params.sigma_e_r_yy",
        sigma_e_r_zz:  "params.sigma_e_r_zz",
        sigma_m_r_xx:  "params.sigma_m_r_xx",
        sigma_m_r_yy:  "params.sigma_m_r_yy",
        sigma_m_r_zz:  "params.sigma_m_r_zz",

        mu_0:          "params.mu_0",
        epsilon_0:     "params.epsilon_0",
    }

    # 3) заменяем на dummy-символы
    dummy_subs = {}
    dummy_to_cfield = {}
    for sage_sym, cfield in symbol_to_cfield.items():
        # Имя dummy: например, dummy_c, dummy_omega, dummy_epsilon_l_xx → безопасное в C
        dummy_name = 'dummy_' + str(sage_sym)#.replace('_', '__')  # защищаем от коллизий
        dummy_sym = sp.Symbol(dummy_name)
        sage_sym_sp = sp.Symbol(str(sage_sym))
        dummy_subs[sage_sym_sp] = dummy_sym
        dummy_to_cfield[dummy_name] = cfield

    # 4) делаем замену
    sym_expr_replaced = sym_expr.xreplace(dummy_subs)

    # 4. Генерация кода
    if precision == 'mpfr_512':
        # Для MPFR генерируем код с явными операциями над mpfr_t
        c_str = _generate_mpfr_code(sym_expr, assign_to)
    else:
        # Для float128/long double — стандартный C-код
        c_str = sp.ccode(sym_expr, assign_to=f"__tmp_{assign_to}__")

        # 5) ccode
        c_str = sp.ccode(sym_expr_replaced, assign_to=assign_to)

        # 6) заменяем dummy_X → params.X
        for dummy_name, cfield in dummy_to_cfield.items():
            c_str = c_str.replace(dummy_name, cfield)

        # Замена тригонометрических функций
        c_str = c_str.replace('sqrt(', 'MENDRIVE_COMPLEX_SQRT(')
        c_str = c_str.replace('exp(',  'MENDRIVE_COMPLEX_EXP(')
        c_str = c_str.replace('sin(',  'MENDRIVE_COMPLEX_SIN(')
        c_str = c_str.replace('cos(',  'MENDRIVE_COMPLEX_COS(')
        c_str = c_str.replace('pow(',  'MENDRIVE_COMPLEX_POW(')

    return c_str

def create_det_c( M_det, M_det_K,
                  K_E_v, K_H_v,
                  K_E_l, K_E_r,
                  K_H_l, K_H_r,
                  precision = MENDRIVE_LIB_PRECISION
                ):
    c_expr_det   = sage_expr_to_c_with_params(M_det,   precision=precision, assign_to="det_v"     )
    c_expr_det_K = sage_expr_to_c_with_params(M_det_K, precision=precision, assign_to="det_K_v"   )
    c_expr_K_E_v = sage_expr_to_c_with_params(K_E_v,   precision=precision, assign_to="(*K_E_v_v)")
    c_expr_K_H_v = sage_expr_to_c_with_params(K_H_v,   precision=precision, assign_to="(*K_H_v_v)")
    c_expr_K_E_l = sage_expr_to_c_with_params(K_E_l,   precision=precision, assign_to="(*K_E_l_v)")
    c_expr_K_E_r = sage_expr_to_c_with_params(K_E_r,   precision=precision, assign_to="(*K_E_r_v)")
    c_expr_K_H_l = sage_expr_to_c_with_params(K_H_l,   precision=precision, assign_to="(*K_H_l_v)")
    c_expr_K_H_r = sage_expr_to_c_with_params(K_H_r,   precision=precision, assign_to="(*K_H_r_v)")

    # выражения для производной комплексной функции по каждой компоненте аргумента
    M_det_diff_kz = M_det.diff(kz)
    M_det_diff_sz = M_det.diff(sz)

    # выражения для отношения каждой компоненты комплексной функции к ее производной по каждой компоненте аргумента
    M_det_div_diff_kz = (M_det / M_det_diff_kz)
    M_det_div_diff_sz = (M_det / M_det_diff_sz)

    # выражения для производной комплексной функции по каждой компоненте аргумента
    c_expr_det_diff_kz = sage_expr_to_c_with_params( M_det_diff_kz,
                                                     precision=precision,
                                                     assign_to="df_dkz" )
    c_expr_det_diff_sz = sage_expr_to_c_with_params( M_det_diff_sz,
                                                     precision=precision,
                                                     assign_to="df_dsz" )

    # выражения для отношения каждой компоненты комплексной функции к ее производной по каждой компоненте аргумента
    c_expr_det_div_diff_kz = sage_expr_to_c_with_params( M_det_div_diff_kz,
                                                         precision=precision,
                                                         assign_to="det_div_diff_kz" )
    c_expr_det_div_diff_sz = sage_expr_to_c_with_params( M_det_div_diff_sz,
                                                         precision=precision,
                                                         assign_to="det_div_diff_sz" )

    mendrive_det_c="""
        #include <stdio.h>
        #include <math.h>

        #include "mendrive_det.h"
        #include "mendrive_log.h"

        // Статическая копия параметров — инициализируется один раз
        static mendrive_params_t params;

        static int sign_K_H_l = 1;
        static int sign_K_E_l = 1;
        static int sign_K_H_r = 1;
        static int sign_K_E_r = 1;

        void get_sign_K( int * sign_K_H_l_d, int * sign_K_E_l_d,
                         int * sign_K_H_r_d, int * sign_K_E_r_d ) {{
            *sign_K_H_l_d = sign_K_H_l;
            *sign_K_E_l_d = sign_K_E_l;
            *sign_K_H_r_d = sign_K_H_r;
            *sign_K_E_r_d = sign_K_E_r;
        }}

        void minus_sign_K_H_l() {{
            sign_K_H_l = -1;
        }}
        void plus_sign_K_H_l() {{
            sign_K_H_l = +1;
        }}
        void minus_sign_K_H_r() {{
            sign_K_H_r = -1;
        }}
        void plus_sign_K_H_r() {{
            sign_K_H_r = +1;
        }}
        void minus_sign_K_E_l() {{
            sign_K_E_l = -1;
        }}
        void plus_sign_K_E_l() {{
            sign_K_E_l = +1;
        }}
        void minus_sign_K_E_r() {{
            sign_K_E_r = -1;
        }}
        void plus_sign_K_E_r() {{
            sign_K_E_r = +1;
        }}

        void K_E_v_eval(mendrive_scalar_t kz, mendrive_scalar_t sz,
            mendrive_complex_t *K_E_v_v) {{

            {K_E_v}
        }}

        void K_H_v_eval(mendrive_scalar_t kz, mendrive_scalar_t sz,
            mendrive_complex_t *K_H_v_v) {{

            {K_H_v}
        }}

        void K_E_l_eval(mendrive_scalar_t kz, mendrive_scalar_t sz,
            mendrive_complex_t *K_E_l_v) {{

            {K_E_l}

            if (MENDRIVE_COMPLEX_IMAG(*K_E_l_v) < 0) {{
                minus_sign_K_E_l();
            }}
            else {{
                plus_sign_K_E_l();
            }}
        }}

        void K_H_l_eval(mendrive_scalar_t kz, mendrive_scalar_t sz,
            mendrive_complex_t *K_H_l_v) {{

            {K_H_l}

            if (MENDRIVE_COMPLEX_IMAG(*K_H_l_v) < 0) {{
                minus_sign_K_H_l();
            }}
            else {{
                plus_sign_K_H_l();
            }}
        }}

        void K_E_r_eval(mendrive_scalar_t kz, mendrive_scalar_t sz,
            mendrive_complex_t *K_E_r_v) {{

            {K_E_r}

            if (MENDRIVE_COMPLEX_IMAG(*K_E_r_v) < 0) {{
                minus_sign_K_E_r();
            }}
            else {{
                plus_sign_K_E_r();
            }}
        }}

        void K_H_r_eval(mendrive_scalar_t kz, mendrive_scalar_t sz,
            mendrive_complex_t *K_H_r_v) {{

            {K_H_r}

            if (MENDRIVE_COMPLEX_IMAG(*K_H_r_v) < 0) {{
                minus_sign_K_H_r();
            }}
            else {{
                plus_sign_K_H_r();
            }}
        }}

        void det_eval(mendrive_scalar_t kz, mendrive_scalar_t sz,
            mendrive_scalar_t *det_re, mendrive_scalar_t *det_im) {{

            mendrive_complex_t K_E_vacuum;
            K_E_v_eval(kz, sz, &K_E_vacuum);

            mendrive_complex_t K_H_vacuum;
            K_H_v_eval(kz, sz, &K_H_vacuum);

            mendrive_complex_t K_E_left_conductor;
            K_E_l_eval(kz, sz, &K_E_left_conductor);

            mendrive_complex_t K_H_left_conductor;
            K_H_l_eval(kz, sz, &K_H_left_conductor);

            mendrive_complex_t K_E_right_conductor;
            K_E_r_eval(kz, sz, &K_E_right_conductor);

            mendrive_complex_t K_H_right_conductor;
            K_H_r_eval(kz, sz, &K_H_right_conductor);

            mendrive_complex_t {det_K}

            *det_re = MENDRIVE_COMPLEX_REAL(det_K_v);
            *det_im = MENDRIVE_COMPLEX_IMAG(det_K_v);
        }}

        void det_init(const mendrive_params_t* p) {{
            params = *p;
            #ifdef LOGGING
            MPREC_LOG_DEBUG("params.c=" MPREC_LOG_FMT_SCALAR "\\n", params.c);
#ifndef QNM
            MPREC_LOG_DEBUG("params.omega=" MPREC_LOG_FMT_SCALAR "\\n", params.omega);
#else
            MPREC_LOG_DEBUG("params.L_z=" MPREC_LOG_FMT_SCALAR "\\n", params.L_z);
#endif
            MPREC_LOG_DEBUG("params.a=" MPREC_LOG_FMT_SCALAR "\\n", params.a);
#ifdef KY
            MPREC_LOG_DEBUG("params.b=" MPREC_LOG_FMT_SCALAR "\\n", params.b);
            MPREC_LOG_DEBUG("params.m=%d\\n",  params.m);
#endif
            MPREC_LOG_DEBUG("params.mu_l_xx=" MPREC_LOG_FMT_SCALAR "\\n", params.mu_l_xx);
            MPREC_LOG_DEBUG("params.mu_l_yy=" MPREC_LOG_FMT_SCALAR "\\n", params.mu_l_yy);
            MPREC_LOG_DEBUG("params.mu_l_zz=" MPREC_LOG_FMT_SCALAR "\\n", params.mu_l_zz);
            MPREC_LOG_DEBUG("params.mu_r_xx=" MPREC_LOG_FMT_SCALAR "\\n", params.mu_r_xx);
            MPREC_LOG_DEBUG("params.mu_r_yy=" MPREC_LOG_FMT_SCALAR "\\n", params.mu_r_yy);
            MPREC_LOG_DEBUG("params.mu_r_zz=" MPREC_LOG_FMT_SCALAR "\\n", params.mu_r_zz);
            MPREC_LOG_DEBUG("params.mu_l_yz=" MPREC_LOG_FMT_SCALAR "\\n", params.mu_l_yz);
            MPREC_LOG_DEBUG("params.mu_l_zy=" MPREC_LOG_FMT_SCALAR "\\n", params.mu_l_zy);
            MPREC_LOG_DEBUG("params.mu_r_yz=" MPREC_LOG_FMT_SCALAR "\\n", params.mu_r_yz);
            MPREC_LOG_DEBUG("params.mu_r_zy=" MPREC_LOG_FMT_SCALAR "\\n", params.mu_r_zy);
            MPREC_LOG_DEBUG("params.sigma_e_l_xx=" MPREC_LOG_FMT_SCALAR "\\n", params.sigma_e_l_xx);
            MPREC_LOG_DEBUG("params.sigma_e_l_yy=" MPREC_LOG_FMT_SCALAR "\\n", params.sigma_e_l_yy);
            MPREC_LOG_DEBUG("params.sigma_e_l_zz=" MPREC_LOG_FMT_SCALAR "\\n", params.sigma_e_l_zz);
            MPREC_LOG_DEBUG("params.sigma_e_r_xx=" MPREC_LOG_FMT_SCALAR "\\n", params.sigma_e_r_xx);
            MPREC_LOG_DEBUG("params.sigma_e_r_yy=" MPREC_LOG_FMT_SCALAR "\\n", params.sigma_e_r_yy);
            MPREC_LOG_DEBUG("params.sigma_e_r_zz=" MPREC_LOG_FMT_SCALAR "\\n", params.sigma_e_r_zz);
            MPREC_LOG_DEBUG("params.sigma_m_l_xx=" MPREC_LOG_FMT_SCALAR "\\n", params.sigma_m_l_xx);
            MPREC_LOG_DEBUG("params.sigma_m_l_yy=" MPREC_LOG_FMT_SCALAR "\\n", params.sigma_m_l_yy);
            MPREC_LOG_DEBUG("params.sigma_m_l_zz=" MPREC_LOG_FMT_SCALAR "\\n", params.sigma_m_l_zz);
            MPREC_LOG_DEBUG("params.sigma_m_r_xx=" MPREC_LOG_FMT_SCALAR "\\n", params.sigma_m_r_xx);
            MPREC_LOG_DEBUG("params.sigma_m_r_yy=" MPREC_LOG_FMT_SCALAR "\\n", params.sigma_m_r_yy);
            MPREC_LOG_DEBUG("params.sigma_m_r_zz=" MPREC_LOG_FMT_SCALAR "\\n", params.sigma_m_r_zz);
            MPREC_LOG_DEBUG("params.eps_l_xx=" MPREC_LOG_FMT_SCALAR "\\n", params.eps_l_xx);
            MPREC_LOG_DEBUG("params.eps_l_yy=" MPREC_LOG_FMT_SCALAR "\\n", params.eps_l_yy);
            MPREC_LOG_DEBUG("params.eps_l_zz=" MPREC_LOG_FMT_SCALAR "\\n", params.eps_l_zz);
            MPREC_LOG_DEBUG("params.eps_r_xx=" MPREC_LOG_FMT_SCALAR "\\n", params.eps_r_xx);
            MPREC_LOG_DEBUG("params.eps_r_yy=" MPREC_LOG_FMT_SCALAR "\\n", params.eps_r_yy);
            MPREC_LOG_DEBUG("params.eps_r_zz=" MPREC_LOG_FMT_SCALAR "\\n", params.eps_r_zz);
            MPREC_LOG_DEBUG("params.mu_0=" MPREC_LOG_FMT_SCALAR "\\n", params.mu_0);
            MPREC_LOG_DEBUG("params.epsilon_0=" MPREC_LOG_FMT_SCALAR "\\n", params.epsilon_0);
            fflush(stdout);
            #endif
        }}

        void det_eval_old(mendrive_scalar_t kz, mendrive_scalar_t sz,
            mendrive_scalar_t *det_re, mendrive_scalar_t *det_im) {{
            mendrive_complex_t {det}

            *det_re = MENDRIVE_COMPLEX_REAL(det_v);
            *det_im = MENDRIVE_COMPLEX_IMAG(det_v);
        }}

        void det_derivatives(mendrive_scalar_t kz, mendrive_scalar_t sz,
                             mendrive_scalar_t* dfdkz_re, mendrive_scalar_t* dfdkz_im,
                             mendrive_scalar_t* dfdsz_re, mendrive_scalar_t* dfdsz_im) {{

            mendrive_complex_t K_E_left_conductor;
            mendrive_complex_t K_H_left_conductor;
            mendrive_complex_t K_E_right_conductor;
            mendrive_complex_t K_H_right_conductor;
            K_E_l_eval(kz, sz, &K_E_left_conductor);
            K_H_l_eval(kz, sz, &K_H_left_conductor);
            K_E_r_eval(kz, sz, &K_E_right_conductor);
            K_H_r_eval(kz, sz, &K_H_right_conductor);

            // ∂f/∂kz
            mendrive_complex_t {dfdkz}

            // ∂f/∂sz
            mendrive_complex_t {dfdsz}

            *dfdkz_re = MENDRIVE_COMPLEX_REAL(df_dkz);
            *dfdkz_im = MENDRIVE_COMPLEX_IMAG(df_dkz);
            *dfdsz_re = MENDRIVE_COMPLEX_REAL(df_dsz);
            *dfdsz_im = MENDRIVE_COMPLEX_IMAG(df_dsz);
        }}

        void det_div_diff_kz_eval(mendrive_scalar_t kz, mendrive_scalar_t sz,
            mendrive_scalar_t *div_re, mendrive_scalar_t *div_im) {{

            mendrive_complex_t K_E_left_conductor;
            mendrive_complex_t K_H_left_conductor;
            mendrive_complex_t K_E_right_conductor;
            mendrive_complex_t K_H_right_conductor;
            K_E_l_eval(kz, sz, &K_E_left_conductor);
            K_H_l_eval(kz, sz, &K_H_left_conductor);
            K_E_r_eval(kz, sz, &K_E_right_conductor);
            K_H_r_eval(kz, sz, &K_H_right_conductor);

            mendrive_complex_t {divkz}

            *div_re = MENDRIVE_COMPLEX_REAL(det_div_diff_kz);
            *div_im = MENDRIVE_COMPLEX_IMAG(det_div_diff_kz);
        }}

        void det_div_diff_sz_eval(mendrive_scalar_t kz, mendrive_scalar_t sz,
            mendrive_scalar_t *div_re, mendrive_scalar_t *div_im) {{

            mendrive_complex_t K_E_left_conductor;
            mendrive_complex_t K_H_left_conductor;
            mendrive_complex_t K_E_right_conductor;
            mendrive_complex_t K_H_right_conductor;
            K_E_l_eval(kz, sz, &K_E_left_conductor);
            K_H_l_eval(kz, sz, &K_H_left_conductor);
            K_E_r_eval(kz, sz, &K_E_right_conductor);
            K_H_r_eval(kz, sz, &K_H_right_conductor);

            mendrive_complex_t {divsz}

            *div_re = MENDRIVE_COMPLEX_REAL(det_div_diff_sz);
            *div_im = MENDRIVE_COMPLEX_IMAG(det_div_diff_sz);
        }}

        """.format(det = c_expr_det,
                   det_K = c_expr_det_K,
                   K_E_v = c_expr_K_E_v, K_H_v = c_expr_K_H_v,
                   K_E_l = c_expr_K_E_l, K_E_r = c_expr_K_E_r,
                   K_H_l = c_expr_K_H_l, K_H_r = c_expr_K_H_r,
                   dfdkz = c_expr_det_diff_kz,
                   dfdsz = c_expr_det_diff_sz,
                   divkz = c_expr_det_div_diff_kz,
                   divsz = c_expr_det_div_diff_sz)

    # print(mendrive_det_c)

    with open("mendrive_det_real.c", "w") as file:
        file.write(mendrive_det_c)
import ctypes
from ctypes.util import find_library

def unload_lib(lib):
    """Выгрузить загруженную CDLL-библиотеку из памяти."""
    if not hasattr(lib, "_handle"):
        raise ValueError("Объект не является загруженной CDLL-библиотекой.")

    # Ищем libdl
    libdl_path = find_library("dl")
    if not libdl_path:
        raise RuntimeError("libdl не найдена — выгрузка .so невозможна.")

    libdl = ctypes.CDLL(libdl_path)

    dlclose = libdl.dlclose
    dlclose.argtypes = [ctypes.c_void_p]
    dlclose.restype = ctypes.c_int

    handle = lib._handle
    ret = dlclose(handle)

    if ret == 0:
        print(f"✅ Библиотека {lib._name} выгружена из памяти.")
        return True
    else:
        # можно вызвать ctypes.get_errno(), но dlclose не устанавливает errno на всех системах
        print(f"❌ dlclose вернул код {ret}. Библиотека, возможно, не выгружена.")
        return False
def compile_lib(name, precision=MENDRIVE_LIB_PRECISION, loglevel='DEBUG'):
    # Command to execute -lgsl -lgslcblas -lm
    # cmd = "gcc -shared -fPIC -O3 -o {name}.so {name}.c -lm".format(name=name)
    """
    Компилирует библиотеку с заданной арифметикой.

    Args:
        name: имя библиотеки (без расширения)
        precision: 'long_double' | 'float128' | 'mpfr_512'
    """

    # Базовые файлы проекта
    source_files = [
        "mendrive_det_real.c",
        "mendrive_isolines.c",
        "mendrive_isolines_traced.c",
        "mendrive_contour_intersections.c",
        "mendrive_contour_intersections_sharp.c",
        "mendrive_characteristic_roots.c",
        "mendrive_root.c",
        "mendrive_newton.c",
        "mendrive_newton_classic/mendrive_newton_classic.c",
        "mendrive_linsolve.c"
    ]

    # Выбор флагов компиляции
    if precision == 'float128':
        define_flag = "-DARITH_FLOAT128"
        extra_libs = "-lquadmath"
    elif precision == 'mpfr_512':
        define_flag = "-DARITH_MPFR_512"
        extra_libs = "-lmpfr -lgmp"
    else:  # long_double
        define_flag = "-DARITH_LONG_DOUBLE"
        extra_libs = "-lm"

    if 'DEBUG' == loglevel:
        define_flag += ' -DMPREC_LOG_LEVEL=MPREC_LOG_LEVEL_DEBUG'
    elif 'TRACE' == loglevel:
        define_flag += ' -DMPREC_LOG_LEVEL=MPREC_LOG_LEVEL_TRACE'
    else:
        define_flag += ' -DMPREC_LOG_LEVEL=MPREC_LOG_LEVEL_OFF'

    cmd = f"gcc -shared -fPIC -O3 {define_flag} -I../vimaniks/gsl/local/include " \
          f"-L../vimaniks/gsl/local/lib -Wl,-rpath='$ORIGIN/../vimaniks/gsl/local/lib' " \
          f"{' '.join(source_files)} -o {name}.so {extra_libs} " \
          f"-DKY -Wl,-rpath='$ORIGIN/../vimaniks/gsl/local/lib'-lgsl -lgslcblas -lm"

    cmd = f"""gcc -shared -fPIC -O3 -o {name}.so \
        mendrive_det_real.c \
        mendrive_isolines.c \
        mendrive_isolines_traced.c \
        mendrive_contour_intersections.c \
        mendrive_contour_intersections_sharp.c \
        mendrive_characteristic_roots.c \
        mendrive_root.c \
        mendrive_newton.c \
        mendrive_newton_classic/mendrive_newton_classic.c \
        mendrive_linsolve.c \
        -I../vimaniks/gsl/local/include \
        -L../vimaniks/gsl/local/lib \
        {define_flag} \
        -DKY \
        -Wl,-rpath='$ORIGIN/../vimaniks/gsl/local/lib' \
        -lgsl -lgslcblas -lm {extra_libs}
    """#.format(name=name, precision=precision)

    import os
    # Using os.system() method


    print(f"🔧 Компиляция ({precision}):")
    print(f"   {cmd}")
    ret = os.system(cmd)

    if ret == 0:
        print(f"✅ Библиотека {name}.so скомпилирована ({precision})")
    else:
        raise RuntimeError(f"❌ Ошибка компиляции (код {ret})")

    return ret == 0
import struct
from ctypes import Structure, c_longdouble, c_void_p, c_byte, c_ubyte, c_double

# === Глобальное определение Float128 ===
class Float128(Structure):
    """Структура для IEEE 754 binary128 (quad precision)"""
    # Create a byte array type of size 16
    # This ensures no extra padding is added, making the total size exactly 16 bytes
    _pack_ = 1
    _fields_ = [
        ("field_0", ctypes.c_ubyte), # 1 byte
        ("field_1", ctypes.c_ubyte), # 1 byte
        ("field_2", ctypes.c_ubyte), # 1 byte
        ("field_3", ctypes.c_ubyte), # 1 byte
        ("field_4", ctypes.c_ubyte), # 1 byte
        ("field_5", ctypes.c_ubyte), # 1 byte
        ("field_6", ctypes.c_ubyte), # 1 byte
        ("field_7", ctypes.c_ubyte), # 1 byte
        ("field_8", ctypes.c_ubyte), # 1 byte
        ("field_9", ctypes.c_ubyte), # 1 byte
        ("field_a", ctypes.c_ubyte), # 1 byte
        ("field_b", ctypes.c_ubyte), # 1 byte
        ("field_c", ctypes.c_ubyte), # 1 byte
        ("field_d", ctypes.c_ubyte), # 1 byte
        ("field_e", ctypes.c_ubyte), # 1 byte
        ("field_f", ctypes.c_ubyte), # 1 byte
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

def load_lib(name, precision=MENDRIVE_LIB_PRECISION):
    import ctypes
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

def init_lib(lib, digit_values, precision=MENDRIVE_LIB_PRECISION):

    #from ctypes import Structure, c_longdouble, POINTER, CDLL, byref, c_int

    num_type = get_numeric_type(precision)

    class MendriveParams(Structure):
        _fields_ = [
            ("c",     num_type),
            ("omega", num_type),
            ("a",     num_type),
            ("b",     num_type),
            ("m",     c_int),

            ("eps_l_xx",     num_type), ("eps_l_yy",     num_type), ("eps_l_zz",     num_type),
            ("mu_l_xx",      num_type), ("mu_l_yy",      num_type), ("mu_l_zz",      num_type),
            ("mu_l_yz",      num_type), ("mu_l_zy",      num_type),
            ("sigma_e_l_xx", num_type), ("sigma_e_l_yy", num_type), ("sigma_e_l_zz", num_type),
            ("sigma_m_l_xx", num_type), ("sigma_m_l_yy", num_type), ("sigma_m_l_zz", num_type),

            ("eps_r_xx",     num_type), ("eps_r_yy",     num_type), ("eps_r_zz",     num_type),
            ("mu_r_xx",      num_type), ("mu_r_yy",      num_type), ("mu_r_zz",      num_type),
            ("mu_r_yz",      num_type), ("mu_r_zy",      num_type),
            ("sigma_e_r_xx", num_type), ("sigma_e_r_yy", num_type), ("sigma_e_r_zz", num_type),
            ("sigma_m_r_xx", num_type), ("sigma_m_r_yy", num_type), ("sigma_m_r_zz", num_type),

            ("mu_0",         num_type),
            ("epsilon_0",    num_type),
        ]

    # Привязка функций
    lib.det_init.argtypes = [POINTER(MendriveParams)]
    lib.det_init.restype  = None

    # === Формируем параметры из digit_values ===
    p = MendriveParams()
    for k_v in digit_values:
        name = str(k_v.lhs())
        v = to_numeric(real128(k_v.rhs()), precision=precision)
        if name == 'c':        p.c = v
        elif name == 'omega':  p.omega = v
        elif name == 'a':      p.a = v
        elif name == 'b':      p.b = v
        elif name == 'm':      p.m = int(k_v.rhs())
        # left conductor
        elif name == 'epsilon_l_xx': p.eps_l_xx = v
        elif name == 'epsilon_l_yy': p.eps_l_yy = v
        elif name == 'epsilon_l_zz': p.eps_l_zz = v
        elif name == 'mu_l_xx':      p.mu_l_xx = v
        elif name == 'mu_l_yy':      p.mu_l_yy = v
        elif name == 'mu_l_zz':      p.mu_l_zz = v
        elif name == 'mu_l_yz_c':    p.mu_l_yz = v
        elif name == 'mu_l_zy_c':    p.mu_l_zy = v
        elif name == 'sigma_e_l_xx': p.sigma_e_l_xx = v
        elif name == 'sigma_e_l_yy': p.sigma_e_l_yy = v
        elif name == 'sigma_e_l_zz': p.sigma_e_l_zz = v
        elif name == 'sigma_m_l_xx': p.sigma_m_l_xx = v
        elif name == 'sigma_m_l_yy': p.sigma_m_l_yy = v
        elif name == 'sigma_m_l_zz': p.sigma_m_l_zz = v
        # right conductor
        elif name == 'epsilon_r_xx': p.eps_r_xx = v
        elif name == 'epsilon_r_yy': p.eps_r_yy = v
        elif name == 'epsilon_r_zz': p.eps_r_zz = v
        elif name == 'mu_r_xx':      p.mu_r_xx = v
        elif name == 'mu_r_yy':      p.mu_r_yy = v
        elif name == 'mu_r_zz':      p.mu_r_zz = v
        elif name == 'mu_r_yz_c':    p.mu_r_yz = v
        elif name == 'mu_r_zy_c':    p.mu_r_zy = v
        elif name == 'sigma_e_r_xx': p.sigma_e_r_xx = v
        elif name == 'sigma_e_r_yy': p.sigma_e_r_yy = v
        elif name == 'sigma_e_r_zz': p.sigma_e_r_zz = v
        elif name == 'sigma_m_r_xx': p.sigma_m_r_xx = v
        elif name == 'sigma_m_r_yy': p.sigma_m_r_yy = v
        elif name == 'sigma_m_r_zz': p.sigma_m_r_zz = v
        elif name == 'mu_0': p.mu_0 = v
        elif name == 'epsilon_0': p.epsilon_0 = v

    # Инициализация один раз:
    lib.det_init(byref(p))


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
