#ifndef MENDDRIVE_DET_H
#define MENDDRIVE_DET_H

// ============================================================================
// СИСТЕМА ВЫБОРА АРИФМЕТИКИ
// ============================================================================
// Выберите ОДИН режим перед компиляцией:
//   -DARITH_LONG_DOUBLE    : стандартный long double (80-bit, ~19 цифр)
//   -DARITH_FLOAT128       : __float128 (128-bit, ~34 цифры) [GCC/Clang]
//   -DARITH_MPFR_512       : MPFR с 512 битами (~154 цифры)
// ============================================================================

// Тип данных для скалярных величин
#ifdef ARITH_FLOAT128
    #include <complex.h>
    #include <quadmath.h>
    typedef __float128 mendrive_scalar_t;
    typedef __complex128 mendrive_complex_t;

    #define MENDRIVE_SET_D(z, val) (z) = (__float128)(val)
    #define MENDRIVE_COMPLEX_SET_D(z, re_val, im_val) (z) = (__float128)(re_val) + (__float128)(im_val) * (__complex128)1.0qi
    #define MENDRIVE_ABS(result, x) (result) = fabsq(x)
    #define MENDRIVE_COMPLEX_ABS(result, z) (result) = cabsq(z)
    #define MENDRIVE_COMPLEX_REAL crealq
    #define MENDRIVE_COMPLEX_IMAG cimagq

    #define MENDRIVE_SQRT sqrtq
    #define MENDRIVE_FABS fabsq
    #define MENDRIVE_POW powq
    #define MENDRIVE_SIN sinq
    #define MENDRIVE_COS cosq
    #define MENDRIVE_EXP expq
    #define MENDRIVE_LOG logq

    #define MENDRIVE_COMPLEX_SQRT csqrtq
    #define MENDRIVE_COMPLEX_FABS cfabsq
    #define MENDRIVE_COMPLEX_POW cpowq
    #define MENDRIVE_COMPLEX_SIN csinq
    #define MENDRIVE_COMPLEX_COS ccosq
    #define MENDRIVE_COMPLEX_EXP cexpq
    #define MENDRIVE_COMPLEX_LOG clogq

    #define MENDRIVE_SNPRINTF(buf, size, fmt, x) snprintfq(buf, size, fmt, x)

#elif defined(ARITH_LONG_DOUBLE)
    #include <complex.h>
    typedef long double mendrive_scalar_t;
    // typedef long double complex mendrive_complex_t;
    typedef long double _Complex mendrive_complex_t;

    #define MENDRIVE_SET_D(z, val) (z) = (long double)(val)
    #define MENDRIVE_COMPLEX_SET_D(z, re_val, im_val) (z) = (long double)(re_val) + (long double)(im_val) * I
    #define MENDRIVE_ABS(result, x) (result) = fabsl(x)
    #define MENDRIVE_COMPLEX_ABS(result, z) (result) = cabsl(z)

    #define MENDRIVE_COMPLEX_REAL creall
    #define MENDRIVE_COMPLEX_IMAG cimagl

    #define MENDRIVE_SQRT sqrtl
    #define MENDRIVE_FABS fabsl
    #define MENDRIVE_POW powl
    #define MENDRIVE_SIN sinl
    #define MENDRIVE_COS cosl
    #define MENDRIVE_EXP expl
    #define MENDRIVE_LOG logl

    #define MENDRIVE_COMPLEX_SQRT csqrtl
    #define MENDRIVE_COMPLEX_FABS cfabsl
    #define MENDRIVE_COMPLEX_POW cpowl
    #define MENDRIVE_COMPLEX_SIN csinl
    #define MENDRIVE_COMPLEX_COS ccosl
    #define MENDRIVE_COMPLEX_EXP cexpl
    #define MENDRIVE_COMPLEX_LOG clogl

    #define MENDRIVE_SNPRINTF(buf, size, fmt, x) snprintf(buf, size, fmt, (long double)(x))

#elif defined(ARITH_MPFR_512)

    #include <mpfr.h>
    #include <gmp.h>
    #define MPFR_PREC 512
    typedef mpfr_t mendrive_scalar_t;  // mpfr_t — это массив, передаётся по ссылке
    // Для MPFR все операции — через функции, макросы ниже не используются напрямую

    // Для комплексных чисел используем пару mpfr_t
    typedef struct {
        mpfr_t re;
        mpfr_t im;
    } mprec_complex_t;

    // Инициализация комплексного числа
    #define MENDRIVE_COMPLEX_INIT(z) \
        mpfr_init2((z).re, MENDRIVE_BITS); \
        mpfr_init2((z).im, MENDRIVE_BITS)

    // Очистка
    #define MENDRIVE_COMPLEX_CLEAR(z) \
        mpfr_clear((z).re); \
        mpfr_clear((z).im)

    // Установка из double
    #define MENDRIVE_SET_D(z, val) mpfr_set_d(z, val, MPFR_RNDN)
    #define MENDRIVE_COMPLEX_SET_D(z, re_val, im_val) \
        MENDRIVE_SET_D((z).re, re_val); \
        MENDRIVE_SET_D((z).im, im_val)

    // Абсолютное значение
    #define MENDRIVE_ABS(result, x) mpfr_abs(result, x, MPFR_RNDN)
    #define MENDRIVE_COMPLEX_ABS(result, z) \
        do { \
            mpfr_t tmp1, tmp2; \
            mpfr_init2(tmp1, MENDRIVE_BITS); \
            mpfr_init2(tmp2, MENDRIVE_BITS); \
            mpfr_sqr(tmp1, (z).re, MPFR_RNDN); \
            mpfr_sqr(tmp2, (z).im, MPFR_RNDN); \
            mpfr_add(result, tmp1, tmp2, MPFR_RNDN); \
            mpfr_sqrt(result, result, MPFR_RNDN); \
            mpfr_clear(tmp1); mpfr_clear(tmp2); \
        } while(0)

#endif

// ============================================================================
// Структура параметров (единая для всех режимов точности)
// ============================================================================
typedef struct {
    mendrive_scalar_t c;
#ifndef QNM
    mendrive_scalar_t omega;
#else
    mendrive_scalar_t L_z;
#endif
    mendrive_scalar_t a;
#ifdef KY
    mendrive_scalar_t b;
    int m;
#endif

    // Левый проводник
    mendrive_scalar_t eps_l_xx, eps_l_yy, eps_l_zz;
#ifndef GYRO_TENSOR
    mendrive_scalar_t mu_l_xx, mu_l_yy, mu_l_zz;
    mendrive_scalar_t mu_l_yz, mu_l_zy;
#else
    mendrive_scalar_t mu_l_perp;
    mendrive_scalar_t kappa_l;
    mendrive_scalar_t mu_l_parallel;
#endif
    mendrive_scalar_t sigma_e_l_xx, sigma_e_l_yy, sigma_e_l_zz;
#ifndef GYRO_TENSOR
    mendrive_scalar_t sigma_m_l_xx, sigma_m_l_yy, sigma_m_l_zz;
#else
    mendrive_scalar_t sigma_m_l_perp;
    mendrive_scalar_t sigma_m_l_gyro;
#endif

    // Правый проводник
    mendrive_scalar_t eps_r_xx, eps_r_yy, eps_r_zz;
#ifndef GYRO_TENSOR
    mendrive_scalar_t mu_r_xx, mu_r_yy, mu_r_zz;
    mendrive_scalar_t mu_r_yz, mu_r_zy;
#else
    mendrive_scalar_t mu_r_perp;
    mendrive_scalar_t kappa_r;
    mendrive_scalar_t mu_r_parallel;
#endif
    mendrive_scalar_t sigma_e_r_xx, sigma_e_r_yy, sigma_e_r_zz;
#ifndef GYRO_TENSOR
    mendrive_scalar_t sigma_m_r_xx, sigma_m_r_yy, sigma_m_r_zz;
#else
    mendrive_scalar_t sigma_m_r_perp;
    mendrive_scalar_t sigma_m_r_gyro;
#endif
    mendrive_scalar_t mu_0;
    mendrive_scalar_t epsilon_0;
} mendrive_params_t;

// ============================================================================
// ГЛОБАЛЬНЫЕ ПЕРЕМЕННЫЕ И ФУНКЦИИ
// ============================================================================
#ifdef ARITH_MPFR_512

    void det_init(const mendrive_params_t* p);
    void det_clear(void);  // Очистка глобальных переменных MPFR

    void det_eval(
        const mpfr_t kz, const mpfr_t sz,
        mpfr_t det_re, mpfr_t det_im
    );

    void det_div_diff_kz_eval(
        const mpfr_t kz, const mpfr_t sz,
        mpfr_t div_re, mpfr_t div_im
    );

    void det_div_diff_sz_eval(
        const mpfr_t kz, const mpfr_t sz,
        mpfr_t div_re, mpfr_t div_im
    );

    void det_derivatives(const mpfr_t kz, const mpfr_t sz,
                     mpfr_t dfdkz_re_out, mpfr_t dfdkz_im_out,
                     mpfr_t dfdsz_re_out, mpfr_t dfdsz_im_out);

#else
    // LONG DOUBLE и FLOAT128 используют одинаковый интерфейс

    void det_init(const mendrive_params_t* p);

    void det_eval_old(
        mendrive_scalar_t kz, mendrive_scalar_t sz,
        mendrive_scalar_t* det_re, mendrive_scalar_t* det_im
    );

    void K_E_v_eval(mendrive_scalar_t kz, mendrive_scalar_t sz,
        mendrive_complex_t *K_E_vacuum);
    void K_H_v_eval(mendrive_scalar_t kz, mendrive_scalar_t sz,
        mendrive_complex_t *K_H_vacuum);

    void K_E_l_eval(mendrive_scalar_t kz, mendrive_scalar_t sz,
        mendrive_complex_t *K_E_left_conductor);
    void K_H_l_eval(mendrive_scalar_t kz, mendrive_scalar_t sz,
        mendrive_complex_t *K_H_left_conductor);

    void K_E_r_eval(mendrive_scalar_t kz, mendrive_scalar_t sz,
        mendrive_complex_t *K_E_right_conductor) ;
    void K_H_r_eval(mendrive_scalar_t kz, mendrive_scalar_t sz,
        mendrive_complex_t *K_H_right_conductor);

    void det_eval(mendrive_scalar_t kz, mendrive_scalar_t sz,
        mendrive_scalar_t *det_re, mendrive_scalar_t *det_im);

    void det_div_diff_kz_eval(
        mendrive_scalar_t kz, mendrive_scalar_t sz,
        mendrive_scalar_t* div_re, mendrive_scalar_t* div_im
    );

    void det_div_diff_sz_eval(
        mendrive_scalar_t kz, mendrive_scalar_t sz,
        mendrive_scalar_t* div_re, mendrive_scalar_t* div_im
    );

    void det_derivatives(
        mendrive_scalar_t kz, mendrive_scalar_t sz,
        mendrive_scalar_t* dfdkz_re, mendrive_scalar_t* dfdkz_im,
        mendrive_scalar_t* dfdsz_re, mendrive_scalar_t* dfdsz_im
    );

#endif

// ============================================================================
// BLAS/LAPACK ИНТЕРФЕЙСЫ (для решения линейных систем в методе Ньютона)
// ============================================================================
#ifdef ARITH_FLOAT128
// Собственная мини-реализация базовых операций (нет стандартного __float128 BLAS)
void blas_dgemv_float128(char trans, int m, int n,
                         __float128 alpha, const __float128* A, int lda,
                         const __float128* x, int incx,
                         __float128 beta, __float128* y, int incy);
int blas_dgesv_float128(int n, int nrhs,
                        __float128* A, int lda,
                        int* ipiv,
                        __float128* B, int ldb,
                        int* info);

#elif defined(ARITH_MPFR_512)
// MPFR-версия решения 2×2 системы для метода Ньютона
int mpfr_linsolve_2x2(const mpfr_t a11, const mpfr_t a12,
                      const mpfr_t a21, const mpfr_t a22,
                      const mpfr_t b1,  const mpfr_t b2,
                      mpfr_t x1_out, mpfr_t x2_out,
                      mpfr_rnd_t rnd);

#else
// LONG DOUBLE: используем стандартный LAPACK через -llapack -lblas
// (long double → double преобразование внутри обёртки)
extern void dgesv_(int* n, int* nrhs, double* A, int* lda,
                   int* ipiv, double* B, int* ldb, int* info);
#endif

#endif // MENDDRIVE_DET_H
