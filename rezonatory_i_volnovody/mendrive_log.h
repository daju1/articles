#ifndef MENDDRIVE_LOG_H
#define MENDDRIVE_LOG_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "mendrive_det.h"

// ============================================================================
// УРОВНИ ЛОГГИРОВАНИЯ
// ============================================================================
#define MPREC_LOG_LEVEL_OFF    0
#define MPREC_LOG_LEVEL_ERROR  1
#define MPREC_LOG_LEVEL_WARN   2
#define MPREC_LOG_LEVEL_INFO   3
#define MPREC_LOG_LEVEL_DEBUG  4
#define MPREC_LOG_LEVEL_TRACE  5

// Выберите уровень логгирования при компиляции:
//   -DMPREC_LOG_LEVEL=MPREC_LOG_LEVEL_DEBUG  (по умолчанию)
//   -DMPREC_LOG_LEVEL=MPREC_LOG_LEVEL_INFO
//   -DMPREC_LOG_LEVEL=MPREC_LOG_LEVEL_OFF    (полностью отключено)
#ifndef MPREC_LOG_LEVEL
    #define MPREC_LOG_LEVEL MPREC_LOG_LEVEL_INFO
#endif

// ============================================================================
// ЦЕЛЕВОЙ ВЫВОД
// ============================================================================
#ifndef MPREC_LOG_TARGET
    #define MPREC_LOG_TARGET stderr  // stderr по умолчанию
    // Альтернативы:
    // #define MPREC_LOG_TARGET stdout
    // #define MPREC_LOG_TARGET fopen("/tmp/mendrive.log", "a")
#endif

// ============================================================================
// ФОРМАТИРОВАНИЕ ЧИСЕЛ В ЗАВИСИМОСТИ ОТ РЕЖИМА ТОЧНОСТИ
// ============================================================================
#ifdef ARITH_MPFR_512
    #include <mpfr.h>
    #define MPREC_LOG_SNPRINTF(buf, size, fmt, ...) \
        mpfr_snprintf(buf, size, fmt, ##__VA_ARGS__)
    #define MPREC_LOG_FMT_SCALAR "%.30Rf"
    #define MPREC_LOG_FMT_COMPLEX "(%.30Rf, %.30Rf)"
#elif defined(ARITH_FLOAT128)
    #include <quadmath.h>
    // quadmath_snprintf — правильная функция для __float128 (НЕ snprintfq!)
    #define MPREC_LOG_SNPRINTF(buf, size, fmt, ...) \
        quadmath_snprintf(buf, size, fmt, ##__VA_ARGS__)
    #define MPREC_LOG_FMT_SCALAR "%.34Qg"
    #define MPREC_LOG_FMT_COMPLEX "(%.34Qg, %.34Qg)"
#else  // ARITH_LONG_DOUBLE
    #define MPREC_LOG_SNPRINTF(buf, size, fmt, ...) \
        snprintf(buf, size, fmt, ##__VA_ARGS__)
    #define MPREC_LOG_FMT_SCALAR "%.19Lg"
    #define MPREC_LOG_FMT_COMPLEX "(%.19Lg, %.19Lg)"
#endif

// ============================================================================
// ВСПОМОГАТЕЛЬНЫЕ ФУНКЦИИ
// ============================================================================
static inline void mprec_log_timestamp(char* buf, size_t size) {
    time_t t = time(NULL);
    struct tm* tm_info = localtime(&t);
    strftime(buf, size, "%Y-%m-%d %H:%M:%S", tm_info);
}

static inline const char* mprec_log_level_str(int level) {
    switch (level) {
        case MPREC_LOG_LEVEL_TRACE: return "TRACE";
        case MPREC_LOG_LEVEL_DEBUG: return "DEBUG";
        case MPREC_LOG_LEVEL_INFO:  return "INFO";
        case MPREC_LOG_LEVEL_WARN:  return "WARN";
        case MPREC_LOG_LEVEL_ERROR: return "ERROR";
        default: return "UNKNOWN";
    }
}

// ============================================================================
// ОСНОВНОЙ МАКРОС ЛОГГИРОВАНИЯ
// ============================================================================
#define MPREC_LOG(level, fmt, ...) \
    do { \
        if (level <= MPREC_LOG_LEVEL) { \
            char ts_buf[20]; \
            char log_buf[1024]; \
            mprec_log_timestamp(ts_buf, sizeof(ts_buf)); \
            MPREC_LOG_SNPRINTF(log_buf, sizeof(log_buf), fmt, ##__VA_ARGS__); \
            fprintf(MPREC_LOG_TARGET, \
                "[%s] [%s] [%s:%d] " fmt "\n", \
                ts_buf, \
                mprec_log_level_str(level), \
                __FILE__, __LINE__, \
                ##__VA_ARGS__); \
            fflush(MPREC_LOG_TARGET); \
        } \
    } while(0)

// Удобные алиасы
#define MPREC_LOG_TRACE(fmt, ...) MPREC_LOG(MPREC_LOG_LEVEL_TRACE, fmt, ##__VA_ARGS__)
#define MPREC_LOG_DEBUG(fmt, ...) MPREC_LOG(MPREC_LOG_LEVEL_DEBUG, fmt, ##__VA_ARGS__)
#define MPREC_LOG_INFO(fmt, ...)  MPREC_LOG(MPREC_LOG_LEVEL_INFO,  fmt, ##__VA_ARGS__)
#define MPREC_LOG_WARN(fmt, ...)  MPREC_LOG(MPREC_LOG_LEVEL_WARN,  fmt, ##__VA_ARGS__)
#define MPREC_LOG_ERROR(fmt, ...) MPREC_LOG(MPREC_LOG_LEVEL_ERROR, fmt, ##__VA_ARGS__)

// ============================================================================
// СПЕЦИАЛИЗИРОВАННЫЕ ФУНКЦИИ ДЛЯ ЧИСЕЛ (безопасные для всех режимов)
// ============================================================================
#ifdef ARITH_MPFR_512
    static inline void mprec_log_scalar(const char* name, mpfr_t x, int level) {
        if (level <= MPREC_LOG_LEVEL) {
            char buf[256];
            mpfr_snprintf(buf, sizeof(buf), MPREC_LOG_FMT_SCALAR, x);
            MPREC_LOG(level, "%s = %s", name, buf);
        }
    }

    static inline void mprec_log_complex(const char* name, mpfr_t re, mpfr_t im, int level) {
        if (level <= MPREC_LOG_LEVEL) {
            char buf_re[256], buf_im[256];
            mpfr_snprintf(buf_re, sizeof(buf_re), MPREC_LOG_FMT_SCALAR, re);
            mpfr_snprintf(buf_im, sizeof(buf_im), MPREC_LOG_FMT_SCALAR, im);
            MPREC_LOG(level, "%s = (%s, %s)", name, buf_re, buf_im);
        }
    }

#elif defined(ARITH_FLOAT128)
    static inline void mprec_log_scalar(const char* name, __float128 x, int level) {
        if (level <= MPREC_LOG_LEVEL) {
            char buf[256];
            quadmath_snprintf(buf, sizeof(buf), MPREC_LOG_FMT_SCALAR, x);
            MPREC_LOG(level, "%s = %s", name, buf);
        }
    }

    static inline void mprec_log_complex(const char* name, __complex128 z, int level) {
        if (level <= MPREC_LOG_LEVEL) {
            char buf_re[256], buf_im[256];
            quadmath_snprintf(buf_re, sizeof(buf_re), MPREC_LOG_FMT_SCALAR, crealq(z));
            quadmath_snprintf(buf_im, sizeof(buf_im), MPREC_LOG_FMT_SCALAR, cimagq(z));
            MPREC_LOG(level, "%s = (%s, %s)", name, buf_re, buf_im);
        }
    }

#else  // ARITH_LONG_DOUBLE
    static inline void mprec_log_scalar(const char* name, long double x, int level) {
        if (level <= MPREC_LOG_LEVEL) {
            char buf[256];
            snprintf(buf, sizeof(buf), MPREC_LOG_FMT_SCALAR, x);
            MPREC_LOG(level, "%s = %s", name, buf);
        }
    }

    static inline void mprec_log_complex(const char* name, long double _Complex z, int level) {
        if (level <= MPREC_LOG_LEVEL) {
            char buf_re[256], buf_im[256];
            snprintf(buf_re, sizeof(buf_re), MPREC_LOG_FMT_SCALAR, creall(z));
            snprintf(buf_im, sizeof(buf_im), MPREC_LOG_FMT_SCALAR, cimagl(z));
            MPREC_LOG(level, "%s = (%s, %s)", name, buf_re, buf_im);
        }
    }
#endif

#define MPREC_LOG_SCALAR_DEBUG(name, x) mprec_log_scalar(name, x, MPREC_LOG_LEVEL_DEBUG)
#define MPREC_LOG_SCALAR_INFO(name, x)  mprec_log_scalar(name, x, MPREC_LOG_LEVEL_INFO)
#define MPREC_LOG_COMPLEX_DEBUG(name, re, im) mprec_log_complex(name, re, im, MPREC_LOG_LEVEL_DEBUG)
#define MPREC_LOG_COMPLEX_INFO(name, re, im)  mprec_log_complex(name, re, im, MPREC_LOG_LEVEL_INFO)

// ============================================================================
// УДОБНЫЙ МАКРОС ДЛЯ ЛОГГИРОВАНИЯ ПАРАМЕТРОВ В det_init()
// ============================================================================
#ifdef ARITH_FLOAT128
    #define MPREC_LOG_PARAM(name, value) \
        do { \
            char buf[256]; \
            quadmath_snprintf(buf, sizeof(buf), "%.34Qg", value); \
            fprintf(MPREC_LOG_TARGET, "  %s = %s\n", name, buf); \
        } while(0)
#elif defined(ARITH_MPFR_512)
    #define MPREC_LOG_PARAM(name, value) \
        do { \
            char buf[256]; \
            mpfr_snprintf(buf, sizeof(buf), "%.30Rf", value); \
            fprintf(MPREC_LOG_TARGET, "  %s = %s\n", name, buf); \
        } while(0)
#else
    #define MPREC_LOG_PARAM(name, value) \
        fprintf(MPREC_LOG_TARGET, "  %s = %.19Lg\n", name, (long double)(value))
#endif

#endif  // MENDDRIVE_LOG_H
