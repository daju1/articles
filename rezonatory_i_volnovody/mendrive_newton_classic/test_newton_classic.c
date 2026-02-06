/**
 * @file test_newton_classic.c
 * @brief Test program for classic Newton solver
 *
 * Compile:
 *   gcc -O3 -DARITH_LONG_DOUBLE test_newton_classic.c -lm -o test_newton
 *   gcc -O3 -DARITH_FLOAT128 test_newton_classic.c -lquadmath -lm -o test_newton
 *
 * This test uses simple polynomial roots, not the full MenDrive det function.
 */

#include "../mendrive_det.h"
#include "newton_complex_classic.h"
#include <stdio.h>

/* ═══════════════════════════════════════════════════════════════════════════
 *  Test function 1: f(z) = z² - 1
 *  Roots: z = ±1
 * ═══════════════════════════════════════════════════════════════════════════ */

void test_func1(REAL x, REAL y, REAL* f_re, REAL* f_im, void* udata) {
    (void)udata;
    /* f(z) = z² - 1 = (x + iy)² - 1 = (x² - y² - 1) + i(2xy) */
    *f_re = x * x - y * y - REAL_LITERAL(1.0);
    *f_im = REAL_LITERAL(2.0) * x * y;
}

void test_jac1(REAL x, REAL y,
               REAL* J11, REAL* J12, REAL* J21, REAL* J22,
               void* udata) {
    (void)udata;
    /* ∂f_re/∂x = 2x, ∂f_re/∂y = -2y */
    /* ∂f_im/∂x = 2y, ∂f_im/∂y = 2x */
    *J11 = REAL_LITERAL(2.0) * x;
    *J12 = REAL_LITERAL(-2.0) * y;
    *J21 = REAL_LITERAL(2.0) * y;
    *J22 = REAL_LITERAL(2.0) * x;
}

/* ═══════════════════════════════════════════════════════════════════════════
 *  Test function 2: f(z) = z³ - 1
 *  Roots: z = 1, e^(2πi/3), e^(4πi/3)
 * ═══════════════════════════════════════════════════════════════════════════ */

void test_func2(REAL x, REAL y, REAL* f_re, REAL* f_im, void* udata) {
    (void)udata;
    /* f(z) = z³ - 1 */
    REAL x2 = x * x;
    REAL y2 = y * y;
    REAL x3 = x * x2;
    REAL y3 = y * y2;

    *f_re = x3 - REAL_LITERAL(3.0) * x * y2 - REAL_LITERAL(1.0);
    *f_im = REAL_LITERAL(3.0) * x2 * y - y3;
}

/* ═══════════════════════════════════════════════════════════════════════════
 *  Test function 3: f(z) = sin(z) (complex sine)
 *  Roots: z = nπ for integer n
 * ═══════════════════════════════════════════════════════════════════════════ */

void test_func_sin(REAL x, REAL y, REAL* f_re, REAL* f_im, void* udata) {
    (void)udata;
    /* sin(z) = sin(x)cosh(y) + i cos(x)sinh(y) */
#ifdef ARITH_FLOAT128
    *f_re = sinq(x) * coshq(y);
    *f_im = cosq(x) * sinhq(y);
#elif defined(ARITH_LONG_DOUBLE)
    *f_re = sinl(x) * coshl(y);
    *f_im = cosl(x) * sinhl(y);
#else
    *f_re = sin(x) * cosh(y);
    *f_im = cos(x) * sinh(y);
#endif
}

/* ═══════════════════════════════════════════════════════════════════════════
 *  Main test
 * ═══════════════════════════════════════════════════════════════════════════ */

int main() {
    newton_params_t params;
    newton_state_t state;

#ifdef ARITH_FLOAT128
    printf("Testing Newton solver with __float128 precision\n");
    newton_params_highprec(&params);
#elif defined(ARITH_LONG_DOUBLE)
    printf("Testing Newton solver with long double precision\n");
    newton_params_highprec(&params);
#else
    printf("Testing Newton solver with double precision\n");
    newton_params_default(&params);
#endif

    printf("─────────────────────────────────────────────────────\n");

    /* Test 1: z² - 1 = 0, starting from (0.5, 0.5) */
    {
        REAL x = REAL_LITERAL(0.5);
        REAL y = REAL_LITERAL(0.5);

        printf("\nTest 1: f(z) = z² - 1, starting at (%.2f, %.2f)\n", (double)x, (double)y);

        int status = newton_complex_solve(&x, &y, test_func1, test_jac1,
                                          &params, NULL, &state);

        printf("  Result: x=" REAL_FMT ", y=" REAL_FMT "\n", x, y);
        printf("  |f| = " REAL_FMT ", iterations = %d\n", state.f_abs, state.iteration);
        printf("  Status: %d (expected root: 1 + 0i)\n", status);
    }

    /* Test 2: z² - 1 = 0, starting from (-0.5, 0.5) */
    {
        REAL x = REAL_LITERAL(-0.5);
        REAL y = REAL_LITERAL(0.5);

        printf("\nTest 2: f(z) = z² - 1, starting at (%.2f, %.2f)\n", (double)x, (double)y);

        int status = newton_complex_solve(&x, &y, test_func1, test_jac1,
                                          &params, NULL, &state);

        printf("  Result: x=" REAL_FMT ", y=" REAL_FMT "\n", x, y);
        printf("  |f| = " REAL_FMT ", iterations = %d\n", state.f_abs, state.iteration);
        printf("  Status: %d (expected root: -1 + 0i)\n", status);
    }

    /* Test 3: z³ - 1 = 0 with numeric Jacobian */
    {
        REAL x = REAL_LITERAL(0.5);
        REAL y = REAL_LITERAL(0.8);

        printf("\nTest 3: f(z) = z³ - 1, numeric Jacobian, starting at (%.2f, %.2f)\n",
               (double)x, (double)y);

        int status = newton_complex_solve(&x, &y, test_func2, NULL,  /* NULL = numeric */
                                          &params, NULL, &state);

        printf("  Result: x=" REAL_FMT ", y=" REAL_FMT "\n", x, y);
        printf("  |f| = " REAL_FMT ", iterations = %d\n", state.f_abs, state.iteration);
        printf("  Status: %d (expected root: e^(2πi/3) ≈ -0.5 + 0.866i)\n", status);
    }

    /* Test 4: sin(z) = 0, find π */
    {
        REAL x = REAL_LITERAL(3.0);
        REAL y = REAL_LITERAL(0.1);

        printf("\nTest 4: f(z) = sin(z), numeric Jacobian, starting at (%.2f, %.2f)\n",
               (double)x, (double)y);

        int status = newton_complex_solve(&x, &y, test_func_sin, NULL,
                                          &params, NULL, &state);

        printf("  Result: x=" REAL_FMT ", y=" REAL_FMT "\n", x, y);
        printf("  |f| = " REAL_FMT ", iterations = %d\n", state.f_abs, state.iteration);
#ifdef ARITH_FLOAT128
        printf("  Status: %d (expected root: π ≈ %.20Qf)\n", status, M_PIq);
#else
        printf("  Status: %d (expected root: π ≈ %.15f)\n", status, (double)3.14159265358979323846L);
#endif
    }

    /* Test 5: Step-by-step iteration */
    {
        REAL x = REAL_LITERAL(2.0);
        REAL y = REAL_LITERAL(0.0);

        printf("\nTest 5: Step-by-step iteration for z² - 1 = 0\n");
        printf("  Starting at (%.2f, %.2f)\n", (double)x, (double)y);

        newton_complex_init(&state, x, y, test_func1, &params, NULL);
        printf("  Initial |f| = " REAL_FMT "\n", state.f_abs);

        for (int i = 0; i < 10; i++) {
            int status = newton_complex_step(&state, test_func1, test_jac1,
                                              &params, NULL);
            printf("  [%d] x=" REAL_FMT ", y=" REAL_FMT ", |f|=" REAL_FMT "\n",
                   i+1, state.x, state.y, state.f_abs);

            if (status != NEWTON_CONTINUE) {
                printf("  Converged with status %d\n", status);
                break;
            }
        }
    }

    printf("\n─────────────────────────────────────────────────────\n");
    printf("All tests completed.\n");

    return 0;
}
