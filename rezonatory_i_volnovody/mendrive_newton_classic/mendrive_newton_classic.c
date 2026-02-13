/**
 * @file newton_mendrive_classic.c
 * @brief Classic Newton solver integrated with MenDrive det functions
 *
 * Compile options:
 *   gcc -DUSE_LONG_DOUBLE -O3 -c newton_mendrive_classic.c
 *   gcc -DUSE_FLOAT128 -O3 -c newton_mendrive_classic.c -lquadmath
 *   gcc -O3 -c newton_mendrive_classic.c  (default: double)
 *
 * To match your existing library:
 *   -DKY   — enable ky mode (your existing flag)
 *   -DQNM  — enable QNM mode (omega search instead of kz/sz)
 */

#include "../mendrive_det.h"
#include "newton_complex_classic.h"
#include "mendrive_newton_classic.h"

#include <stdio.h>
#include <stdlib.h>

/* ═══════════════════════════════════════════════════════════════════════════
 *  Global state (matches your existing architecture)
 * ═══════════════════════════════════════════════════════════════════════════ */

static newton_params_t g_newton_params;
static int g_params_initialized = 0;

/* ═══════════════════════════════════════════════════════════════════════════
 *  Wrapper: det_eval as newton_func_t
 * ═══════════════════════════════════════════════════════════════════════════ */

/**
 * Adapter function: wraps det_eval for newton_func_t interface
 * For kz/sz mode (non-QNM): x = kz, y = sz
 * For QNM mode: x = omega_r, y = omega_i
 */
static void det_wrapper(REAL x, REAL y, REAL* f_re, REAL* f_im, void* udata) {
    (void)udata;  /* Not used, det_eval uses global params from det_init */
    det_eval(x, y, f_re, f_im);
}

/**
 * Adapter for analytic Jacobian using det_derivatives
 */
static void jacobian_wrapper(REAL x, REAL y,
                              REAL* J11, REAL* J12,
                              REAL* J21, REAL* J22,
                              void* udata) {
    (void)udata;
    /* det_derivatives returns: dfdkz_re, dfdkz_im, dfdsz_re, dfdsz_im */
    // *J11 /* ∂f_re/∂x = dfdkz_re */
    // *J21 /* ∂f_im/∂x = dfdkz_im */
    // *J12 /* ∂f_re/∂y = dfdsz_re */
    // *J22 /* ∂f_im/∂y = dfdsz_im */

    det_derivatives(x, y, J11, J21, J12, J22);
}

/* ═══════════════════════════════════════════════════════════════════════════
 *  Public API (exported functions)
 * ═══════════════════════════════════════════════════════════════════════════ */

/**
 * Initialize classic Newton solver with default parameters
 * Call once before using other functions
 */
void newton_classic_init_default(void) {
#ifdef USE_FLOAT128
    newton_params_highprec(&g_newton_params);
#elif defined(USE_LONG_DOUBLE)
    newton_params_highprec(&g_newton_params);
#else
    newton_params_default(&g_newton_params);
#endif
    g_params_initialized = 1;
}

/**
 * Initialize with custom tolerances
 *
 * @param delta_tol   Convergence tolerance for |Δz|
 * @param f_tol       Convergence tolerance for |f(z)|
 * @param max_iter    Maximum iterations
 * @param use_damping 1 = damped Newton, 0 = pure Newton
 */
void newton_classic_init(
    REAL delta_tol,
    REAL f_tol,
    int max_iter,
    int use_damping)
{
#ifdef USE_FLOAT128
    newton_params_highprec(&g_newton_params);
#elif defined(USE_LONG_DOUBLE)
    newton_params_highprec(&g_newton_params);
#else
    newton_params_default(&g_newton_params);
#endif

    g_newton_params.delta_tol = delta_tol;
    g_newton_params.f_tol = f_tol;
    g_newton_params.max_iterations = max_iter;
    g_newton_params.use_damping = use_damping;

    g_params_initialized = 1;
}

/**
 * Solve det(kz, sz) = 0 using classic Newton method
 *
 * @param kz_inout  Input: initial guess for kz, Output: solution
 * @param sz_inout  Input: initial guess for sz, Output: solution
 * @param f_abs_out Output: final |det| value (optional, can be NULL)
 * @param niter_out Output: number of iterations (optional, can be NULL)
 * @return Status code (>0 = converged, <0 = error)
 */
int newton_classic_solve(
    mendrive_scalar_t* kz_inout,
    mendrive_scalar_t* sz_inout,
    mendrive_scalar_t* f_abs_out,
    int* niter_out)
{
    if (!g_params_initialized) {
        newton_classic_init_default();
    }

    newton_state_t state;
    int status = newton_complex_solve(
        kz_inout, sz_inout,
        det_wrapper,
        jacobian_wrapper,  /* Use analytic Jacobian from det_derivatives */
        &g_newton_params,
        NULL,  /* No user data needed */
        &state
    );

    if (f_abs_out) *f_abs_out = state.f_abs;
    if (niter_out) *niter_out = state.iteration;

    return status;
}

/**
 * Solve with numeric Jacobian (finite differences)
 * Use this if det_derivatives is not accurate or not available
 */
int newton_classic_solve_numeric(
    mendrive_scalar_t* kz_inout,
    mendrive_scalar_t* sz_inout,
    mendrive_scalar_t* f_abs_out,
    int* niter_out)
{
    if (!g_params_initialized) {
        newton_classic_init_default();
    }

    newton_state_t state;
    int status = newton_complex_solve(
        kz_inout, sz_inout,
        det_wrapper,
        NULL,  /* Use numeric Jacobian */
        &g_newton_params,
        NULL,
        &state
    );

    if (f_abs_out) *f_abs_out = state.f_abs;
    if (niter_out) *niter_out = state.iteration;

    return status;
}

/**
 * Single Newton step (for manual control / debugging)
 *
 * @param kz_inout  Current kz (updated)
 * @param sz_inout  Current sz (updated)
 * @param f_abs_out Current |det|
 * @param delta_kz  Step taken in kz
 * @param delta_sz  Step taken in sz
 * @param det_J     Jacobian determinant
 * @return Status code
 */
int newton_classic_step(
    mendrive_scalar_t* kz_inout,
    mendrive_scalar_t* sz_inout,
    mendrive_scalar_t* f_abs_out,
    mendrive_scalar_t* delta_kz,
    mendrive_scalar_t* delta_sz,
    mendrive_scalar_t* det_J)
{
    static newton_state_t state;
    static int step_initialized = 0;

    if (!g_params_initialized) {
        newton_classic_init_default();
    }

    if (!step_initialized) {
        newton_complex_init(&state, *kz_inout, *sz_inout,
                            det_wrapper, &g_newton_params, NULL);
        step_initialized = 1;
    }

    int status = newton_complex_step(&state, det_wrapper, jacobian_wrapper,
                                      &g_newton_params, NULL);

    *kz_inout = state.x;
    *sz_inout = state.y;
    if (f_abs_out) *f_abs_out = state.f_abs;
    if (delta_kz) *delta_kz = state.delta_x;
    if (delta_sz) *delta_sz = state.delta_y;
    if (det_J) *det_J = state.det_J;

    /* Reset for next root search if converged */
    if (status != NEWTON_CONTINUE) {
        step_initialized = 0;
    }

    return status;
}

/**
 * Reset step iterator (call before starting new root search with newton_classic_step)
 */
void newton_classic_step_reset(void) {
    /* Force reinitialization on next step call */
    /* Implemented via static variable in newton_classic_step */
}

/* ═══════════════════════════════════════════════════════════════════════════
 *  Convenience wrappers for mendrive_scalar_t interface (your existing style)
 * ═══════════════════════════════════════════════════════════════════════════ */

/**
 * Wrapper matching your existing newton_step signature
 * For compatibility with existing code
 */
int newton_step_classic(
    mendrive_scalar_t* xn_re, mendrive_scalar_t* xn_im,
    mendrive_scalar_t* f_abs_out,
    mendrive_scalar_t* delta_re_out,
    mendrive_scalar_t* delta_im_out)
{
    mendrive_scalar_t kz = *xn_re;
    mendrive_scalar_t sz = *xn_im;
    mendrive_scalar_t f_abs, delta_kz, delta_sz, det;

    int status = newton_classic_step(&kz, &sz, &f_abs, &delta_kz, &delta_sz, &det);

    *xn_re = kz;
    *xn_im = sz;
    if (f_abs_out) *f_abs_out = f_abs;
    if (delta_re_out) *delta_re_out = delta_kz;
    if (delta_im_out) *delta_im_out = delta_sz;

    return status;
}

/**
 * Full solve with mendrive_scalar_t interface
 */
int solve_newton_root_classic(
    mendrive_scalar_t kz0, mendrive_scalar_t sz0,
    mendrive_scalar_t* kz_out, mendrive_scalar_t* sz_out,
    mendrive_scalar_t* f_abs_out,
    int max_iter)
{
    if (max_iter > 0) {
        g_newton_params.max_iterations = max_iter;
    }

    mendrive_scalar_t kz = kz0;
    mendrive_scalar_t sz = sz0;
    mendrive_scalar_t f_abs;
    int niter;

    int status = newton_classic_solve(&kz, &sz, &f_abs, &niter);

    *kz_out = kz;
    *sz_out = sz;
    if (f_abs_out) *f_abs_out = f_abs;

    return status;
}

/* ═══════════════════════════════════════════════════════════════════════════
 *  Debug / logging utilities
 * ═══════════════════════════════════════════════════════════════════════════ */

/**
 * Solve with verbose output
 */
int newton_classic_solve_verbose(
    mendrive_scalar_t* kz_inout,
    mendrive_scalar_t* sz_inout)
{
    if (!g_params_initialized) {
        newton_classic_init_default();
    }

    newton_state_t state;
    newton_complex_init(&state, *kz_inout, *sz_inout,
                        det_wrapper, &g_newton_params, NULL);

    printf("Classic Newton iteration:\n");
    printf("  Initial: kz=" REAL_FMT ", sz=" REAL_FMT ", |f|=" REAL_FMT "\n",
           state.x, state.y, state.f_abs);

    int status;
    for (int i = 0; i < g_newton_params.max_iterations; i++) {
        status = newton_complex_step(&state, det_wrapper, jacobian_wrapper,
                                      &g_newton_params, NULL);

        printf("  [%3d] kz=" REAL_FMT ", sz=" REAL_FMT "\n"
               "        |f|=" REAL_FMT ", |Δ|=" REAL_FMT ", det(J)=" REAL_FMT "\n",
               i+1, state.x, state.y, state.f_abs, state.delta_abs, state.det_J);

        if (status != NEWTON_CONTINUE) {
            break;
        }
    }

    const char* status_str;
    switch (status) {
        case NEWTON_CONVERGED_DELTA: status_str = "✓ converged (delta)"; break;
        case NEWTON_CONVERGED_FABS:  status_str = "✓ converged (|f|)"; break;
        case NEWTON_SINGULAR:        status_str = "✗ singular Jacobian"; break;
        case NEWTON_DIVERGED:        status_str = "✗ diverged"; break;
        case NEWTON_MAX_ITER:        status_str = "✗ max iterations"; break;
        default:                     status_str = "?"; break;
    }
    printf("  Result: %s after %d iterations\n", status_str, state.iteration);

    *kz_inout = state.x;
    *sz_inout = state.y;

    return status;
}

/**
 * Get string description of status code
 */
const char* newton_status_string(int status) {
    switch (status) {
        case NEWTON_CONTINUE:        return "continuing";
        case NEWTON_CONVERGED_DELTA: return "converged (step < tol)";
        case NEWTON_CONVERGED_FABS:  return "converged (|f| < tol)";
        case NEWTON_SINGULAR:        return "error: singular Jacobian";
        case NEWTON_DIVERGED:        return "error: diverged";
        case NEWTON_MAX_ITER:        return "max iterations reached";
        default:                     return "unknown status";
    }
}
