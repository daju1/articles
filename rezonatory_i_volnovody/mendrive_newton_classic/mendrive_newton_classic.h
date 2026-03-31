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



/* ═══════════════════════════════════════════════════════════════════════════
 *  Public API (exported functions)
 * ═══════════════════════════════════════════════════════════════════════════ */

/**
 * Initialize classic Newton solver with default parameters
 * Call once before using other functions
 */
void newton_classic_init_default(void);

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
    int use_damping);

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
    int* niter_out);

/**
 * Solve with numeric Jacobian (finite differences)
 * Use this if det_derivatives is not accurate or not available
 */
int newton_classic_solve_numeric(
    mendrive_scalar_t* kz_inout,
    mendrive_scalar_t* sz_inout,
    mendrive_scalar_t* f_abs_out,
    int* niter_out);

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
    mendrive_scalar_t* det_J);
/**
 * Reset step iterator (call before starting new root search with newton_classic_step)
 */
void newton_classic_step_reset(void);

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
    mendrive_scalar_t* delta_im_out);

/**
 * Full solve with mendrive_scalar_t interface
 */
int solve_newton_root_classic(
    mendrive_scalar_t kz0, mendrive_scalar_t sz0,
    mendrive_scalar_t* kz_out, mendrive_scalar_t* sz_out,
    mendrive_scalar_t* f_abs_out,
    int max_iter);

/* ═══════════════════════════════════════════════════════════════════════════
 *  Debug / logging utilities
 * ═══════════════════════════════════════════════════════════════════════════ */

/**
 * Solve with verbose output
 */
int newton_classic_solve_verbose(
    mendrive_scalar_t* kz_inout,
    mendrive_scalar_t* sz_inout);

/**
 * Get string description of status code
 */
const char* newton_status_string(int status);
