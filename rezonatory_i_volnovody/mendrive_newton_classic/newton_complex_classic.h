/**
 * @file newton_complex_classic.h
 * @brief Classic complex Newton method for finding roots of f(z)=0
 *
 * Configurable precision via compile-time macros:
 *   -DUSE_FLOAT128    — __float128 (quadruple precision, requires libquadmath)
 *   -DUSE_LONG_DOUBLE — long double (extended precision, 80-bit on x86)
 *   default           — double (64-bit IEEE 754)
 *
 * Usage:
 *   1. Define your function: void my_func(REAL x, REAL y, REAL* f_re, REAL* f_im)
 *   2. Call newton_complex_init() with parameters
 *   3. Call newton_complex_solve() to find root
 *   4. Or use newton_complex_step() for manual iteration control
 *
 * Author: Generated for MenDrive project
 */

#ifndef NEWTON_COMPLEX_CLASSIC_H
#define NEWTON_COMPLEX_CLASSIC_H

#include <math.h>

/* ═══════════════════════════════════════════════════════════════════════════
 *  Precision selection
 * ═══════════════════════════════════════════════════════════════════════════ */

#ifdef ARITH_FLOAT128
    #include <quadmath.h>
    typedef __float128 REAL;
    #define REAL_SQRT  sqrtq
    #define REAL_FABS  fabsq
    #define REAL_EPS   FLT128_EPSILON
    #define REAL_FMT   "%.36Qe"
    #define REAL_LITERAL(x) (x##Q)
#elif defined(ARITH_LONG_DOUBLE)
    typedef long double REAL;
    #define REAL_SQRT  sqrtl
    #define REAL_FABS  fabsl
    #define REAL_EPS   LDBL_EPSILON
    #define REAL_FMT   "%.21Le"
    #define REAL_LITERAL(x) (x##L)
#else
    typedef double REAL;
    #define REAL_SQRT  sqrt
    #define REAL_FABS  fabs
    #define REAL_EPS   2.2204460492503131e-16
    #define REAL_FMT   "%.15e"
    #define REAL_LITERAL(x) (x)
#endif

/* ═══════════════════════════════════════════════════════════════════════════
 *  Return codes
 * ═══════════════════════════════════════════════════════════════════════════ */

#define NEWTON_CONTINUE        0   /* Iteration continues */
#define NEWTON_CONVERGED_DELTA 1   /* |Δz| < delta_tol */
#define NEWTON_CONVERGED_FABS  2   /* |f(z)| < f_tol */
#define NEWTON_SINGULAR       -1   /* Jacobian singular (det ≈ 0) */
#define NEWTON_DIVERGED       -2   /* |f| increased too much */
#define NEWTON_MAX_ITER       -3   /* Maximum iterations reached */

/* ═══════════════════════════════════════════════════════════════════════════
 *  Structures
 * ═══════════════════════════════════════════════════════════════════════════ */

/**
 * Parameters for the Newton solver
 */
typedef struct {
    /* Tolerances */
    REAL delta_tol;      /**< Convergence tolerance for |Δz| */
    REAL f_tol;          /**< Convergence tolerance for |f(z)| */
    REAL jacobian_tol;   /**< Threshold for singular Jacobian detection */

    /* Adaptive step parameters (optional) */
    int  use_damping;    /**< 1 = use damped Newton, 0 = pure Newton */
    REAL damping_init;   /**< Initial damping factor (0 < λ ≤ 1) */
    REAL damping_min;    /**< Minimum damping factor */
    REAL damping_decrease; /**< Factor to decrease damping on bad step */
    REAL damping_increase; /**< Factor to increase damping on good step */

    /* Numerical differentiation step (if analytic Jacobian not provided) */
    REAL diff_h;         /**< Step for finite differences */

    /* Limits */
    int max_iterations;  /**< Maximum number of iterations */
    int max_backsteps;   /**< Maximum backsteps in damped mode */
} newton_params_t;

/**
 * State of the Newton iteration
 */
typedef struct {
    /* Current point */
    REAL x;              /**< Real part of z = x + iy */
    REAL y;              /**< Imaginary part of z = x + iy */

    /* Function value at current point */
    REAL f_re;           /**< Re(f(z)) */
    REAL f_im;           /**< Im(f(z)) */
    REAL f_abs;          /**< |f(z)| */

    /* Last step */
    REAL delta_x;        /**< Δx from last iteration */
    REAL delta_y;        /**< Δy from last iteration */
    REAL delta_abs;      /**< |Δz| from last iteration */

    /* Jacobian at current point */
    REAL J11, J12;       /**< ∂f_re/∂x, ∂f_re/∂y */
    REAL J21, J22;       /**< ∂f_im/∂x, ∂f_im/∂y */
    REAL det_J;          /**< Determinant of Jacobian */

    /* Damping state */
    REAL damping;        /**< Current damping factor */

    /* Iteration counter */
    int iteration;       /**< Current iteration number */
    int status;          /**< Last return code */
} newton_state_t;

/**
 * Function type: computes f(z) = f_re + i*f_im at point z = x + i*y
 *
 * @param x      Real part of z
 * @param y      Imaginary part of z
 * @param f_re   Output: real part of f(z)
 * @param f_im   Output: imaginary part of f(z)
 * @param udata  User data pointer (passed through)
 */
typedef void (*newton_func_t)(REAL x, REAL y, REAL* f_re, REAL* f_im, void* udata);

/**
 * Jacobian function type (optional, for analytic Jacobian)
 * Computes the 2x2 Jacobian matrix:
 *   J = [ ∂f_re/∂x  ∂f_re/∂y ]
 *       [ ∂f_im/∂x  ∂f_im/∂y ]
 *
 * @param x      Real part of z
 * @param y      Imaginary part of z
 * @param J11    Output: ∂f_re/∂x
 * @param J12    Output: ∂f_re/∂y
 * @param J21    Output: ∂f_im/∂x
 * @param J22    Output: ∂f_im/∂y
 * @param udata  User data pointer
 */
typedef void (*newton_jacobian_t)(REAL x, REAL y,
                                   REAL* J11, REAL* J12,
                                   REAL* J21, REAL* J22,
                                   void* udata);

/* ═══════════════════════════════════════════════════════════════════════════
 *  Default parameters
 * ═══════════════════════════════════════════════════════════════════════════ */

/**
 * Initialize parameters with reasonable defaults
 */
static inline void newton_params_default(newton_params_t* p) {
    p->delta_tol       = REAL_LITERAL(1e-12);
    p->f_tol           = REAL_LITERAL(1e-12);
    p->jacobian_tol    = REAL_LITERAL(1e-30);

    p->use_damping     = 1;
    p->damping_init    = REAL_LITERAL(1.0);
    p->damping_min     = REAL_LITERAL(1e-6);
    p->damping_decrease = REAL_LITERAL(0.5);
    p->damping_increase = REAL_LITERAL(1.2);

    p->diff_h          = REAL_LITERAL(1e-8);

    p->max_iterations  = 100;
    p->max_backsteps   = 20;
}

/**
 * High-precision defaults (for long double / float128)
 */
static inline void newton_params_highprec(newton_params_t* p) {
    p->delta_tol       = REAL_LITERAL(1e-28);
    p->f_tol           = REAL_LITERAL(1e-28);
    p->jacobian_tol    = REAL_LITERAL(1e-60);

    p->use_damping     = 1;
    p->damping_init    = REAL_LITERAL(1.0);
    p->damping_min     = REAL_LITERAL(1e-10);
    p->damping_decrease = REAL_LITERAL(0.5);
    p->damping_increase = REAL_LITERAL(1.2);

    p->diff_h          = REAL_LITERAL(1e-15);

    p->max_iterations  = 200;
    p->max_backsteps   = 30;
}

/* ═══════════════════════════════════════════════════════════════════════════
 *  Core functions
 * ═══════════════════════════════════════════════════════════════════════════ */

/**
 * Initialize Newton state with starting point
 *
 * @param state  Output state structure
 * @param x0     Initial guess: real part
 * @param y0     Initial guess: imaginary part
 * @param func   Function to find root of
 * @param params Solver parameters
 * @param udata  User data for function
 */
static inline void newton_complex_init(
    newton_state_t* state,
    REAL x0, REAL y0,
    newton_func_t func,
    const newton_params_t* params,
    void* udata)
{
    state->x = x0;
    state->y = y0;
    state->delta_x = 0;
    state->delta_y = 0;
    state->delta_abs = 0;
    state->damping = params->damping_init;
    state->iteration = 0;
    state->status = NEWTON_CONTINUE;

    /* Evaluate f at initial point */
    func(x0, y0, &state->f_re, &state->f_im, udata);
    state->f_abs = REAL_SQRT(state->f_re * state->f_re + state->f_im * state->f_im);
}

/**
 * Compute Jacobian via central finite differences
 */
static inline void newton_jacobian_numeric(
    REAL x, REAL y,
    REAL* J11, REAL* J12, REAL* J21, REAL* J22,
    newton_func_t func,
    REAL h,
    void* udata)
{
    REAL f_p_re, f_p_im, f_m_re, f_m_im;
    REAL inv_2h = REAL_LITERAL(1.0) / (REAL_LITERAL(2.0) * h);

    /* ∂f/∂x via central difference */
    func(x + h, y, &f_p_re, &f_p_im, udata);
    func(x - h, y, &f_m_re, &f_m_im, udata);
    *J11 = (f_p_re - f_m_re) * inv_2h;
    *J21 = (f_p_im - f_m_im) * inv_2h;

    /* ∂f/∂y via central difference */
    func(x, y + h, &f_p_re, &f_p_im, udata);
    func(x, y - h, &f_m_re, &f_m_im, udata);
    *J12 = (f_p_re - f_m_re) * inv_2h;
    *J22 = (f_p_im - f_m_im) * inv_2h;
}

/**
 * Single iteration of the Newton method
 *
 * Solves the linear system:
 *   J · Δz = -f(z)
 * where J is the 2x2 Jacobian matrix, then updates:
 *   z_new = z + λ·Δz  (λ = damping factor)
 *
 * @param state    Current state (updated in place)
 * @param func     Function f(z)
 * @param jacobian Analytic Jacobian (NULL to use finite differences)
 * @param params   Solver parameters
 * @param udata    User data
 * @return Status code (see NEWTON_* constants)
 */
static inline int newton_complex_step(
    newton_state_t* state,
    newton_func_t func,
    newton_jacobian_t jacobian,
    const newton_params_t* params,
    void* udata)
{
    state->iteration++;

    /* ─────────────────────────────────────────────────────────────────────
     *  1. Compute Jacobian
     * ───────────────────────────────────────────────────────────────────── */
    if (jacobian) {
        jacobian(state->x, state->y,
                 &state->J11, &state->J12,
                 &state->J21, &state->J22,
                 udata);
    } else {
        newton_jacobian_numeric(state->x, state->y,
                                &state->J11, &state->J12,
                                &state->J21, &state->J22,
                                func, params->diff_h, udata);
    }

    /* ─────────────────────────────────────────────────────────────────────
     *  2. Compute determinant and check for singularity
     * ───────────────────────────────────────────────────────────────────── */
    state->det_J = state->J11 * state->J22 - state->J12 * state->J21;

    if (REAL_FABS(state->det_J) < params->jacobian_tol) {
        state->status = NEWTON_SINGULAR;
        return NEWTON_SINGULAR;
    }

    /* ─────────────────────────────────────────────────────────────────────
     *  3. Solve 2x2 linear system: J · Δ = -f
     *     Using Cramer's rule:
     *       Δx = (-f_re·J22 + f_im·J12) / det
     *       Δy = (-f_im·J11 + f_re·J21) / det
     * ───────────────────────────────────────────────────────────────────── */
    REAL inv_det = REAL_LITERAL(1.0) / state->det_J;
    REAL raw_dx = (-state->f_re * state->J22 + state->f_im * state->J12) * inv_det;
    REAL raw_dy = (-state->f_im * state->J11 + state->f_re * state->J21) * inv_det;

    /* ─────────────────────────────────────────────────────────────────────
     *  4. Apply damping (line search)
     * ───────────────────────────────────────────────────────────────────── */
    REAL f_abs_old = state->f_abs;
    REAL x_new, y_new, f_re_new, f_im_new, f_abs_new;
    int backsteps = 0;
    REAL lambda = state->damping;

    if (params->use_damping) {
        /* Damped Newton with backtracking */
        while (backsteps < params->max_backsteps) {
            x_new = state->x + lambda * raw_dx;
            y_new = state->y + lambda * raw_dy;

            func(x_new, y_new, &f_re_new, &f_im_new, udata);
            f_abs_new = REAL_SQRT(f_re_new * f_re_new + f_im_new * f_im_new);

            /* Accept step if |f| decreased */
            if (f_abs_new <= f_abs_old) {
                /* Increase damping for next iteration (up to 1.0) */
                state->damping = lambda * params->damping_increase;
                if (state->damping > REAL_LITERAL(1.0)) {
                    state->damping = REAL_LITERAL(1.0);
                }
                break;
            }

            /* Reduce step size */
            lambda *= params->damping_decrease;
            backsteps++;
        }

        if (backsteps >= params->max_backsteps) {
            /* Failed to decrease |f|, take the step anyway but record divergence */
            x_new = state->x + lambda * raw_dx;
            y_new = state->y + lambda * raw_dy;
            func(x_new, y_new, &f_re_new, &f_im_new, udata);
            f_abs_new = REAL_SQRT(f_re_new * f_re_new + f_im_new * f_im_new);

            state->damping = params->damping_min;
        }
    } else {
        /* Pure Newton step */
        lambda = REAL_LITERAL(1.0);
        x_new = state->x + raw_dx;
        y_new = state->y + raw_dy;
        func(x_new, y_new, &f_re_new, &f_im_new, udata);
        f_abs_new = REAL_SQRT(f_re_new * f_re_new + f_im_new * f_im_new);
    }

    /* ─────────────────────────────────────────────────────────────────────
     *  5. Update state
     * ───────────────────────────────────────────────────────────────────── */
    state->delta_x = lambda * raw_dx;
    state->delta_y = lambda * raw_dy;
    state->delta_abs = REAL_SQRT(state->delta_x * state->delta_x +
                                  state->delta_y * state->delta_y);

    state->x = x_new;
    state->y = y_new;
    state->f_re = f_re_new;
    state->f_im = f_im_new;
    state->f_abs = f_abs_new;

    /* ─────────────────────────────────────────────────────────────────────
     *  6. Check convergence
     * ───────────────────────────────────────────────────────────────────── */
    if (state->delta_abs < params->delta_tol) {
        state->status = NEWTON_CONVERGED_DELTA;
        return NEWTON_CONVERGED_DELTA;
    }

    if (state->f_abs < params->f_tol) {
        state->status = NEWTON_CONVERGED_FABS;
        return NEWTON_CONVERGED_FABS;
    }

    state->status = NEWTON_CONTINUE;
    return NEWTON_CONTINUE;
}

/**
 * Full Newton solver: iterate until convergence or failure
 *
 * @param x0       Initial guess: real part (updated to solution)
 * @param y0       Initial guess: imaginary part (updated to solution)
 * @param func     Function f(z)
 * @param jacobian Analytic Jacobian (NULL for numeric)
 * @param params   Solver parameters
 * @param udata    User data
 * @param state_out Optional: final state (can be NULL)
 * @return Status code
 */
static inline int newton_complex_solve(
    REAL* x0, REAL* y0,
    newton_func_t func,
    newton_jacobian_t jacobian,
    const newton_params_t* params,
    void* udata,
    newton_state_t* state_out)
{
    newton_state_t state;
    newton_complex_init(&state, *x0, *y0, func, params, udata);

    /* Check if already converged */
    if (state.f_abs < params->f_tol) {
        if (state_out) *state_out = state;
        return NEWTON_CONVERGED_FABS;
    }

    int status;
    for (int i = 0; i < params->max_iterations; i++) {
        status = newton_complex_step(&state, func, jacobian, params, udata);

        if (status != NEWTON_CONTINUE) {
            *x0 = state.x;
            *y0 = state.y;
            if (state_out) *state_out = state;
            return status;
        }
    }

    /* Max iterations reached */
    *x0 = state.x;
    *y0 = state.y;
    state.status = NEWTON_MAX_ITER;
    if (state_out) *state_out = state;
    return NEWTON_MAX_ITER;
}

int solve_newton_root_classic(
    mendrive_scalar_t kz0, mendrive_scalar_t sz0,
    mendrive_scalar_t* kz_out, mendrive_scalar_t* sz_out,
    mendrive_scalar_t* f_abs_out,
    int max_iter);

#endif /* NEWTON_COMPLEX_CLASSIC_H */
