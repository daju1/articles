// mendrive_root.c
#include <math.h>
#include <complex.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "mendrive_root.h"
#include "mendrive_det.h"  // ← содержит det_eval()

static mendrive_params_t g_params;

int det_system(const gsl_vector *x, void *params, gsl_vector *f) {
    long double kz = (long double) gsl_vector_get(x, 0);
    long double sz = (long double) gsl_vector_get(x, 1);

    mendrive_scalar_t re, im;
    det_eval(kz, sz, &re, &im);

    gsl_vector_set(f, 0, (double)re);  // f₁ = Re(det)
    gsl_vector_set(f, 1, (double)im);  // f₂ = Im(det)
    return GSL_SUCCESS;
}

int solve_newton_root(double *kz_guess, double *sz_guess,
                      double *kz_sol, double *sz_sol,
                      double epsabs, int max_iter) {
    const gsl_multiroot_fsolver_type *T = gsl_multiroot_fsolver_hybrids;
    gsl_multiroot_fsolver *s = gsl_multiroot_fsolver_alloc(T, 2);
    gsl_vector *x = gsl_vector_alloc(2);
    gsl_multiroot_function F = {&det_system, 2, NULL};

    gsl_vector_set(x, 0, *kz_guess);
    gsl_vector_set(x, 1, *sz_guess);

    gsl_multiroot_fsolver_set(s, &F, x);

    int status, iter = 0;
    do {
        iter++;
        status = gsl_multiroot_fsolver_iterate(s);
        if (status) break;

        status = gsl_multiroot_test_residual(s->f, epsabs);
    } while (status == GSL_CONTINUE && iter < max_iter);

    *kz_sol = gsl_vector_get(s->x, 0);
    *sz_sol = gsl_vector_get(s->x, 1);

    gsl_multiroot_fsolver_free(s);
    gsl_vector_free(x);
    return status;
}