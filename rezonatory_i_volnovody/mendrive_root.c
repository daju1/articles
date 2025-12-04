// mendrive_root.c
#include <math.h>
#include <complex.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>
#include "mendrive_root.h"
#include "mendrive_det.h"  // ← содержит det_eval()

static mendrive_params_t g_params;

void set_params(double c, double omega, double a,
                double eps_l_xx, double mu_l_yy, double sigma_e_l_xx, double sigma_m_l_yy,
                double eps_r_xx, double mu_r_yy, double sigma_e_r_xx, double sigma_m_r_yy) {
    g_params.c = c;
    g_params.omega = omega;
    g_params.a = a;
    g_params.eps_l_xx = eps_l_xx; g_params.mu_l_yy = mu_l_yy;
    g_params.sigma_e_l_xx = sigma_e_l_xx; g_params.sigma_m_l_yy = sigma_m_l_yy;
    g_params.eps_r_xx = eps_r_xx; g_params.mu_r_yy = mu_r_yy;
    g_params.sigma_e_r_xx = sigma_e_r_xx; g_params.sigma_m_r_yy = sigma_m_r_yy;
    det_init(&g_params);
}

int det_system(const gsl_vector *x, void *params, gsl_vector *f) {
    double kz = gsl_vector_get(x, 0);
    double sz = gsl_vector_get(x, 1);

    double re, im;
    det_eval(kz, sz, &re, &im);

    gsl_vector_set(f, 0, re);  // f₁ = Re(det)
    gsl_vector_set(f, 1, im);  // f₂ = Im(det)
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