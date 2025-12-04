// mendrive_root.h
#ifndef MENDDRIVE_ROOT_H
#define MENDDRIVE_ROOT_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

// ——— параметры ———
void set_params(double c, double omega, double a,
                double eps_l_xx, double mu_l_yy, double sigma_e_l_xx, double sigma_m_l_yy,
                double eps_r_xx, double mu_r_yy, double sigma_e_r_xx, double sigma_m_r_yy);

// ——— система ———
int det_system(const gsl_vector *x, void *params, gsl_vector *f);

// ——— решение ———
int solve_newton_root(double *kz_guess, double *sz_guess,
                      double *kz_sol, double *sz_sol,
                      double epsabs, int max_iter);

#endif