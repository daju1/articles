// mendrive_root.h
#ifndef MENDDRIVE_ROOT_H
#define MENDDRIVE_ROOT_H

#include <gsl/gsl_vector.h>
#include <gsl/gsl_multiroots.h>

// ——— система ———
int det_system(const gsl_vector *x, void *params, gsl_vector *f);

// ——— решение ———
int solve_newton_root(double *kz_guess, double *sz_guess,
                      double *kz_sol, double *sz_sol,
                      double epsabs, int max_iter);

#endif