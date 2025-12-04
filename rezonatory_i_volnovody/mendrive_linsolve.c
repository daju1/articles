// mendrive_linsolve.c
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

int solve_complex_least_squares(
    const double *A_re, const double *A_im, int n, int p,
    const double *b_re, const double *b_im,
    double *x_re, double *x_im,
    double *residual_norm) {

    // Формируем действительную (2n × 2p) систему из комплексной:
    // [Re(A)  -Im(A)] [Re(x)]   [Re(b)]
    // [Im(A)   Re(A)] [Im(x)] = [Im(b)]
    gsl_matrix *X = gsl_matrix_alloc(2*n, 2*p);
    gsl_vector *y = gsl_vector_alloc(2*n);
    gsl_vector *c = gsl_vector_alloc(2*p);
    gsl_matrix *cov = gsl_matrix_alloc(2*p, 2*p);
    double chisq;

    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < p; ++j) {
            int idx = i * p + j;
            gsl_matrix_set(X, i,          j,         A_re[idx]);  //  Re(A)
            gsl_matrix_set(X, i,          j + p,    -A_im[idx]);  // -Im(A)
            gsl_matrix_set(X, i + n,      j,         A_im[idx]);  //  Im(A)
            gsl_matrix_set(X, i + n,      j + p,     A_re[idx]);  //  Re(A)
        }
        gsl_vector_set(y, i,          b_re[i]);  // Re(b)
        gsl_vector_set(y, i + n,      b_im[i]);  // Im(b)
    }

    gsl_multifit_linear_workspace *work = gsl_multifit_linear_alloc(2*n, 2*p);
    int status = gsl_multifit_linear(X, y, c, cov, &chisq, work);
    gsl_multifit_linear_free(work);

    if (status == GSL_SUCCESS) {
        for (int j = 0; j < p; ++j) {
            x_re[j] = gsl_vector_get(c, j);
            x_im[j] = gsl_vector_get(c, j + p);
        }
        *residual_norm = sqrt(chisq);
    }

    gsl_matrix_free(X); gsl_vector_free(y); gsl_vector_free(c); gsl_matrix_free(cov);
    return status;
}