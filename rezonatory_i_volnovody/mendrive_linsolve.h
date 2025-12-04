// mendrive_linsolve.h
int solve_complex_least_squares(
    const double *A_re, const double *A_im, int n, int p,
    const double *b_re, const double *b_im,
    double *x_re, double *x_im,
    double *residual_norm);