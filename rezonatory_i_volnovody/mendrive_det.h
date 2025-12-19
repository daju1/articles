#ifndef MENDDRIVE_DET_H
#define MENDDRIVE_DET_H

typedef struct {
    long double c;
    long double omega;
    long double a;

    // left conductor
    long double eps_l_xx, eps_l_yy, eps_l_zz;
    long double mu_l_xx, mu_l_yy, mu_l_zz;
    long double mu_l_yz, mu_l_zy;
    long double sigma_e_l_xx, sigma_e_l_yy, sigma_e_l_zz;
    long double sigma_m_l_xx, sigma_m_l_yy, sigma_m_l_zz;

    // right conductor
    long double eps_r_xx, eps_r_yy, eps_r_zz;
    long double mu_r_xx, mu_r_yy, mu_r_zz;
    long double mu_r_yz, mu_r_zy;
    long double sigma_e_r_xx, sigma_e_r_yy, sigma_e_r_zz;
    long double sigma_m_r_xx, sigma_m_r_yy, sigma_m_r_zz;
} mendrive_params_t;

void det_init( const mendrive_params_t* params );
void det_eval( long double kz, long double sz, long double *det_re, long double *det_im );

void det_div_diff_kz_eval(long double kz, long double sz,
    long double *div_re, long double *div_im);

void det_div_diff_sz_eval(long double kz, long double sz,
    long double *div_re, long double *div_im);

// Численные производные (по-умолчанию — центральная разность)
void det_derivatives(long double kz, long double sz,
                     long double* dfdkz_re, long double* dfdkz_im,
                     long double* dfdsz_re, long double* dfdsz_im);

// Шаг метода Ньютона
int newton_step(long double* xn_re, long double* xn_im,
                long double* f_abs_out,
                long double* delta_re_re_out,
                long double* delta_re_im_out,
                long double* delta_im_re_out,
                long double* delta_im_im_out);
#endif
