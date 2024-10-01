long double sx(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);
long double sy(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);
long double sz(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);

long double vx(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);
long double vy(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);
long double vz(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);

long double wx(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);
long double wy(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);
long double wz(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);

int ccalc_sum_F_t(int N, long double t_i,
                    long double alpha0_l,
                    long double alpha0_r,
                    long double * Fx,
                    long double * Fy,
                    long double * Fz,
                    long double * F_alpha_l,
                    long double * F_alpha_r,
                    long double * sum_rlagerror_sqare,
                    _Bool to_log);
int ccalc_Maxwells_stress_tensor(long double X_a, long double Y_a, long double Z_a, long double t_i,
                                 int N,
                                 long double cos_nx, long double cos_ny, long double cos_nz,
                                 long double alpha0_l,
                                 long double alpha0_r,
                                 long double * T_xn,
                                 long double * T_yn,
                                 long double * T_zn,
                                 long double * N_x,
                                 long double * N_y,
                                 long double * N_z,
                                 long double * S_n,
                                 long double * E_n,
                                 long double * H_n,
                                 long double * A_n,
                                 long double * sum_rlagerror_sqare
                                 );

int spherical_ccalc_Maxwells_stress_tensor(
                                 long double r, long double theta, long double varphi, long double t,
                                 long double * Txn, long double * Tyn, long double * Tzn,
                                 long double * Nx, long double * Ny, long double * Nz,
                                 long double * Sn,
                                 long double * En,
                                 long double * Hn,
                                 long double * An,
                                 long double * sum_rlagerror_sqare);

int spherical_ccalc_Maxwells_stress_tensor_R_t(
                                 long double theta, long double varphi, long double t,
                                 long double * pTxn, long double * pTyn, long double * pTzn,
                                 long double * pNx, long double * pNy, long double * pNz,
                                 long double * pSn,
                                 long double * pEn,
                                 long double * pHn,
                                 long double * pAn,
                                 long double * sum_rlagerror_sqare);

void cset_sphere_R(long double R);
long double cget_sphere_R();

void cset_vc(long double vc);
long double cget_vc();
long double cget_omega(); 

long double cget_S();
long double cget_R_l();
long double cget_R_r();

long double cget_xc_l();
long double cget_xc_r();

long double cget_yc_l();
long double cget_yc_r();

long double cget_zc_l();
long double cget_zc_r();
