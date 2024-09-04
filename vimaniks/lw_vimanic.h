long double sx(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);

long double sy(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);

long double sz(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);

long double vx(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);

long double vy(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);

long double vz(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);

long double wx(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);
long double wy(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);
long double wz(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);

long double ccalc_sum_Fy_t(int N, long double t_i,
                  long double alpha0_l,
                  long double alpha0_r, _Bool to_log);
long double ccalc_Maxwells_stress_tensor(long double X_a, long double Y_a, long double Z_a, long double t_i,
                                 int N,
                                long double cos_nx, long double cos_ny, long double cos_nz,
                                long double alpha0_l,
                                long double alpha0_r
                               );

long double spherical_ccalc_Maxwells_stress_tensor(
    long double r, long double theta, long double varphi, long double t);

long double spherical_ccalc_Maxwells_stress_tensor_R_t(
    long double theta, long double varphi, long double t);

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
