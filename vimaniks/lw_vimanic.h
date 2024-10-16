long double xc_l();
long double xc_r();

long double yc_l();
long double yc_r();

long double zc_l();
long double zc_r();


long double sx(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);
long double sy(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);
long double sz(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);

long double vx(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);
long double vy(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);
long double vz(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);

long double wx(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);
long double wy(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);
long double wz(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);

long double dot_wx(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);
long double dot_wy(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);
long double dot_wz(long double t, long double xc, long double yc, long double zc, long double R, long double omega, long double alpha);


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


//#define USE_LEFT_CHARGE_ONLY
//#define USE_LEFT_CENTERED_CHARGE_ONLY
//#define TEST_WITH_X_DIRECTED_DIPOLE
