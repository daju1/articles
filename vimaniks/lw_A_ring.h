void cset_vc(long double vc);
long double cget_vc();

long double cget_omega();

long double sx(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha);

long double sy(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha);

long double sz(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha);

long double vx(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha);

long double vy(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha);

long double vz(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha);

long double wx(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha);

long double wy(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha);

long double wz(long double t, long double xc, long double yc, long double zc,
               long double R, long double omega, long double alpha);

int ccalc_A_ring(long double X_a, long double Y_a, long double Z_a, long double t,
                 charge q, long double alpha0, long double *, coordinate * rlagerror);

long double cget_R_ring();
void cset_R_ring(long double r);

long double cget_xc_ring();
long double cget_yc_ring();
long double cget_zc_ring();

void cset_xc_ring(long double x);
void cset_yc_ring(long double y);
void cset_zc_ring(long double z);
