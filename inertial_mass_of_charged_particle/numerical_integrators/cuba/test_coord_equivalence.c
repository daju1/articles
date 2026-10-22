#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

typedef enum { M_Z_Z = 0, M_Z_RHO = 1, M_Z_PHI = 2 } Mode;

/* (1) прямое отображение: локальные сферические -> лабораторные декартовы (строки из ваших веток #ifdef) */
static void local_to_lab(Mode m, double rho0, double r, double th, double ps,
                         double *x, double *y, double *z)
{
    switch (m) {
    case M_Z_Z:   *x = rho0 + r*sin(th)*cos(ps); *y = r*sin(th)*sin(ps); *z = r*cos(th);        break;
    case M_Z_RHO: *x = rho0 + r*cos(th);         *y = r*sin(th)*sin(ps); *z = r*sin(th)*cos(ps); break;
    case M_Z_PHI: *x = rho0 + r*sin(th)*cos(ps); *y = r*cos(th);         *z = r*sin(th)*sin(ps); break;
    }
}

/* (2) обратное отображение: смещение от центра (dx,dy,dz) -> (r,theta,psi) в заданной ориентации */
static void lab_to_local(Mode m, double dx, double dy, double dz,
                         double *r, double *th, double *ps)
{
    *r = sqrt(dx*dx + dy*dy + dz*dz);
    double ax, ay, az;
    switch (m) {
    case M_Z_Z:   az = dz; ax = dx; ay = dy; break;
    case M_Z_RHO: az = dx; ax = dz; ay = dy; break;
    default:      az = dy; ax = dx; ay = dz; break;
    }
    *th = acos(az / *r);
    *ps = atan2(ay, ax);
}

typedef struct {
    double x_a,y_a,z_a, x_q,y_q,z_q;
    double rho_a,phi_a,rho_q,phi_q;
    double R_rho,R_phi,R_z,R, v_phi,a_rho;
} Block;

/* (3) первый блок compute_electric_field() после ветвления по дефайну */
static Block first_block(Mode m, double rho0, double omega, double c, int lorentz,
                         double ra, double ta, double pa, double rq, double tq, double pq)
{
    Block b;
    local_to_lab(m, rho0, ra, ta, pa, &b.x_a, &b.y_a, &b.z_a);
    local_to_lab(m, rho0, rq, tq, pq, &b.x_q, &b.y_q, &b.z_q);
    b.rho_a = sqrt(b.x_a*b.x_a + b.y_a*b.y_a); b.phi_a = atan2(b.y_a, b.x_a); if (b.phi_a < 0) b.phi_a += 2*M_PI;
    b.rho_q = sqrt(b.x_q*b.x_q + b.y_q*b.y_q); b.phi_q = atan2(b.y_q, b.x_q); if (b.phi_q < 0) b.phi_q += 2*M_PI;
    b.v_phi = b.rho_q*omega;  b.a_rho = -b.rho_q*omega*omega;
    if (lorentz) {
        double gamma = 1.0/sqrt(1.0 - b.v_phi*b.v_phi/(c*c));
        b.y_q /= gamma;  b.y_a /= gamma;
        b.rho_a = sqrt(b.x_a*b.x_a + b.y_a*b.y_a); b.phi_a = atan2(b.y_a, b.x_a); if (b.phi_a < 0) b.phi_a += 2*M_PI;
        b.rho_q = sqrt(b.x_q*b.x_q + b.y_q*b.y_q); b.phi_q = atan2(b.y_q, b.x_q); if (b.phi_q < 0) b.phi_q += 2*M_PI;
        b.v_phi = b.rho_q*omega;  b.a_rho = -b.rho_q*omega*omega;
    }
    b.R_rho = b.rho_a*cos(b.phi_a - b.phi_q) - b.rho_q;
    b.R_phi = b.rho_a*sin(b.phi_a - b.phi_q);
    b.R_z   = b.z_a - b.z_q;
    b.R     = sqrt(b.R_rho*b.R_rho + b.R_phi*b.R_phi + b.R_z*b.R_z);
    return b;
}

/* get_source_position(): копия вашей (без Лоренца). fixed_version=0 - как сейчас, 1 - x=rho*cos, y=rho*sin */
static void get_source_position(Mode m, double tp, double rq, double tq, double pq,
                                double rho0, double omega, int fixed_version,
                                double *x, double *y, double *z)
{
    double x0,y0,z0;  local_to_lab(m, rho0, rq, tq, pq, &x0,&y0,&z0);
    double rho = sqrt(x0*x0 + y0*y0);
    double phi = atan2(y0, x0); if (phi < 0) phi += 2*M_PI;
    phi += omega*tp;
    if (fixed_version) { *x = rho*cos(phi); *y = rho*sin(phi); }
    else               { *x = rho*sin(phi); *y = rho*cos(phi); }
    *z = z0;
}

static double dphi(double a, double b) { double d = fmod(a-b, 2*M_PI); if (d > M_PI) d -= 2*M_PI; if (d < -M_PI) d += 2*M_PI; return d; }
static double urand(void) { return rand()/(RAND_MAX+1.0); }
static void rand_ball(double R0, double *dx, double *dy, double *dz) {
    double r = R0*cbrt(urand()), t = acos(1-2*urand()), p = 2*M_PI*urand();
    *dx = r*sin(t)*cos(p); *dy = r*sin(t)*sin(p); *dz = r*cos(t);
}

int main(void)
{
    const double R0 = 1.0, c = 1.0, tol = 1e-11;
    double rho0s[3] = { 0.92, 2.0, 100.0 };
    srand(12345);

    printf("=== T1+T2: одна физическая пара точек в Z_Z / Z_RHO / Z_PHI (эталон: Z_Z) ===\n");
    for (int ir = 0; ir < 3; ir++) for (int lor = 0; lor < 2; lor++) {
        double rho0 = rho0s[ir], omega = sqrt(5*R0*c*c/(3*2*4*rho0*rho0*rho0));
        double worst[3] = {0,0,0}; int bad = 0;
        for (int k = 0; k < 20000; k++) {
            double ax,ay,az,qx,qy,qz, ra,ta,pa,rq,tq,pq;
            rand_ball(R0,&ax,&ay,&az);  rand_ball(R0,&qx,&qy,&qz);
            lab_to_local(M_Z_Z, ax,ay,az,&ra,&ta,&pa);  lab_to_local(M_Z_Z, qx,qy,qz,&rq,&tq,&pq);
            Block ref = first_block(M_Z_Z, rho0, omega, c, lor, ra,ta,pa, rq,tq,pq);
            for (int m = 0; m < 3; m++) {
                lab_to_local((Mode)m, ax,ay,az,&ra,&ta,&pa);  lab_to_local((Mode)m, qx,qy,qz,&rq,&tq,&pq);
                Block b = first_block((Mode)m, rho0, omega, c, lor, ra,ta,pa, rq,tq,pq);
                double d = 0;
                d = fmax(d, fabs(b.rho_a-ref.rho_a)); d = fmax(d, fabs(dphi(b.phi_a,ref.phi_a))*fmin(1,ref.rho_a));
                d = fmax(d, fabs(b.z_a-ref.z_a));
                d = fmax(d, fabs(b.rho_q-ref.rho_q)); d = fmax(d, fabs(dphi(b.phi_q,ref.phi_q))*fmin(1,ref.rho_q));
                d = fmax(d, fabs(b.z_q-ref.z_q));
                d = fmax(d, fabs(b.R_rho-ref.R_rho)); d = fmax(d, fabs(b.R_phi-ref.R_phi));
                d = fmax(d, fabs(b.R_z-ref.R_z));     d = fmax(d, fabs(b.R-ref.R));
                d = fmax(d, fabs(b.v_phi-ref.v_phi)); d = fmax(d, fabs(b.a_rho-ref.a_rho)/fabs(ref.a_rho));
                if (d > worst[m]) worst[m] = d;
                if (d > tol) bad++;
            }
        }
        printf(" rho0=%6.2f lorentz=%d : max|diff|  Z_Z=%.1e  Z_RHO=%.1e  Z_PHI=%.1e   %s\n",
               rho0, lor, worst[0], worst[1], worst[2], bad ? "FAIL" : "OK");
    }

    printf("\n=== T3: get_source_position(t'=0) против первого блока ===\n");
    {
        double rho0 = 2.0, omega = 0.01, dx,dy,dz, r,t,p, xs,ys,zs;
        rand_ball(R0,&dx,&dy,&dz);  lab_to_local(M_Z_RHO,dx,dy,dz,&r,&t,&p);
        Block b = first_block(M_Z_RHO, rho0, omega, c, 0, r,t,p, r,t,p);
        printf(" блок:                           (%+.6f, %+.6f, %+.6f)\n", b.x_q,b.y_q,b.z_q);
        get_source_position(M_Z_RHO, 0.0, r,t,p, rho0, omega, 0, &xs,&ys,&zs);
        printf(" get_source_position (как есть): (%+.6f, %+.6f, %+.6f)\n", xs,ys,zs);
        get_source_position(M_Z_RHO, 0.0, r,t,p, rho0, omega, 1, &xs,&ys,&zs);
        printf(" get_source_position (исправл.): (%+.6f, %+.6f, %+.6f)\n", xs,ys,zs);
    }

    printf("\n=== T4: направление вращения центра за dt'=1 (omega=+0.01) ===\n");
    {
        double rho0 = 2.0, omega = 0.01, x0,y0,z0, x1,y1,z1;
        for (int fixed = 0; fixed < 2; fixed++) {
            get_source_position(M_Z_RHO, 0.0, 0,0,0, rho0,omega, fixed, &x0,&y0,&z0);
            get_source_position(M_Z_RHO, 1.0, 0,0,0, rho0,omega, fixed, &x1,&y1,&z1);
            printf(" %-9s: сдвиг угла %+.5f рад (ожидается %+.5f)\n",
                   fixed ? "исправл." : "как есть", dphi(atan2(y1,x1), atan2(y0,x0)), omega);
        }
    }

    printf("\n=== T5: инвариантность g=R_rho^2*a_rho/R^3 к (psi_a,psi_q)->(+d,+d); 0 = ось z есть ось симметрии ===\n");
    for (int ir = 0; ir < 3; ir++) {
        double rho0 = rho0s[ir], omega = sqrt(5*R0*c*c/(3*2*4*rho0*rho0*rho0));
        double sd[3] = {0,0,0}, sg[3] = {0,0,0};
        for (int k = 0; k < 200000; k++) {
            double ra = R0*cbrt(urand()), rq = R0*cbrt(urand()), ta = acos(1-2*urand()), tq = acos(1-2*urand());
            double pa = 2*M_PI*urand(), pq = 2*M_PI*urand(), d = 2*M_PI*urand();
            for (int m = 0; m < 3; m++) {
                Block b1 = first_block((Mode)m, rho0, omega, c, 0, ra,ta,pa,   rq,tq,pq);
                Block b2 = first_block((Mode)m, rho0, omega, c, 0, ra,ta,pa+d, rq,tq,pq+d);
                double g1 = b1.R_rho*b1.R_rho*b1.a_rho/pow(b1.R,3), g2 = b2.R_rho*b2.R_rho*b2.a_rho/pow(b2.R,3);
                sd[m] += (g1-g2)*(g1-g2);  sg[m] += g1*g1;
            }
        }
        printf(" rho0=%6.2f : Z_Z=%.3f   Z_RHO=%.3f   Z_PHI=%.3f\n", rho0, sqrt(sd[0]/sg[0]), sqrt(sd[1]/sg[1]), sqrt(sd[2]/sg[2]));
    }
    return 0;
}
