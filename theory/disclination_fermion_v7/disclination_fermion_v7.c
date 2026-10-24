// ===================== disclination_fermion_v7.c =====================
// Фермион = диспирационная петля с ВЕКТОРНЫМ твистом b(θ) в 3D.
// Плоскость твиста + winding + киральность; шесть каналов взаимодействия.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "disclination_fermion_v7.h"

static void vadd(const double a[3],const double b[3],double o[3]){o[0]=a[0]+b[0];o[1]=a[1]+b[1];o[2]=a[2]+b[2];}
static void vsub(const double a[3],const double b[3],double o[3]){o[0]=a[0]-b[0];o[1]=a[1]-b[1];o[2]=a[2]-b[2];}
static void vscale(const double a[3],double s,double o[3]){o[0]=a[0]*s;o[1]=a[1]*s;o[2]=a[2]*s;}
static double vdot(const double a[3],const double b[3]){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
static double vnorm(const double a[3]){return sqrt(vdot(a,a));}
static void vcross(const double a[3],const double b[3],double o[3]){
    o[0]=a[1]*b[2]-a[2]*b[1];o[1]=a[2]*b[0]-a[0]*b[2];o[2]=a[0]*b[1]-a[1]*b[0];}
static void vnormalize(double a[3]){double n=vnorm(a); if(n>1e-30){a[0]/=n;a[1]/=n;a[2]/=n;}}

/* Ортотройка в точке θ кольца (кольцо в плоскости xy). */
static void frame(double theta,double t[3],double r[3],double z[3]){
    t[0]=-sin(theta); t[1]=cos(theta); t[2]=0.0;
    r[0]= cos(theta); r[1]=sin(theta); r[2]=0.0;
    z[0]=0.0; z[1]=0.0; z[2]=1.0;
}
/* Нормаль плоскости твиста (ось Франка, до умножения на n/2). */
static void twist_normal(const Disclination *D,double theta,double nn[3]){
    double t[3],r[3],z[3]; frame(theta,t,r,z);
    switch(D->plane){
        case PLANE_T_R: nn[0]=z[0];nn[1]=z[1];nn[2]=z[2]; break;        /* ⊥ плоскости кольца = ẑ */
        case PLANE_R_Z: nn[0]=t[0];nn[1]=t[1];nn[2]=t[2]; break;        /* нормаль = t̂ */
        case PLANE_T_Z: nn[0]=r[0];nn[1]=r[1];nn[2]=r[2]; break;        /* нормаль = r̂ */
        case PLANE_TILTED: {                                            /* наклон от ẑ к r̂ на β */
            double b=D->normal_beta;
            for(int i=0;i<3;i++) nn[i]=cos(b)*z[i]+sin(b)*r[i];
        } break;
    }
    vnormalize(nn);
}
/* Два орта в плоскости твиста (⊥ нормали). */
static void twist_basis(const Disclination *D,double theta,double e1[3],double e2[3]){
    double t[3],r[3],z[3],nn[3]; frame(theta,t,r,z); twist_normal(D,theta,nn);
    /* e1 = проекция r̂ на плоскость (если r почти || nn, берём t̂) */
    double ref[3]; if(fabs(vdot(nn,r))<0.9){ref[0]=r[0];ref[1]=r[1];ref[2]=r[2];}
                   else {ref[0]=t[0];ref[1]=t[1];ref[2]=t[2];}
    double proj=vdot(ref,nn),tmp[3]; vscale(nn,proj,tmp); vsub(ref,tmp,e1); vnormalize(e1);
    vcross(nn,e1,e2);
}

/* Вектор Бюргерса b(θ): вращается в плоскости твиста с фазой chirality·n·θ/2. */
void burgers_vector(const Disclination *D,double theta,double b_out[3]){
    double e1[3],e2[3]; twist_basis(D,theta,e1,e2);
    double ph=D->chirality*D->winding*theta/2.0, c=cos(ph),s=sin(ph);
    for(int i=0;i<3;i++) b_out[i]=D->b0*(c*e1[i]+s*e2[i]);
}
/* Вектор Франка Ω(θ) = нормаль × (n/2) (скорость накопления поворота b). */
void frank_vector(const Disclination *D,double theta,double O_out[3]){
    double nn[3]; twist_normal(D,theta,nn);
    vscale(nn,D->chirality*D->winding/2.0,O_out);
}

/* Магнитный момент μ: усреднён по 4π от дисклинац. части (Франк Ω).
   ЧИСЛЕННО: μ_z≠0 только когда Ω||ẑ (плоскость T_R). */
void magnetic_moment(const Disclination *D,double mu_out[3]){
    int N=D->n_elem*50; if(N<4000)N=4000;
    double dth=4.0*M_PI/N; mu_out[0]=mu_out[1]=mu_out[2]=0.0;
    for(int k=0;k<N;k++){ double theta=D->phase0+dth*k, O[3];
        frank_vector(D,theta,O); vadd(mu_out,O,mu_out); /* накопление */
    }
    for(int i=0;i<3;i++) mu_out[i]*=dth/(4.0*M_PI);
}
/* Момент количества движения J: орбитальный (прецессия) + спиновый (Франк). */
void angular_momentum(const Disclination *D,double c,double J_out[3]){
    (void)c; int N=D->n_elem*20; if(N<2000)N=2000;
    double dth=4.0*M_PI/N; J_out[0]=J_out[1]=J_out[2]=0.0;
    for(int k=0;k<N;k++){ double theta=D->phase0+dth*k;
        double t[3],r[3],z[3]; frame(theta,t,r,z);
        double pos[3]; vscale(r,D->a,pos);
        double vel[3]; vscale(t,D->a*D->omega,vel);
        double b[3]; burgers_vector(D,theta,b);
        double dens=vdot(b,b);                  /* масса-плотность ~ b² */
        double rv[3]; vcross(pos,vel,rv);
        for(int i=0;i<3;i++) J_out[i]+=dens*rv[i]*dth;
    }
    for(int i=0;i<3;i++) J_out[i]/=(4.0*M_PI);
}
/* Гиромагнитное: возвращает μ_z, S_z (спиновую часть) и отношение. */
double gyromagnetic_g(const Disclination *D,double c,double *mu_z,double *S_z){
    double mu[3],J[3]; magnetic_moment(D,mu); angular_momentum(D,c,J);
    if(mu_z)*mu_z=mu[2]; if(S_z)*S_z=spin_value(D);
    double s=spin_value(D); return (s>1e-12)? fabs(mu[2])/s : 0.0;
}
double spin_value(const Disclination *D){ return fabs(D->winding)/2.0; }
int twist_period_2pi_factor(const Disclination *D){ return (D->winding%2==0)?1:-1; }

/* ---- поле элемента контура (дислокац. T(3), запаздывающее) ---- */
static void elem_kin(const Disclination *D,int k,double t,double pos[3],double vel[3],double acc[3]){
    double theta=D->phase0+2.0*M_PI*k/D->n_elem+D->omega*t;
    double tt[3],r[3],z[3]; frame(theta,tt,r,z);
    vscale(r,D->a,pos); vadd(D->center,pos,pos);
    vscale(tt,D->a*D->omega,vel);
    vscale(r,-D->a*D->omega*D->omega,acc);
}
static double retarded_time(const double fp[3],double tf,const Disclination *D,int k,double c){
    double pos[3],vel[3],acc[3],d[3]; elem_kin(D,k,tf,pos,vel,acc);
    vsub(fp,pos,d); double R=vnorm(d),tp=tf-R/c;
    for(int i=0;i<25;i++){ elem_kin(D,k,tp,pos,vel,acc); vsub(fp,pos,d);
        double Rn=vnorm(d),f=tp+Rn/c-tf,vr=vdot(d,vel)/(Rn+1e-30),df=1.0-vr/c;
        if(fabs(df)<1e-13)break; double dt=-f/df; tp+=dt; if(fabs(dt)<1e-15)break; }
    return tp;
}
typedef struct { double central[3]; double vortex[3]; } Field;
static Field elem_field(const double fp[3],const Disclination *D,int k,double tf,double c){
    Field out; memset(&out,0,sizeof(out));
    double tp=retarded_time(fp,tf,D,k,c),pos[3],vel[3],acc[3];
    elem_kin(D,k,tp,pos,vel,acc);
    double Rv[3]; vsub(fp,pos,Rv); double R=vnorm(Rv); if(R<1e-9)return out;
    double n[3]; vscale(Rv,1.0/R,n);
    double vn=vdot(vel,n),Rstar=R*(1.0-vn/c); if(Rstar<1e-9)Rstar=1e-9;
    double v2=vdot(vel,vel)/(c*c),Ra=vdot(Rv,acc)/(c*c);
    double common=1.0/(Rstar*Rstar*Rstar),br=(1.0+Ra-v2),tc[3],ta[3];
    for(int i=0;i<3;i++) tc[i]=br*Rv[i]-br*R*vel[i]/c;
    for(int i=0;i<3;i++) ta[i]=acc[i]*Rstar*R/(c*c);
    for(int i=0;i<3;i++) out.central[i]=common*(tc[i]-ta[i]);
    vcross(n,out.central,out.vortex);
    return out;
}
/* Мгновенная сила по каналу ch. */
static void pair_force_inst(const Disclination *F1,const Disclination *F2,
                            double t,double c,Channel ch,double out[3]){
    out[0]=out[1]=out[2]=0.0;
    for(int i=0;i<F1->n_elem;i++){
        double fpos[3],fvel[3],facc[3]; elem_kin(F1,i,t,fpos,fvel,facc);
        double th1=F1->phase0+2.0*M_PI*i/F1->n_elem+F1->omega*t, b1v[3];
        burgers_vector(F1,th1,b1v); double b1=vnorm(b1v);
        for(int j=0;j<F2->n_elem;j++){
            double th2=F2->phase0+2.0*M_PI*j/F2->n_elem+F2->omega*t, b2v[3];
            burgers_vector(F2,th2,b2v); double b2=vnorm(b2v);
            Field fld=elem_field(fpos,F2,j,t,c);
            if(ch==CH_COULOMB){
                for(int k=0;k<3;k++) out[k]+=b1*b2*fld.central[k];
            } else if(ch==CH_GRAVITY){
                double W=vnorm(facc)/(c*c);
                for(int k=0;k<3;k++) out[k]-=(b1*b1)*(b2*b2)*W*fld.central[k];
            } else if(ch==CH_EXCHANGE){
                double vxB[3]; vcross(fvel,fld.vortex,vxB);
                for(int k=0;k<3;k++) out[k]+=b1*b1*b2*b2*vxB[k]/c;
            } else if(ch==CH_MAGNETIC){
                /* МАГНИТНОЕ (Осипов/Кочетов arXiv:1001.4450): взаимодействие
                   потоков дисклинаций Ааронов-Бома. Поток Φ_i=2π·ν_i, ν=winding/2
                   (индекс Франка). Сила диполь-дипольная ∝ Φ1Φ2/R³ вдоль линии
                   центров. НЕ Франк×vortex (та усреднялась в ноль). */
                double nu1=F1->winding/2.0, nu2=F2->winding/2.0;
                double Phi1=2.0*M_PI*nu1, Phi2=2.0*M_PI*nu2;
                double Rv[3]; vsub(F1->center,F2->center,Rv); double Rm=vnorm(Rv);
                if(Rm>1e-9){ double coef=Phi1*Phi2/(Rm*Rm*Rm)/(F1->n_elem*F2->n_elem);
                    for(int k=0;k<3;k++) out[k]+=coef*Rv[k]/Rm; }
            } else if(ch==CH_WEAK){
                /* слабое: ЗАВИСИТ ОТ ЗНАКА киральности (P-нарушение) */
                double chir=(double)(F1->chirality*F2->chirality);
                double proj=(1.0+chir)/2.0;   /* 1 если одинаковая кирал., 0 если разная */
                for(int k=0;k<3;k++) out[k]+=proj*b1*b2*fld.vortex[k]/c;
            } else { /* CH_STRONG: заготовка — линейный потенциал (лента SF) */
                double Rv[3]; vsub(fpos,F2->center,Rv); double Rmag=vnorm(Rv);
                double frac=(F1->b0<0.99)?1.0:0.0;  /* активен только для дробного b (кварк) */
                if(Rmag>1e-9) for(int k=0;k<3;k++) out[k]-=frac*b1*b2*Rv[k]/Rmag; /* const натяжение */
            }
        }
    }
}
/* Сила, усреднённая по периоду прецессии. */
void pair_force_avg(const Disclination *F1,const Disclination *F2,
                    double c,Channel ch,int nT,double out[3]){
    out[0]=out[1]=out[2]=0.0; double T=2.0*M_PI/F1->omega;
    for(int it=0;it<nT;it++){ double t=T*it/nT,f[3];
        pair_force_inst(F1,F2,t,c,ch,f); out[0]+=f[0];out[1]+=f[1];out[2]+=f[2]; }
    out[0]/=nT; out[1]/=nT; out[2]/=nT;
}
