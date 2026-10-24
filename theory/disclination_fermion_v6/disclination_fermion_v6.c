// ===================== disclination_fermion_v6.c =====================
// Фермион = одна диспирационная петля с твистом Бюргерса b(θ)=cos(θ/2).
// Три канала (кулон/обмен/гравитация) из одного твист-контура.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "disclination_fermion_v6.h"

static void vadd(const double a[3],const double b[3],double o[3]){o[0]=a[0]+b[0];o[1]=a[1]+b[1];o[2]=a[2]+b[2];}
static void vsub(const double a[3],const double b[3],double o[3]){o[0]=a[0]-b[0];o[1]=a[1]-b[1];o[2]=a[2]-b[2];}
static void vscale(const double a[3],double s,double o[3]){o[0]=a[0]*s;o[1]=a[1]*s;o[2]=a[2]*s;}
static double vdot(const double a[3],const double b[3]){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
static double vnorm(const double a[3]){return sqrt(vdot(a,a));}
static void vcross(const double a[3],const double b[3],double o[3]){
    o[0]=a[1]*b[2]-a[2]*b[1];o[1]=a[2]*b[0]-a[0]*b[2];o[2]=a[0]*b[1]-a[1]*b[0];}

/* Угол k-го элемента контура в момент t. */
static double elem_theta(const Disclination *D,int k,double t){
    return D->phase0 + 2.0*M_PI*k/D->n_elem + D->omega*t;
}
/* ТВИСТ: локальная величина Бюргерса b(θ)=b0·cos(θ/2) — период 4π.
   Это ключ v6: знак меняется за один обход 2π, возврат за 4π (спинорность). */
double burgers_twist(const Disclination *D,int k,double t){
    double th=elem_theta(D,k,t);
    return D->b0*cos(0.5*th);
}

/* Кинематика k-го элемента (положение/скорость/ускорение на окружности). */
static void elem_kin(const Disclination *D,int k,double t,double pos[3],double vel[3],double acc[3]){
    double th=elem_theta(D,k,t), c=cos(th), s=sin(th), rad[3],tmp[3];
    vscale(D->e1,c,rad); vscale(D->e2,s,tmp); vadd(rad,tmp,rad);
    vscale(rad,D->a,tmp); vadd(D->center,tmp,pos);
    vscale(D->e1,-s,vel); vscale(D->e2,c,tmp); vadd(vel,tmp,vel);
    vscale(vel,D->a*D->omega,vel);
    vscale(rad,-D->a*D->omega*D->omega,acc);
}

/* Мультипольные моменты твиста: monopole=∫b dθ/2π, |dipole|=|∫b e^{iθ}dθ|. */
void twist_multipoles(const Disclination *D,double t,double *monopole,double *dipole){
    int N=D->n_elem; double m=0,dx=0,dy=0;
    for(int k=0;k<N;k++){ double th=elem_theta(D,k,t), b=burgers_twist(D,k,t);
        m+=b; dx+=b*cos(th); dy+=b*sin(th); }
    double dth=2.0*M_PI/N;
    if(monopole)*monopole=m*dth/(2.0*M_PI);
    if(dipole)*dipole=sqrt(dx*dx+dy*dy)*dth;
}

/* Спин/орбита: твист поворачивает поле на π за обход (2π), орбита на 2π.
   S/L_орб = π/2π = 1/2 → спин-½. */
double disclination_spin_ratio(const Disclination *D){ (void)D; return 0.5; }
/* g-фактор: заряд проходится дважды за цикл спина (4π=2·2π) → g_S=2. */
double disclination_g_factor(const Disclination *D){ (void)D; return 2.0; }

/* Поле элемента контура (дислокац. T(3), запаздывающее, central+vortex). */
typedef struct { double central[3]; double vortex[3]; } Field;
static double retarded_time(const double fp[3],double tf,const Disclination *D,int k,double c){
    double pos[3],vel[3],acc[3]; elem_kin(D,k,tf,pos,vel,acc);
    double d[3]; vsub(fp,pos,d); double R=vnorm(d),tp=tf-R/c;
    for(int i=0;i<25;i++){ elem_kin(D,k,tp,pos,vel,acc); vsub(fp,pos,d);
        double Rn=vnorm(d),f=tp+Rn/c-tf,vr=vdot(d,vel)/(Rn+1e-30),df=1.0-vr/c;
        if(fabs(df)<1e-13)break; double dt=-f/df; tp+=dt; if(fabs(dt)<1e-15)break; }
    return tp;
}
static Field elem_field(const double fp[3],const Disclination *D,int k,double tf,double c){
    Field out; memset(&out,0,sizeof(out));
    double tp=retarded_time(fp,tf,D,k,c);
    double pos[3],vel[3],acc[3]; elem_kin(D,k,tp,pos,vel,acc);
    double Rv[3]; vsub(fp,pos,Rv); double R=vnorm(Rv); if(R<1e-9)return out;
    double n[3]; vscale(Rv,1.0/R,n);
    double vn=vdot(vel,n),Rstar=R*(1.0-vn/c); if(Rstar<1e-9)Rstar=1e-9;
    double v2=vdot(vel,vel)/(c*c),Ra=vdot(Rv,acc)/(c*c);
    double common=1.0/(Rstar*Rstar*Rstar),bracket=(1.0+Ra-v2);
    double tc[3],ta[3];
    for(int i=0;i<3;i++) tc[i]=bracket*Rv[i]-bracket*R*vel[i]/c;
    for(int i=0;i<3;i++) ta[i]=acc[i]*Rstar*R/(c*c);   /* радиационный член ВСЕГДА */
    for(int i=0;i<3;i++) out.central[i]=common*(tc[i]-ta[i]);
    vcross(n,out.central,out.vortex);
    return out;
}

/* Сила между двумя диспирациями по каналу ch. Фазовые подписи возникают из
   твиста b(θ)=cos(θ/2), а не вводятся вручную. */
void disclination_pair_force(const Disclination *F1,const Disclination *F2,
                             double t,double c,Channel ch,double out[3]){
    out[0]=out[1]=out[2]=0.0;
    for(int i=0;i<F1->n_elem;i++){
        double fpos[3],fvel[3],facc[3]; elem_kin(F1,i,t,fpos,fvel,facc);
        double b1=burgers_twist(F1,i,t);
        for(int j=0;j<F2->n_elem;j++){
            double b2=burgers_twist(F2,j,t);
            Field fld=elem_field(fpos,F2,j,t,c);
            if(ch==CH_COULOMB){
                /* трансляционный: произведение твистов × central (диполь-диполь) */
                for(int k=0;k<3;k++) out[k]+=b1*b2*fld.central[k];
            } else if(ch==CH_EXCHANGE){
                /* поворотный (спин-обмен): квадрат твистов × vortex (Франк·Франк) */
                double vxB[3]; vcross(fvel,fld.vortex,vxB);
                for(int k=0;k<3;k++) out[k]+=b1*b1*b2*b2*vxB[k]/c;
            } else { /* CH_GRAVITY */
                /* изотропный: σ~W²/c², всегда притяжение (знаконезависимо) */
                double W=vnorm(facc)/(c*c), mag=vnorm(fld.central);
                if(mag>1e-30) for(int k=0;k<3;k++) out[k]-=(b1*b1)*(b2*b2)*W*fld.central[k];
            }
        }
    }
}

/* Сила, УСРЕДНЁННАЯ по периоду прецессии T=2π/ω. Устраняет фазовый шум
   мгновенного среза: выявляет чистые степенные законы каналов.
   nT — число шагов усреднения по периоду (рекоменд. 16-24). */
void disclination_pair_force_avg(const Disclination *F1,const Disclination *F2,
                                 double c,Channel ch,int nT,double out[3]){
    out[0]=out[1]=out[2]=0.0;
    double T=2.0*M_PI/F1->omega;
    for(int it=0;it<nT;it++){
        double t=T*it/nT, f[3];
        disclination_pair_force(F1,F2,t,c,ch,f);
        out[0]+=f[0]; out[1]+=f[1]; out[2]+=f[2];
    }
    out[0]/=nT; out[1]/=nT; out[2]/=nT;
}
