// ================= unified_channels.c =================
// Расчёт всех каналов взаимодействия для всех типов частиц.
// Строго по UNIFIED_THEORY_consistent_2026-10-24.md (упругий базис).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "unified_channels.h"

/* ---------- векторная алгебра ---------- */
static double vdot(const double a[3],const double b[3]){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
static double vnorm(const double a[3]){return sqrt(vdot(a,a));}

/* ---------- инициализация частицы ---------- */
Particle make_particle(ParticleType type,const double center[3],
                       double a,double omega,int chirality){
    Particle p; memset(&p,0,sizeof(p));
    p.center[0]=center[0]; p.center[1]=center[1]; p.center[2]=center[2];
    p.a=a; p.d=a; p.omega=omega; p.phase0=0.0; p.chirality=chirality;
    p.mu_shear=1.0; p.R_weak=0.05; p.sigma_string=0.05;
    p.type=type; p.n_elem=48;
    switch(type){
        case PARTICLE_GLIDE_DIPOLE:  p.b0=1.0; break;
        case PARTICLE_MOBIUS_DIPOLE: p.b0=1.0; break;
        case PARTICLE_SINGLE_RING:   p.b0=1.0; break;
        case PARTICLE_QUARK_PARTIAL: p.b0=1.0/3.0; break; /* дробный Бюргерс */
    }
    return p;
}

double particle_charge(const Particle *p){ return p->b0; }           /* q_m=b [L] */
double particle_magnetic_dipole(const Particle *p){ return p->b0*p->d; } /* q_m*d */
double particle_mass(const Particle *p,double c){
    return p->mu_shear*p->b0*p->b0*p->a/(c*c);                        /* m=mu b² a/c² */
}
double particle_period(const Particle *p){
    /* мёбиус-пара и одиночный мёбиус → 4π; glide-пара → 2π */
    if(p->type==PARTICLE_MOBIUS_DIPOLE) return 4.0*M_PI;
    return 2.0*M_PI;
}
double particle_spin(const Particle *p){
    return (p->type==PARTICLE_MOBIUS_DIPOLE)? 0.5 : 0.0;
}

/* Твист знака Бюргерса: мёбиус cos(chir·θ/2) период 4π; glide — постоянный. */
static double burgers_twist(const Particle *p,double theta){
    if(p->type==PARTICLE_MOBIUS_DIPOLE) return cos(p->chirality*theta/2.0);
    return 1.0;
}
/* Фаза особой точки в момент t. */
static double sp_phase(const Particle *p,double t){ return p->phase0+p->omega*t; }
/* Собственное ускорение особой точки (прецессия): W=ω²a. */
static double self_accel(const Particle *p){ return p->omega*p->omega*p->a; }

/* ----- мгновенная сила канала между двумя частицами (ось x, разнос R) ----- */
static double channel_force_inst(const Particle *p1,const Particle *p2,
                                 double R,double c,Channel ch,double t){
    double q1=p1->b0, q2=p2->b0;
    double ph1=sp_phase(p1,t), ph2=sp_phase(p2,t)+M_PI; /* противофаза */
    double dphi=ph1-ph2;
    switch(ch){
        case CH_COULOMB: {
            /* произведение зарядов-Бюргерсов, 1/R², фазовый множитель прецессии */
            double tw=burgers_twist(p1,ph1)*burgers_twist(p2,ph2);
            return fabs(q1*q2*tw)/(R*R);
        }
        case CH_MAGNETIC: {
            /* диполь-диполь 1/R⁴ (статические магн. диполи q_m*d) */
            double mu1=q1*p1->d, mu2=q2*p2->d;
            return fabs(mu1*mu2)/(R*R*R*R);
        }
        case CH_GRAVITY: {
            /* плёнка Ферми: источник σ=mu*b²*⟨W²⟩/c² >0, сила σ1σ2/(4π γ R²) */
            double W1=self_accel(p1), W2=self_accel(p2);
            double s1=p1->mu_shear*q1*q1*W1*W1/(c*c);
            double s2=p2->mu_shear*q2*q2*W2*W2/(c*c);
            double gamma=p1->mu_shear; /* жёсткость плёнки (упругая) */
            return (s1*s2)/(4.0*M_PI*gamma*R*R);
        }
        case CH_EXCHANGE: {
            /* фаза прецессии cos2Δφ, короткодействие 1/R³ */
            return fabs(q1*q1*q2*q2*cos(2.0*dphi))/(R*R*R);
        }
        case CH_WEAK: {
            /* киральность: связь зависит от знака (нарушение чётности); e^-R/R_W */
            double chir=(double)(p1->chirality*p2->chirality);
            double proj=(1.0+chir)/2.0;            /* 1 для одинаковой кирал., 0 разной */
            return proj*fabs(q1*q2)*exp(-R/p1->R_weak)/R;
        }
        case CH_STRONG: {
            /* конфайнмент: активен для дробного Бюргерса (кварк), линейный V=σR */
            double frac=(p1->b0<0.9 && p2->b0<0.9)? 1.0 : 0.0;
            return frac*p1->sigma_string*R;        /* F=σ=const (линейный потенциал) */
        }
    }
    return 0.0;
}

/* Сила канала, усреднённая по периоду прецессии. */
double channel_force(const Particle *p1,const Particle *p2,
                     double R,double c,Channel ch,int nT){
    double T=2.0*M_PI/p1->omega, sum=0.0;
    for(int i=0;i<nT;i++){
        double t=T*i/nT;
        sum+=channel_force_inst(p1,p2,R,c,ch,t);
    }
    return sum/nT;
}

/* Безразмерная константа связи канала на масштабе a (норм. ħc-аналог = mu*a³/... ). */
double channel_coupling(const Particle *p1,const Particle *p2,double c,Channel ch){
    double R=p1->a; /* на характерном масштабе */
    double F=channel_force(p1,p2,R,c,ch,24);
    double hbar_analog=p1->mu_shear*p1->a*p1->a*p1->a/c; /* ħc-аналог в упругих ед. */
    return F*R*R/(hbar_analog*c+1e-300);
}
