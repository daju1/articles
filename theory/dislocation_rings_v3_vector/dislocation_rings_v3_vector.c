// ============== dislocation_rings_v3_vector.c ==============
// Пара скользящих (glide) колец, векторный Бюргерс, q_m=Δb в особой точке.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dislocation_rings_v3_vector.h"

static void vsub(const double a[3],const double b[3],double o[3]){o[0]=a[0]-b[0];o[1]=a[1]-b[1];o[2]=a[2]-b[2];}
static void vadd(const double a[3],const double b[3],double o[3]){o[0]=a[0]+b[0];o[1]=a[1]+b[1];o[2]=a[2]+b[2];}
static void vscale(const double a[3],double s,double o[3]){o[0]=a[0]*s;o[1]=a[1]*s;o[2]=a[2]*s;}
static double vdot(const double a[3],const double b[3]){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
static void vcross(const double a[3],const double b[3],double o[3]){
    o[0]=a[1]*b[2]-a[2]*b[1];o[1]=a[2]*b[0]-a[0]*b[2];o[2]=a[0]*b[1]-a[1]*b[0];}

/* Особая точка: бежит по кольцу (плоскость неподвижна!), угол phase0+ωt. */
void special_point(const RingV *R,double t,double pos[3],double vel[3]){
    double ph=R->phase0+R->omega*t, c=cos(ph), s=sin(ph);
    pos[0]=R->center[0]+R->a*c;
    pos[1]=R->center[1]+R->a*s;
    pos[2]=R->center[2];
    vel[0]=-R->a*R->omega*s;
    vel[1]= R->a*R->omega*c;
    vel[2]=0.0;
}

/* Дипольный момент пары: p = q_top·r_top + q_bot·r_bot, q=±delta_b. */
void pair_dipole(const RingV *top,const RingV *bot,double t,double p_out[3]){
    double pt[3],vt[3],pb[3],vb[3];
    special_point(top,t,pt,vt); special_point(bot,t,pb,vb);
    double qt=top->delta_b, qb=-bot->delta_b; /* противоположные знаки */
    /* bot->delta_b хранит |Δb|; знак минус задаёт противоположный заряд */
    for(int i=0;i<3;i++) p_out[i]=qt*pt[i]+qb*pb[i];
}

/* μ и J пары, усреднённые по периоду.
   ВАЖНО (установлено 24 окт 2026): осевой магнитный момент КОЛЬЦЕВОГО ТОКА
   пары = 0! При противофазе оба кольца бегут в одну сторону (r×v одинаков),
   но заряды противоположны (±Δb) → вклады μ_top=+0.5, μ_bot=-0.5 сокращаются.
   Поэтому mu_out здесь — это ТОКОВЫЙ момент (≈0 по оси); физический магнитный
   момент электрона в модели = ДИПОЛЬНЫЙ (см. pair_magnetic_dipole ниже).
   J_out — механический момент (складывается, ≠0). */
void pair_moments(const RingV *top,const RingV *bot,double c,
                  double mu_out[3],double J_out[3]){
    (void)c; int N=720; double T=2.0*M_PI/top->omega;
    for(int k=0;k<3;k++){mu_out[k]=0;J_out[k]=0;}
    const RingV *rings[2]; rings[0]=top; rings[1]=bot;
    double qsign[2]; qsign[0]=+1.0; qsign[1]=-1.0;
    for(int r=0;r<2;r++){
        for(int i=0;i<N;i++){
            double t=T*i/N,pos[3],vel[3];
            special_point(rings[r],t,pos,vel);
            double rxv[3]; vcross(pos,vel,rxv);
            double qm=qsign[r]*rings[r]->delta_b;
            for(int k=0;k<3;k++){ mu_out[k]+=0.5*qm*rxv[k]; J_out[k]+=rxv[k]; }
        }
    }
    for(int k=0;k<3;k++){ mu_out[k]/=N; J_out[k]/=N; }
}

/* ФИЗИЧЕСКИЙ магнитный момент пары = ДИПОЛЬНЫЙ (не токовый).
   Осевая часть p_z=q_m·d постоянна; |μ_dipole| = q_m·d (магнитный диполь ±q_m). */
double pair_magnetic_dipole(const RingV *top,const RingV *bot){
    double d=top->center[2]-bot->center[2];      /* разнос по оси */
    return top->delta_b*d;                        /* q_m·d = осевой диполь */
}

/* Масса из самодействия: m = k·μ_shear·q_m²·a/c². */
double loop_mass(const RingV *R,double c,double k){
    return k*R->mu_shear*R->delta_b*R->delta_b*R->a/(c*c);
}

/* Энергия пары E(d): притяжение Пича-Кёлера -A/d + ЦЕНТРОБЕЖНЫЙ барьер
   от сохранения момента импульса J²/(2 m d²).
   ИСПРАВЛЕНО 24 окт 2026: равновесие даёт сохранение J (как боровская орбита),
   а НЕ вихревой член силы Пича-Кёлера (central и vortex оба ~1/d², не дают
   равновесия на конечном d). d_eq = J²/(m·A). */
double pair_energy(double d,double Jcons,double m_eff,double A){
    return -A/d + Jcons*Jcons/(2.0*m_eff*d*d);
}

/* Суммарный Бюргерс (=0 при аннигиляции Френкеля). */
void total_burgers(const RingV *top,const RingV *bot,double out[3]){
    /* top интерстициал (+b), bot вакансия (-b) */
    for(int i=0;i<3;i++) out[i]=top->b_vec[i]-bot->b_vec[i];
}
