// ============== mobius_dipole_v3B.c ==============
// Ветка B: диполь из пары мёбиус-колец (кольцо N + анти-мёбиус S).
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "mobius_dipole_v3B.h"

void mb_special_point(const MobiusRing *R,double t,double pos[3],double vel[3]){
    double ph=R->phase0+R->omega*t, c=cos(ph), s=sin(ph);
    pos[0]=R->center[0]+R->a*c;
    pos[1]=R->center[1]+R->a*s;
    pos[2]=R->center[2];
    vel[0]=-R->a*R->omega*s;
    vel[1]= R->a*R->omega*c;
    vel[2]=0.0;
}
/* Твист знака Бюргерса: cos(chirality·(phase0+ωt)/2), период 4π. */
double mb_burgers_sign(const MobiusRing *R,double t){
    double ph=R->phase0+R->omega*t;
    return cos(R->chirality*ph/2.0);
}
/* Период пары: мёбиус-звенья возвращают знак за 4π. */
double mb_pair_period(const MobiusRing *N,const MobiusRing *S){
    (void)N;(void)S; return 4.0*M_PI;
}
/* Магнитный диполь = q_m·d (осевой, статический). */
double mb_pair_magnetic_dipole(const MobiusRing *N,const MobiusRing *S){
    double d=N->center[2]-S->center[2];
    return N->delta_b*d;
}
/* Энергия пары: притяжение -A/d + центробежный барьер (защита прецессией). */
double mb_pair_energy(double d,double Jcons,double m_eff,double A){
    return -A/d + Jcons*Jcons/(2.0*m_eff*d*d);
}
/* Барьер разворота твиста ∫(db/dθ)² dθ по 4π (двухступенчатая аннигиляция). */
double mb_twist_unwind_barrier(const MobiusRing *R){
    int N=1000; double dth=4.0*M_PI/N, sum=0.0;
    for(int i=0;i<N;i++){
        double th=dth*i;
        double dcard=-0.5*R->chirality*sin(R->chirality*th/2.0); /* d/dθ cos(chir·θ/2) */
        sum+=dcard*dcard*dth;
    }
    return sum;
}
