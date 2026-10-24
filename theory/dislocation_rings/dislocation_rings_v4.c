// ===================== dislocation_rings_v4.c =====================
// v4: векторный заряд Бюргерса b^a + мёбиусово расслоение знака.
// Обобщает dislocation_rings_v3.c: скалярный total_charge -> вектор |b|,theta_b.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dislocation_rings_v4.h"

/* ---- Векторная алгебра в R^3 (компактные inline-хелперы) ---- */
static void vec_add(const double a[3], const double b[3], double o[3]){o[0]=a[0]+b[0];o[1]=a[1]+b[1];o[2]=a[2]+b[2];}  /* o=a+b */
static void vec_scale(const double a[3], double s, double o[3]){o[0]=a[0]*s;o[1]=a[1]*s;o[2]=a[2]*s;}
static double vec_dot(const double a[3], const double b[3]){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
static void vec_sub(const double a[3], const double b[3], double o[3]){o[0]=a[0]-b[0];o[1]=a[1]-b[1];o[2]=a[2]-b[2];}
static double vec_norm(const double a[3]){return sqrt(vec_dot(a,a));}
static void vec_cross(const double a[3], const double b[3], double o[3]){
    o[0]=a[1]*b[2]-a[2]*b[1]; o[1]=a[2]*b[0]-a[0]*b[2]; o[2]=a[0]*b[1]-a[1]*b[0];
}

/* Мёбиусов знак величины заряда вдоль контура: период 4pi, один обход 2pi
   меняет знак. cos(theta/2): антипериодические граничные условия. */
static double mobius_sign(const RingConfig *R, double theta){
    /* Мёбиусов знак величины заряда вдоль контура.
       cos(theta/2): период 4pi, НЕ 2pi. За один обход (theta: 0->2pi) знак
       меняется на противоположный, возврат к исходному лишь за ДВА обхода.
       Это антипериодические граничные условия -> спинорная двузначность
       (поворот на 4pi для возврата волновой функции фермиона).
       Если mobius_enable=0, обычный ориентируемый контур (знак постоянен). */
    if(!R->mobius_enable) return 1.0;
    return cos(0.5*theta);
}

/* Кинематика k-го элемента (как в v3). */
/* Кинематика k-го элемента кольца в момент t:
   theta -- угловое положение элемента (нач. фаза + k-я позиция + вращение omega*t).
   pos = центр + a*(e1 cos + e2 sin); vel -- касательная; acc -- центростремительное. */
static void ring_element_kinematics(const RingConfig *R,int k,double t,double pos[3],double vel[3],double acc[3]){
    double theta=R->phase0+2.0*M_PI*k/R->n_elem+R->omega*t;
    double c1=cos(theta),s1=sin(theta);
    double radial[3],tang[3],tmp[3];
    vec_scale(R->e1,c1,radial); vec_scale(R->e2,s1,tmp); vec_add(radial,tmp,radial);
    vec_scale(radial,R->radius,radial); vec_add(R->center,radial,pos);
    vec_scale(R->e1,-s1,tang); vec_scale(R->e2,c1,tmp); vec_add(tang,tmp,tang);
    vec_scale(tang,R->radius*R->omega,vel);
    vec_scale(radial,-R->omega*R->omega,acc);
}

/* Полный вектор Бюргерса k-го элемента в лаборатории: ориентация theta_b
   задана в базисе (e1,e2) и поворачивается вместе с кольцом (фаза theta). */
void ring_element_burgers(const RingConfig *R,int k,double t,double b_out[3]){
    double theta=R->phase0+2.0*M_PI*k/R->n_elem+R->omega*t;
    double ang=theta+R->burgers_angle;   /* поворот ориентации b вместе с петлёй */
    double mag=R->burgers_magnitude*R->charge_sign*mobius_sign(R,theta);
    double tmp[3];
    vec_scale(R->e1,mag*cos(ang),b_out); vec_scale(R->e2,mag*sin(ang),tmp); vec_add(b_out,tmp,b_out);
}

/* Скалярная "величина заряда" k-го элемента для central-канала:
   проекция Бюргерса на направление, релевантное кулону (его модуль со знаком).
   В пределе burgers_angle=0, mobius_enable=0 сводится к v3: |b|*charge_sign. */
/* Скалярная "величина заряда" k-го элемента для central-канала:
   |b| * знак * мёбиусов_множитель / число_элементов.
   В пределе mobius_enable=0 сводится к v3 (|b| играет роль total_charge). */
static double element_charge_scalar(const RingConfig *R,int k,double t){
    double theta=R->phase0+2.0*M_PI*k/R->n_elem+R->omega*t;
    return R->burgers_magnitude*R->charge_sign*mobius_sign(R,theta)/R->n_elem;
}

static double compute_retarded_time(const double fp[3],double tf,const RingConfig *R,int k,double c){
    double pos[3],vel[3],acc[3]; ring_element_kinematics(R,k,tf,pos,vel,acc);
    double diff[3]; vec_sub(fp,pos,diff); double Rd=vec_norm(diff); double tp=tf-Rd/c;
    for(int i=0;i<20;i++){
        ring_element_kinematics(R,k,tp,pos,vel,acc); vec_sub(fp,pos,diff);
        double Rn=vec_norm(diff); double f=tp+Rn/c-tf;
        double vr=vec_dot(diff,vel)/(Rn+1e-30); double df=1.0-vr/c;
        if(fabs(df)<1e-12) break; double d=-f/df; tp+=d; if(fabs(d)<1e-14) break;
    }
    return tp;
}

typedef struct { double central[3]; double vortex[3]; } DislField;

static DislField compute_field_from_element(const double fp[3],const RingConfig *R,int k,
    double tf,double c,double a_extra[3],FermiMode mode){
    DislField out; memset(&out,0,sizeof(out));
    double tprime=compute_retarded_time(fp,tf,R,k,c);
    double pos[3],vel[3],acc[3]; ring_element_kinematics(R,k,tprime,pos,vel,acc);
    double Rvec[3]; vec_sub(fp,pos,Rvec); double Rd=vec_norm(Rvec);
    if(Rd<1e-9) return out;
    double n[3]; vec_scale(Rvec,1.0/Rd,n);
    /* R* = R(1 - v.n/c): запаздывающий радиус Лиенара-Вихерта (защита R*->0) */
    double vdotn=vec_dot(vel,n); double Rstar=Rd*(1.0-vdotn/c); if(Rstar<1e-9) Rstar=1e-9;
    double a_eff[3]={acc[0],acc[1],acc[2]};
    if(mode!=FERMI_NONE && a_extra!=NULL){a_eff[0]+=a_extra[0];a_eff[1]+=a_extra[1];a_eff[2]+=a_extra[2];}
    /* фактор Ферми (1 + a_eff.R/c^2): a_eff = ускорение источника + наклон плёнки */
    double fermi=1.0; if(mode!=FERMI_NONE) fermi=1.0+vec_dot(a_eff,Rvec)/(c*c);
    double v2c2=vec_dot(vel,vel)/(c*c); double Rn_a=vec_dot(Rvec,acc)/(c*c);
    double common=1.0/(Rstar*Rstar*Rstar); double bracket=fermi*(1.0+Rn_a-v2c2);
    double term_c[3],term_a[3];
    for(int i=0;i<3;i++) term_c[i]=bracket*Rvec[i]-bracket*Rd*vel[i]/c;
    for(int i=0;i<3;i++) term_a[i]=fermi*acc[i]*Rstar*Rd/(c*c);
    /* central = (1/R*^3)(скоростная term_c - ускорительная/радиационная term_a).
       ВАЖНО (соответствие теории): ПОЛЕ Лиенара-Вихерта ВСЕГДА полное и
       запаздывающее -- с радиационным членом term_a -- во ВСЕХ режимах,
       включая FERMI_FILM_INSTANT. Поле есть результат минимизации САМОГО
       интеграла действия (ур. Максвелла / потенциалы ЛВ / силы Кулона и
       Пича-Кёлера, см. Braginsky_12_2018), оно запаздывает и излучает.
       Мгновенность и нелокальность относятся ТОЛЬКО к поправке Ферми --
       наклону мыльной плёнки границы действия, который собирается в
       a_extra (в compute_ring_dipole_force): для FILM_INSTANT ускорения
       чужих трубок берутся в ОДНОВРЕМЕННЫЙ лабораторный момент t, для
       FILM_RETARDED -- в запаздывающий t'. Поправка Ферми НЕ запрещает
       излучение поля; она лишь мгновенно подстраивает границу интегрирования. */
    for(int i=0;i<3;i++) out.central[i]=common*(term_c[i]-term_a[i]);
    /* vortex = n x central -- аналог B (сила Пича-Кёлера) */
    vec_cross(n,out.central,out.vortex);
    return out;
}

static void film_effective_accel(const double fp[3],double teval,const RingConfig *rings,int n_rings,
    int excl_r,int excl_k,double c,double out_accel[3]){
    (void)c; double wsum=0.0,asum[3]={0,0,0};
    for(int r=0;r<n_rings;r++){ const RingConfig *R=&rings[r];
        for(int k=0;k<R->n_elem;k++){ if(r==excl_r&&k==excl_k) continue;
            double pos[3],vel[3],acc[3]; ring_element_kinematics(R,k,teval,pos,vel,acc);
            double diff[3]; vec_sub(fp,pos,diff); double d=vec_norm(diff); if(d<1e-9) d=1e-9;
            double w=1.0/d; wsum+=w; asum[0]+=w*acc[0]; asum[1]+=w*acc[1]; asum[2]+=w*acc[2];
        } }
    if(wsum>1e-30){out_accel[0]=asum[0]/wsum;out_accel[1]=asum[1]/wsum;out_accel[2]=asum[2]/wsum;}
    else{out_accel[0]=out_accel[1]=out_accel[2]=0.0;}
}

void compute_ring_dipole_force(const RingConfig *field_rings,int n_field_rings,
    const RingConfig *source_rings,int n_source_rings,double t,double c,int mode_int,double out_force[3]){
    FermiMode mode=(FermiMode)mode_int;
    out_force[0]=out_force[1]=out_force[2]=0.0;
    for(int fr=0;fr<n_field_rings;fr++){ const RingConfig *FR=&field_rings[fr];
        for(int fk=0;fk<FR->n_elem;fk++){
            double q_field=element_charge_scalar(FR,fk,t);
            double fpos[3],fvel[3],facc[3]; ring_element_kinematics(FR,fk,t,fpos,fvel,facc);
            for(int sr=0;sr<n_source_rings;sr++){ const RingConfig *SR=&source_rings[sr];
                for(int sk=0;sk<SR->n_elem;sk++){
                    double q_src=element_charge_scalar(SR,sk,t);
                    double a_extra[3]={0,0,0};
                    if(mode==FERMI_FILM_RETARDED){
                        double tp=compute_retarded_time(fpos,t,SR,sk,c);
                        film_effective_accel(fpos,tp,source_rings,n_source_rings,sr,sk,c,a_extra);
                    } else if(mode==FERMI_FILM_INSTANT){
                        film_effective_accel(fpos,t,source_rings,n_source_rings,sr,sk,c,a_extra);
                    }
                    DislField fld=compute_field_from_element(fpos,SR,sk,t,c,a_extra,mode);
                    double vxB[3]; vec_cross(fvel,fld.vortex,vxB); vec_scale(vxB,1.0/c,vxB);
                    for(int i=0;i<3;i++) out_force[i]+=q_field*q_src*fld.central[i]+q_field*vxB[i];
                }
            }
        }
    }
}


/* =====================================================================
   УСТОЙЧИВОСТЬ (добавлено): линейное натяжение дислокации.
   Полная энергия заряженного кольца радиуса a с core r0, зарядом q,
   кулоновской постоянной k и натяжением T=mu*b^2:
       E(a) = (k q^2 /(pi a)) ln(8a/r0)   [кулоновская самоэнергия, ~1/a]
            + T * 2*pi*a                   [натяжение дислокации, ~a]
   Устойчивый минимум: a_eq ~ sqrt(k q^2/(2 pi T)), E''(a_eq)>0.
   Это недостающий отталкивающий механизм (в прежних v1-v3 его не было).
   ===================================================================== */
double ring_total_energy(double a, double q, double r0, double k, double T,
                         double *E_coul_out, double *E_tens_out) {
    double E_coul = (k*q*q/(M_PI*a))*log(8.0*a/r0);
    double E_tens = T*2.0*M_PI*a;
    if (E_coul_out) *E_coul_out = E_coul;
    if (E_tens_out) *E_tens_out = E_tens;
    return E_coul + E_tens;
}

/* Численный поиск равновесного радиуса (минимум полной энергии) на [amin,amax].
   Возвращает a_eq; *stable=1 если это устойчивый минимум (E''>0). */
double ring_equilibrium_radius(double q, double r0, double k, double T,
                               double amin, double amax, int n, int *stable) {
    double best_a=amin, best_E=1e300; int bi=0;
    double *Es=(double*)malloc(sizeof(double)*n);
    for (int i=0;i<n;i++){
        double a=amin+(amax-amin)*i/(n-1);
        Es[i]=ring_total_energy(a,q,r0,k,T,NULL,NULL);
        if (Es[i]<best_E){best_E=Es[i];best_a=a;bi=i;}
    }
    if (stable) *stable = (bi>0 && bi<n-1 && (Es[bi+1]-2*Es[bi]+Es[bi-1])>0.0)?1:0;
    free(Es);
    return best_a;
}
