// ===================== positronium_v5.c =====================
// Электрон = диполь из двух прецессирующих магнитных монополей Дирака.
// Электрического заряда нет; кулон выводится из прецессии диполей.
// Поле — дислокационная T(3)-теория Лазара (central+vortex), запаздывающее.
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "positronium_v5.h"

/* ---- векторная алгебра ---- */
static void vadd(const double a[3],const double b[3],double o[3]){o[0]=a[0]+b[0];o[1]=a[1]+b[1];o[2]=a[2]+b[2];}
static void vsub(const double a[3],const double b[3],double o[3]){o[0]=a[0]-b[0];o[1]=a[1]-b[1];o[2]=a[2]-b[2];}
static void vscale(const double a[3],double s,double o[3]){o[0]=a[0]*s;o[1]=a[1]*s;o[2]=a[2]*s;}
static double vdot(const double a[3],const double b[3]){return a[0]*b[0]+a[1]*b[1]+a[2]*b[2];}
static double vnorm(const double a[3]){return sqrt(vdot(a,a));}
static void vcross(const double a[3],const double b[3],double o[3]){
    o[0]=a[1]*b[2]-a[2]*b[1];o[1]=a[2]*b[0]-a[0]*b[2];o[2]=a[0]*b[1]-a[1]*b[0];}

/* Кинематика монополя m в момент t: положение, скорость, ускорение.
   Монополь бежит по окружности прецессии радиуса a вокруг center. Знак sign
   ставит +q_m и -q_m в противоположные точки окружности (диполь). */
static void mono_kin(const Monopole *m,double t,double pos[3],double vel[3],double acc[3]){
    double th=m->phase0+m->omega*t + (m->sign>0?0.0:M_PI); /* -q_m на π впереди */
    double c=cos(th),s=sin(th),rad[3],tmp[3];
    vscale(m->e1,c,rad); vscale(m->e2,s,tmp); vadd(rad,tmp,rad);     /* орт радиуса */
    vscale(rad,m->a,tmp); vadd(m->center,tmp,pos);                   /* положение */
    vscale(m->e1,-s,vel); vscale(m->e2,c,tmp); vadd(vel,tmp,vel);
    vscale(vel,m->a*m->omega,vel);                                   /* скорость */
    vscale(rad,-m->a*m->omega*m->omega,acc);                         /* центростр. ускорение */
}

/* Запаздывающее время от наблюдателя fp (в момент tf) к монополю m. */
static double retarded_time(const double fp[3],double tf,const Monopole *m,double c){
    double pos[3],vel[3],acc[3]; mono_kin(m,tf,pos,vel,acc);
    double d[3]; vsub(fp,pos,d); double R=vnorm(d); double tp=tf-R/c;
    for(int i=0;i<30;i++){ mono_kin(m,tp,pos,vel,acc); vsub(fp,pos,d);
        double Rn=vnorm(d),f=tp+Rn/c-tf,vr=vdot(d,vel)/(Rn+1e-30),df=1.0-vr/c;
        if(fabs(df)<1e-13)break; double dt=-f/df; tp+=dt; if(fabs(dt)<1e-15)break; }
    return tp;
}

typedef struct { double central[3]; double vortex[3]; } Field;

/* Поле монополя в точке fp (дислокационная T(3), запаздывающее, с излучением).
   central — кулоно-аналог, vortex=n×central — аналог B (Пича-Кёлера).
   emode задаёт способ получения эффективного электрического вклада. */
static Field mono_field(const double fp[3],const Monopole *m,double tf,double c,EffChargeMode emode){
    Field out; memset(&out,0,sizeof(out));
    double tp=retarded_time(fp,tf,m,c);
    double pos[3],vel[3],acc[3]; mono_kin(m,tp,pos,vel,acc);
    double Rv[3]; vsub(fp,pos,Rv); double R=vnorm(Rv); if(R<1e-9)return out;
    double n[3]; vscale(Rv,1.0/R,n);
    double vn=vdot(vel,n), Rstar=R*(1.0-vn/c); if(Rstar<1e-9)Rstar=1e-9;
    double v2=vdot(vel,vel)/(c*c), Ra=vdot(Rv,acc)/(c*c);
    /* фактор Ферми — поправка границы; в режиме ECHARGE_FERMI усиливает
       эффективный заряд (вариация B), в ECHARGE_DARWIN не добавляется. */
    double fermi = (emode==ECHARGE_FERMI) ? (1.0+Ra) : 1.0;
    double common=1.0/(Rstar*Rstar*Rstar), bracket=fermi*(1.0+Ra-v2);
    double tc[3],ta[3];
    for(int i=0;i<3;i++) tc[i]=bracket*Rv[i]-bracket*R*vel[i]/c;      /* скоростная часть */
    for(int i=0;i<3;i++) ta[i]=fermi*acc[i]*Rstar*R/(c*c);           /* радиационная (ВСЕГДА) */
    for(int i=0;i<3;i++) out.central[i]=m->sign*m->qm*common*(tc[i]-ta[i]);
    vcross(n,out.central,out.vortex);
    return out;
}

/* ============ ЗАДАЧА 1: устойчивость фермиона от самоаннигиляции ============ */
double fermion_self_energy(const Fermion *F,double a,double *E_attr,double *E_prec){
    /* ИСПРАВЛЕНО (24 окт): устойчивость через СОХРАНЕНИЕ УГЛОВОГО МОМЕНТА.
       Пара +q_m,-q_m притягивается (аннигиляция): E_attr = -qm²/(2a) ~ -1/a.
       Прецессия несёт СОХРАНЯЮЩИЙСЯ угловой момент L = m_eff·a²·ω = const
       (а НЕ ω=const!). Центробежная энергия E_cf = L²/(2 m_eff a²) ~ 1/a².
       При a→0 барьер 1/a² ПОБЕЖДАЕТ притяжение 1/a → устойчивый минимум при
       a_eq = 2L²/(m_eff qm²). Это аналог того, как угловой момент удерживает
       электрон от падения на ядро. m_eff берём ПОСТОЯННОЙ (не qm²/a), иначе
       барьер тоже ~1/a и минимума нет.
       Угловой момент фиксируется начальными условиями: L = m_eff·a0²·ω,
       где a0,ω — заданные в структуре (радиус и частота прецессии). */
    double qm=F->plus.qm, w=F->plus.omega, a0=F->plus.a, mu=F->plus.mu, g=F->plus.gamma;
    double m_eff=1.0;                          /* инертная масса монополя (постоянная) */
    double L=m_eff*a0*a0*w;                     /* сохраняющийся угловой момент */
    double E1=-qm*qm/(2.0*a);                  /* притяжение пары (стягивание) ~ -1/a */
    double E2=L*L/(2.0*m_eff*a*a);             /* центробежный барьер ~ +1/a² */
    double E3=mu*g*a;                          /* натяжение дислокации (линейное) */
    if(E_attr)*E_attr=E1; if(E_prec)*E_prec=E2;
    return E1+E2+E3;
}
double fermion_equilibrium_radius(const Fermion *F,double amin,double amax,int n,int *stable){
    double best=amin,bestE=1e300; int bi=0; double *E=(double*)malloc(sizeof(double)*n);
    for(int i=0;i<n;i++){ double a=amin+(amax-amin)*i/(n-1);
        E[i]=fermion_self_energy(F,a,NULL,NULL); if(E[i]<bestE){bestE=E[i];best=a;bi=i;} }
    if(stable)*stable=(bi>0&&bi<n-1&&(E[bi+1]-2*E[bi]+E[bi-1])>0.0)?1:0;
    free(E); return best;
}

/* ====== ЗАДАЧИ 2-4: взаимодействие двух фермионов (кулон/гравитация/обмен) ====== */
/* Суммарная сила на фермион F1 со стороны F2 в момент t. Канал — по task. */
void fermion_pair_force(const Fermion *F1,const Fermion *F2,double t,double c,
                        Task task,EffChargeMode emode,double out_force[3]){
    out_force[0]=out_force[1]=out_force[2]=0.0;
    const Monopole *src[2]={&F2->plus,&F2->minus};
    const Monopole *fld[2]={&F1->plus,&F1->minus};
    for(int fi=0;fi<2;fi++){
        double fpos[3],fvel[3],facc[3]; mono_kin(fld[fi],t,fpos,fvel,facc);
        for(int si=0;si<2;si++){
            Field F=mono_field(fpos,src[si],t,c,emode);
            double q1=fld[fi]->sign*fld[fi]->qm, q2=src[si]->sign*src[si]->qm;
            if(task==TASK_COULOMB_PAIR){
                /* кулон: знаковое произведение монополей × central (даёт cosΔφ/R) */
                for(int i=0;i<3;i++) out_force[i]+=q1*q2*F.central[i];
                double vxB[3]; vcross(fvel,F.vortex,vxB);
                for(int i=0;i<3;i++) out_force[i]+=q1/c*vxB[i];
            } else if(task==TASK_GRAVITY_PAIR){
                /* гравитация: изотропный знаконезависимый источник σ~W²/c².
                   Берём |центр.поле| × (ускорение источника)² — всегда притяжение. */
                double Wsrc=vnorm(facc)/(c*c);   /* собственное ускорение/c² */
                double mag=vnorm(F.central);
                for(int i=0;i<3;i++) out_force[i]-= (q1*q1)*(q2*q2)*Wsrc*mag*F.central[i]/(mag+1e-30);
            } else if(task==TASK_EXCHANGE){
                /* обмен: член cos2Δφ — проекция на удвоенную разность фаз.
                   Моделируем как квадрупольную (по фазе) часть central. */
                double dphi=(fld[fi]->phase0 - src[si]->phase0);
                double w2=cos(2.0*dphi);
                for(int i=0;i<3;i++) out_force[i]+=q1*q2*w2*F.central[i];
            }
        }
    }
}

/* ============ ЗАДАЧА 5: боровская орбита монополей (безызлучательность) ============ */
/* Условие стоячей волны на окружности прецессии: n·λ = 2π a, где λ=2πc/ω
   (длина волны поля на частоте прецессии). Резонанс (резидуал=0) =
   безызлучательная орбита: излучение гасится интерференцией по контуру. */
double bohr_orbit_residual(const Fermion *F,double c,int n){
    double a=F->plus.a, w=F->plus.omega;
    double lambda=2.0*M_PI*c/w;            /* длина волны поля */
    double circ=2.0*M_PI*a;                /* длина окружности прецессии */
    return circ - (double)n*lambda;        /* 0 => точный резонанс */
}
