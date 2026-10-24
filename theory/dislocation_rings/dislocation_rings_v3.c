// ===================== dislocation_rings_v3.c =====================
// ---------------------------------------------------------------------
// v3: модель частицы как ДИСЛОКАЦИОННОГО КОЛЬЦА (или диполя из двух колец)
// в калибровочной теории дислокаций Лазара-Брагинского, изоморфной
// электродинамике. Заряд кольца здесь СКАЛЯРНЫЙ (total_charge * charge_sign);
// векторный Бюргерс появляется только в v4.
//
// ПОЛЕ раскладывается на два канала (как у Лиенара-Вихерта):
//   central -- кулоно-аналог (радиальная часть, R*-формула с фактором Ферми);
//   vortex  -- аналог магнитного поля B (= n x central), даёт силу Пича-Кёлера.
// Поле ВСЕГДА запаздывающее (досветовое). Поправка Ферми -- три режима
// (PLAIN / FILM_RETARDED / FILM_INSTANT), см. dislocation_rings.h.
//
// Единицы: c -- скорость света (обычно c=1). t -- лабораторное время.
// ---------------------------------------------------------------------
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "dislocation_rings.h"

/* ---- Элементарная векторная алгебра в R^3 (все векторы -- double[3]) ---- */
static void vec_add(const double a[3], const double b[3], double out[3]) { out[0]=a[0]+b[0]; out[1]=a[1]+b[1]; out[2]=a[2]+b[2]; }   /* out = a + b */
static void vec_scale(const double a[3], double s, double out[3]) { out[0]=a[0]*s; out[1]=a[1]*s; out[2]=a[2]*s; }   /* out = s*a (масштабирование) */
static double vec_dot(const double a[3], const double b[3]) { return a[0]*b[0]+a[1]*b[1]+a[2]*b[2]; }   /* скалярное произведение a.b */
static void vec_sub(const double a[3], const double b[3], double out[3]) { out[0]=a[0]-b[0]; out[1]=a[1]-b[1]; out[2]=a[2]-b[2]; }   /* out = a - b */
static double vec_norm(const double a[3]) { return sqrt(vec_dot(a,a)); }   /* длина |a| */
static void vec_cross(const double a[3], const double b[3], double out[3]) {
    out[0]=a[1]*b[2]-a[2]*b[1]; out[1]=a[2]*b[0]-a[0]*b[2]; out[2]=a[0]*b[1]-a[1]*b[0];
}

/* Положение, скорость и ускорение k-того элемента кольца R в момент t.
   Элементы дискретизируют континуум мировых трубок кольца. */
static void ring_element_kinematics(const RingConfig *R, int k, double t, double pos[3], double vel[3], double acc[3]) {
    double theta = R->phase0 + 2.0*M_PI*k/R->n_elem + R->omega*t;
    double c1 = cos(theta), s1 = sin(theta);
    double radial[3], tang[3];
    vec_scale(R->e1, c1, radial);
    { double tmp[3]; vec_scale(R->e2, s1, tmp); vec_add(radial, tmp, radial); }
    vec_scale(radial, R->radius, radial);
    vec_add(R->center, radial, pos);
    vec_scale(R->e1, -s1, tang);
    { double tmp[3]; vec_scale(R->e2, c1, tmp); vec_add(tang, tmp, tang); }
    vec_scale(tang, R->radius*R->omega, vel);
    vec_scale(radial, -R->omega*R->omega, acc);
}

/* Ньютоновская итерация для запаздывающего времени t' между точкой
   наблюдения field_pos (в лабораторный момент t_field) и k-тым
   элементом кольца R (источник поля). Это ВСЕГДА нужно для самого
   поля (central/vortex часть) -- поле распространяется досветово. */
static double compute_retarded_time(const double field_pos[3], double t_field, const RingConfig *R, int k, double c) {
    double pos[3], vel[3], acc[3];
    ring_element_kinematics(R, k, t_field, pos, vel, acc);
    double diff[3]; vec_sub(field_pos, pos, diff);
    double Rdist = vec_norm(diff);
    double tprime = t_field - Rdist/c;
    for (int i=0;i<20;i++) {
        ring_element_kinematics(R, k, tprime, pos, vel, acc);
        vec_sub(field_pos, pos, diff);
        double Rnew = vec_norm(diff);
        double f = tprime + Rnew/c - t_field;
        double vr = vec_dot(diff, vel)/(Rnew+1e-30);
        double dRdt = -vr;
        double df = 1.0 + dRdt/c;
        if (fabs(df) < 1e-12) break;
        double delta = -f/df;
        tprime += delta;
        if (fabs(delta) < 1e-14) break;
    }
    return tprime;
}

typedef struct { double central[3]; double vortex[3]; } DislField;

/* Поле элемента источника в точке наблюдения. Структура поля (central
   = кулон-аналог, vortex = аналог B для силы Лоренца/Пича-Келлера)
   ВСЕГДА вычисляется в запаздывающем состоянии источника (tprime) --
   это отдельный, всегда досветовой уровень "поля/волны". Различие
   между режимами FILM_RETARDED/FILM_INSTANT входит ТОЛЬКО через
   a_eff_extra -- дополнительный наклон плёнки от других трубок,
   который передаётся снаружи (см. compute_ring_dipole_force). */
static DislField compute_field_from_element(const double field_pos[3], const RingConfig *R, int k,
    double t_field, double c, double a_eff_extra[3], FermiMode mode) {
    DislField out; memset(&out,0,sizeof(out));
    double tprime = compute_retarded_time(field_pos, t_field, R, k, c);
    double pos[3], vel[3], acc[3];
    ring_element_kinematics(R, k, tprime, pos, vel, acc);
    double Rvec[3]; vec_sub(field_pos, pos, Rvec);
    double Rd = vec_norm(Rvec);
    if (Rd < 1e-9) return out;
    double n[3]; vec_scale(Rvec, 1.0/Rd, n);
    double vdotn = vec_dot(vel, n);
    /* R* = R (1 - v.n/c) -- "запаздывающий" радиус Лиенара-Вихерта.
       Именно 1/R*^3 (а не 1/R^3) обеспечивает правильную релятивистскую
       фокусировку поля движущегося источника. Защита от R*->0. */
    double Rstar = Rd*(1.0 - vdotn/c);
    if (Rstar < 1e-9) Rstar = 1e-9;

    double a_eff[3] = {acc[0],acc[1],acc[2]};
    if (mode!=FERMI_NONE && a_eff_extra!=NULL) {
        a_eff[0]+=a_eff_extra[0]; a_eff[1]+=a_eff_extra[1]; a_eff[2]+=a_eff_extra[2];
    }
    /* Фактор Ферми (1 + a.R/c^2): поправка к полю от искривления мировой
       трубки ускоряющегося источника (вариация действия "типа B").
       a_eff = собственное ускорение источника + (опц.) наклон плёнки от
       других трубок (a_eff_extra) в режимах FILM_*. При FERMI_NONE = 1. */
    double fermi_factor = 1.0;
    if (mode!=FERMI_NONE) fermi_factor = 1.0 + vec_dot(a_eff, Rvec)/(c*c);

    double v2c2 = vec_dot(vel,vel)/(c*c);   /* (v/c)^2 -- релятивистское сжатие */
    double Rn_a = vec_dot(Rvec, acc)/(c*c);  /* (R.a)/c^2 -- ускорительный вклад */
    double common = 1.0/(Rstar*Rstar*Rstar); /* общий множитель 1/R*^3 */
    /* bracket = фактор Ферми * (1 + R.a/c^2 - v^2/c^2):
       скобка "скоростной" (велосити) части поля Лиенара-Вихерта. */
    double bracket = fermi_factor*(1.0 + Rn_a - v2c2);
    /* term_central: "скоростная" часть ~ (R - (R/c) v) * bracket */
    double term_central[3];
    for (int i=0;i<3;i++) term_central[i] = bracket*Rvec[i] - bracket*Rd*vel[i]/c;
    /* term_acc: "ускорительная" (радиационная) часть ~ a * R* R / c^2 */
    double term_acc[3];
    for (int i=0;i<3;i++) term_acc[i] = fermi_factor*acc[i]*Rstar*Rd/(c*c);
    /* central (кулоно-аналог) = (1/R*^3)(скоростная - ускорительная) */
    for (int i=0;i<3;i++) out.central[i] = common*(term_central[i] - term_acc[i]);
    /* vortex (аналог B) = n x central -- даёт силу Пича-Кёлера (q v x B)/c */
    double B_like[3]; vec_cross(n, out.central, B_like);
    vec_scale(B_like, 1.0, out.vortex);
    return out;
}

/* Взвешенная по обратному расстоянию интерполяция ускорений ВСЕХ
   остальных активных элементов, оцененная в момент t_eval.

   t_eval = tprime (запаздывающее время главного источника)  -->  режим
            FERMI_FILM_RETARDED: "мейнстримная" причинная версия --
            плёнка знает о других трубках только через сигнал t' <= t.

   t_eval = t_field (лабораторное время точки наблюдения)     -->  режим
            FERMI_FILM_INSTANT: ваша мгновенная нелокальная версия --
            плёнка знает "текущую" (одновременную в лаб. системе)
            конфигурацию всех трубок, независимо от того, успел ли
            сигнал дойти со скоростью c. */
static void film_effective_accel(const double field_pos[3], double t_eval, const RingConfig *rings, int n_rings,
    int exclude_ring, int exclude_k, double c, double out_accel[3]) {
    double wsum=0.0, asum[3]={0,0,0};
    for (int r=0;r<n_rings;r++) {
        const RingConfig *R = &rings[r];
        for (int k=0;k<R->n_elem;k++) {
            if (r==exclude_ring && k==exclude_k) continue;
            double pos[3], vel[3], acc[3];
            ring_element_kinematics(R, k, t_eval, pos, vel, acc);
            double diff[3]; vec_sub(field_pos, pos, diff);
            double d = vec_norm(diff); if (d < 1e-9) d = 1e-9;
            double w = 1.0/d;
            wsum += w;
            asum[0]+=w*acc[0]; asum[1]+=w*acc[1]; asum[2]+=w*acc[2];
        }
    }
    if (wsum>1e-30) { out_accel[0]=asum[0]/wsum; out_accel[1]=asum[1]/wsum; out_accel[2]=asum[2]/wsum; }
    else { out_accel[0]=out_accel[1]=out_accel[2]=0.0; }
}

void compute_ring_dipole_force(const RingConfig *field_rings, int n_field_rings,
    const RingConfig *source_rings, int n_source_rings, double t, double c, int mode_int, double out_force[3]) {
    FermiMode mode = (FermiMode)mode_int;
    out_force[0]=out_force[1]=out_force[2]=0.0;
    for (int fr=0; fr<n_field_rings; fr++) {
        const RingConfig *FR = &field_rings[fr];
        /* заряд одного элемента поля = знак * (полный заряд / число элементов) */
        double q_field = FR->charge_sign * FR->total_charge / FR->n_elem;
        for (int fk=0; fk<FR->n_elem; fk++) {
            double fpos[3], fvel[3], facc[3];
            ring_element_kinematics(FR, fk, t, fpos, fvel, facc);
            for (int sr=0; sr<n_source_rings; sr++) {
                const RingConfig *SR = &source_rings[sr];
                double q_src = SR->charge_sign * SR->total_charge / SR->n_elem;
                for (int sk=0; sk<SR->n_elem; sk++) {
                    double a_extra[3] = {0,0,0};
                    if (mode==FERMI_FILM_RETARDED) {
                        /* плёнка настраивается по ЗАПАЗДЫВАЮЩЕМУ состоянию
                           остальных трубок -- причинная, мейнстримная версия */
                        double tprime = compute_retarded_time(fpos, t, SR, sk, c);
                        film_effective_accel(fpos, tprime, source_rings, n_source_rings, sr, sk, c, a_extra);
                    } else if (mode==FERMI_FILM_INSTANT) {
                        /* плёнка настраивается по ОДНОВРЕМЕННОМУ (в лаб.
                           системе) состоянию -- ваша мгновенная нелокальная
                           версия, независимая от конуса прошлого */
                        film_effective_accel(fpos, t, source_rings, n_source_rings, sr, sk, c, a_extra);
                    }
                    /* поле k-го элемента источника в точке fk-го элемента поля */
                    DislField fld = compute_field_from_element(fpos, SR, sk, t, c, a_extra, mode);
                    /* сила = кулоновская (q_f q_s * central) + магнитная (q_f * v x B) */
                    double vxB[3]; vec_cross(fvel, fld.vortex, vxB);
                    for (int i=0;i<3;i++) out_force[i] += q_field*q_src*fld.central[i] + q_field*vxB[i];
                }
            }
        }
    }
}