// ===================== disclination_fermion_v6.h =====================
#ifndef DISCLINATION_FERMION_V6_H
#define DISCLINATION_FERMION_V6_H
/* =====================================================================
   disclination_fermion_v6.h — фермион как ОДНА ДИСПИРАЦИОННАЯ ПЕТЛЯ.
   ---------------------------------------------------------------------
   РАЗВИТИЕ v5: вместо ДВУХ точечных монополей ±q_m — ОДИН контур с
   ТВИСТОМ вектора Бюргерса b(θ)=cos(θ/2) (период 4π). Диспирация =
   дислокация (Бюргерс b, сдвиг → кулон) + дисклинация (Франк Ω=π за
   обход, поворот → спин-½, g=2).

   ТРИ КАНАЛА взаимодействия возникают ЕСТЕСТВЕННО из твиста (не вручную):
   • CH_COULOMB  — трансляционный: дипольный момент твиста, фаза ~cos(Δφ/2);
   • CH_EXCHANGE — поворотный (Франк): спин-обмен, фаза ~cos(2·...);
   • CH_GRAVITY  — изотропный: σ~⟨W²⟩/c², всегда притяжение.

   ЧИСЛЕННО ПРОВЕРЕНО (part A/B): твист cos(θ/2) имеет нулевой монополь,
   ненулевой диполь (|p|=8/3) ≡ эффективная пара ±q_m; S/L_орб=1/2 (спин-½);
   g_S=2 (заряд проходится дважды за цикл спина 4π). ⟨b1·b2⟩ идёт как
   cos(Δφ/2) — подпись периода 4π.

   Поле — дислокационная T(3) Лазара, запаздывающее с излучением.
   ===================================================================== */

typedef enum { CH_COULOMB=0, CH_EXCHANGE=1, CH_GRAVITY=2 } Channel;

/* Диспирационный фермион: петля радиуса a, прецессия ω, твист b(θ)=cos(θ/2). */
typedef struct {
    double center[3];
    double e1[3], e2[3];     /* базис плоскости петли */
    double a, omega, phase0;
    double b0;                /* амплитуда вектора Бюргерса (|b| до твиста) */
    double mu, nu, gamma;      /* упругие постоянные T(3) */
    int    n_elem;             /* дискретизация контура */
} Disclination;

/* Локальная величина Бюргерса на контуре с твистом: b(θ)=b0·cos(θ/2). */
double burgers_twist(const Disclination *D, int k, double t);

/* Мультипольные моменты твиста (проверка твист≡диполь): monopole, |dipole|. */
void twist_multipoles(const Disclination *D, double t, double *monopole, double *dipole);

/* Спин и гиромагнитное: S/L, g-фактор из периода 4π. */
double disclination_spin_ratio(const Disclination *D);   /* S/L_орб = 1/2 */
double disclination_g_factor(const Disclination *D);     /* g_S = 2 */

/* Сила между двумя диспирациями по каналу ch, момент t. */
void disclination_pair_force(const Disclination *F1, const Disclination *F2,
                             double t, double c, Channel ch, double out[3]);

/* Сила, усреднённая по периоду прецессии (nT шагов) — чистые законы каналов. */
void disclination_pair_force_avg(const Disclination *F1, const Disclination *F2,
                                 double c, Channel ch, int nT, double out[3]);

#endif
