// ===================== disclination_fermion_v7.h =====================
#ifndef DISCLINATION_FERMION_V7_H
#define DISCLINATION_FERMION_V7_H
/* =====================================================================
   disclination_fermion_v7.h — фермион как диспирационная петля с
   ВЕКТОРНЫМ твистом вектора Бюргерса b(θ) в 3D.
   ---------------------------------------------------------------------
   РАЗВИТИЕ v6: b(θ) теперь ВЕКТОР (не скаляр-амплитуда). Задаётся:
   • plane  — плоскость твиста (T_R / R_Z / T_Z / косая по normal_beta);
   • winding n — число полуоборотов за обход (спин = n/2);
   • chirality ±1 — знак закрутки (частица/античастица, слабый заряд).

   ЧИСЛЕННО УСТАНОВЛЕНО (v7 физика, 24 окт 2026):
   • Магнитный момент μ≠0 ТОЛЬКО для плоскости T_R (b в плоскости кольца,
     ось Франка Ω||ẑ). Для R_Z, T_Z — μ=0. Косой наклон: μ_z=0.5·cosβ.
     => электрон = плоскость T_R (β=0). Выбор ОПРЕДЕЛЁН физикой, не постулат.
   • μ_z ∝ спин: 0.5 при n=1, 1.5 при n=3 (μ пропорционален S=n/2).
   • μ от ТРАНСЛЯЦИОННОГО тока b = 0 (осевая симметрия) — μ живёт в
     ДИСКЛИНАЦИОННОЙ части (Франк Ω, момент-напряжение H_a).

   ШЕСТЬ каналов взаимодействия (Пич-Кёлер, Lazar 2010 ур.31):
   CH_COULOMB (1/R²), CH_GRAVITY (1/R², слабее), CH_EXCHANGE (короткий),
   CH_MAGNETIC (μ∝спин, НОВОЕ), CH_WEAK (по знаку киральности, P-нар.),
   CH_STRONG (частичный b + лента SF, линейный — ЗАГОТОВКА).
   ===================================================================== */

typedef enum { PLANE_T_R=0, PLANE_R_Z=1, PLANE_T_Z=2, PLANE_TILTED=3 } TwistPlane;
typedef enum { CH_COULOMB=0, CH_GRAVITY=1, CH_EXCHANGE=2,
               CH_MAGNETIC=3, CH_WEAK=4, CH_STRONG=5 } Channel;

typedef struct {
    double center[3];
    double a, omega, phase0;
    double b0;                 /* амплитуда |b| (дробная для кварков: b0/3) */
    int    winding;            /* n: спин = n/2 */
    int    chirality;          /* +1/-1: частица/античастица, слабый заряд */
    TwistPlane plane;          /* плоскость твиста */
    double normal_beta;        /* наклон нормали от ẑ (только для PLANE_TILTED) */
    double mu, nu, gamma;      /* упругие постоянные T(3) */
    int    n_elem;
} Disclination;

/* Вектор Бюргерса b(θ) в 3D: вращается в выбранной плоскости, фаза n·θ/2. */
void burgers_vector(const Disclination *D, double theta, double b_out[3]);
/* Вектор Франка Ω(θ) = ось поворота b = нормаль плоскости твиста × (n/2). */
void frank_vector(const Disclination *D, double theta, double O_out[3]);

/* Магнитный момент μ (из дисклинац. части, усреднён по 4π). */
void magnetic_moment(const Disclination *D, double mu_out[3]);
/* Момент количества движения J (орбитальный + спиновый). */
void angular_momentum(const Disclination *D, double c, double J_out[3]);
/* Гиромагнитное g = (μ/S)·(нормировка). Возвращает μ_z и S_z отдельно. */
double gyromagnetic_g(const Disclination *D, double c, double *mu_z, double *S_z);

/* Спин и период (из winding). */
double spin_value(const Disclination *D);        /* n/2 */
int    twist_period_2pi_factor(const Disclination *D); /* (-1)^n: -1 фермион */

/* Сила между двумя диспирациями по каналу ch, усреднённая по периоду. */
void pair_force_avg(const Disclination *F1, const Disclination *F2,
                    double c, Channel ch, int nT, double out[3]);

#endif
