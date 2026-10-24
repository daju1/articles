// ================= unified_channels.h =================
#ifndef UNIFIED_CHANNELS_H
#define UNIFIED_CHANNELS_H
/* =====================================================================
   unified_channels.h — расчёт ВСЕХ каналов взаимодействия для ВСЕХ типов
   частиц единой геометродинамической модели (дислокационные кольца).
   Строго по UNIFIED_THEORY_consistent_2026-10-24.md.

   УПРУГИЙ БАЗИС (L,T,M), без СИ-ампер:
     вектор Бюргерса b = смещение [L]; магнитный заряд q_m = разрыв b [L];
     модуль сдвига mu [Па]; c=sqrt(mu/rho); масса m=mu*b^2*a/c^2.

   ТИПЫ ЧАСТИЦ (§2):
     PARTICLE_GLIDE_DIPOLE  — пара обычных glide-колец, противофаза, период 2π (бозон);
     PARTICLE_MOBIUS_DIPOLE — пара мёбиус (кольцо+анти), период 4π (фермион=электрон);
     PARTICLE_SINGLE_RING   — одно кольцо (ссылочное, монополь-подобное);
     PARTICLE_QUARK_PARTIAL — частичный Бюргерс b/3 (сильный сектор, конфайнмент).

   ШЕСТЬ КАНАЛОВ (§5):
     CH_COULOMB   1/R²  (произведение зарядов-Бюргерсов);
     CH_MAGNETIC  1/R⁴  (диполь-диполь q_m*d, статический);
     CH_GRAVITY   1/R²  (плёнка Ферми, σ~⟨W²⟩>0, всегда притяжение);
     CH_EXCHANGE  1/R³  (фаза прецессии cos2Δφ, короткодействие);
     CH_WEAK      e^-R/R_W (киральность твиста, нарушение чётности);
     CH_STRONG    линейный V=σR (частичный b/3 + лента, конфайнмент).
   ===================================================================== */

typedef enum {
    PARTICLE_GLIDE_DIPOLE  = 0,
    PARTICLE_MOBIUS_DIPOLE = 1,
    PARTICLE_SINGLE_RING   = 2,
    PARTICLE_QUARK_PARTIAL = 3
} ParticleType;

typedef enum {
    CH_COULOMB  = 0,
    CH_MAGNETIC = 1,
    CH_GRAVITY  = 2,
    CH_EXCHANGE = 3,
    CH_WEAK     = 4,
    CH_STRONG   = 5
} Channel;

typedef struct {
    double center[3];      /* положение частицы */
    double a;              /* радиус кольца */
    double d;              /* разнос полюсов (плечо диполя) */
    double omega;          /* частота прецессии (бег особой точки) */
    double phase0;         /* нач. фаза */
    double b0;             /* амплитуда вектора Бюргерса (b0/3 для кварка) */
    int    chirality;      /* +1/-1 киральность твиста (частица/античастица) */
    double mu_shear;       /* модуль сдвига эфира */
    double R_weak;         /* радиус слабого (ħ/M_W c) в модельных ед. */
    double sigma_string;   /* натяжение струны (сильное) */
    ParticleType type;
    int    n_elem;         /* дискретизация контура */
} Particle;

/* Инициализация частицы данного типа с разумными умолчаниями. */
Particle make_particle(ParticleType type, const double center[3],
                       double a, double omega, int chirality);

/* Магнитный заряд частицы (разрыв Бюргерса, упругие ед. [L]). */
double particle_charge(const Particle *p);
/* Магнитный дипольный момент (статический, q_m*d). */
double particle_magnetic_dipole(const Particle *p);
/* Масса из самодействия: m = mu*b^2*a/c^2. */
double particle_mass(const Particle *p, double c);
/* Период конфигурации пары (2π glide / 4π mobius) → статистика. */
double particle_period(const Particle *p);
/* Спин (0 для бозона glide-пары, 1/2 для мёбиус-пары). */
double particle_spin(const Particle *p);

/* Сила данного канала между двумя частицами на расстоянии R, усреднённая
   по периоду прецессии (nT шагов). Результат — модуль силы. */
double channel_force(const Particle *p1, const Particle *p2,
                     double R, double c, Channel ch, int nT);

/* Безразмерная константа связи канала (норм. на ħc-аналог) на масштабе a. */
double channel_coupling(const Particle *p1, const Particle *p2,
                        double c, Channel ch);

#endif
