// ===================== dislocation_rings_v4.h =====================
#ifndef DISLOCATION_RINGS_V4_H
#define DISLOCATION_RINGS_V4_H

/* =====================================================================
   dislocation_rings_v4.h  (v4 -- ВЕКТОРНЫЙ заряд Бюргерса + мёбиусово
   расслоение знака)
   ---------------------------------------------------------------------
   ОТЛИЧИЕ ОТ v3:
   В v1--v3 заряд кольца был СКАЛЯРНЫМ: total_charge * charge_sign.
   В v4 заряд кольца -- это ВЕКТОР БЮРГЕРСА b^a = \int_S T^a
   (Lazar 2010): топологический аналог квантованного магнитного заряда.
   Он задаётся:
       burgers_magnitude  |b|   -- модуль вектора Бюргерса,
       burgers_angle      theta_b -- ориентация b в плоскости кольца
                                      (угол относительно e1).
   Скалярный total_charge из v3 восстанавливается как частный случай
   |b| с фиксированной ориентацией; см. совместимость ниже.

   МЁБИУСОВО РАССЛОЕНИЕ (v3.tex, раздел mobius-compatibility):
   Петля остаётся обычной ориентируемой окружностью, но ВЕЛИЧИНЕ заряда
   приписывается мёбиусово 1-расслоение над контуром: при ОДНОКРАТНОМ
   обходе петли знак заряда меняется на противоположный и возвращается
   лишь при ПОВТОРНОМ обходе (антипериодические граничные условия,
   как у фермионных полей на окружности -> спинорная двузначность 4pi).
   Реализуется множителем mobius_sign(theta) = cos(theta/2):
   период 4pi, знак меняется за один обход 2pi. Включается флагом
   mobius_enable (0/1).

   КАНАЛЫ ПОЛЯ (как в v3):
     central -- кулоно-аналог (R*-формула ЛВ x множитель Ферми),
     vortex  -- аналог B (сила Пича-Кёлера), vortex = n x central.
   Поле ВСЕГДА запаздывающее (досветовое). Поправка Ферми -- три
   режима, включая вашу мгновенную нелокальную плёнку.
   ===================================================================== */

typedef enum {
    FERMI_PLAIN         = 0,
    FERMI_FILM_RETARDED = 1,
    FERMI_FILM_INSTANT  = 2,   /* ваша мгновенная нелокальная плёнка */
    FERMI_NONE          = 3
} FermiMode;

typedef struct {
    double center[3];
    double e1[3];
    double e2[3];
    double radius;
    double omega;
    double phase0;
    double burgers_angle;      /* theta_b: ориентация вектора Бюргерса в (e1,e2) */
    double burgers_magnitude;  /* |b|: модуль вектора Бюргерса (бывш. total_charge) */
    int    n_elem;
    int    charge_sign;        /* +1/-1: общий знак петли (N/S диполя) */
    int    mobius_enable;      /* 1: включить мёбиусову двузначность знака */
} RingConfig;

/* Полный вектор Бюргерса k-го элемента кольца R в лаборатории в момент t
   (ориентация поворачивается вместе с кольцом; модуль постоянен). */
void ring_element_burgers(const RingConfig *R, int k, double t, double b_out[3]);

/* Сила между набором field_rings и source_rings в момент t.
   mode_int: см. FermiMode. c: скорость света. */
void compute_ring_dipole_force(
    const RingConfig *field_rings, int n_field_rings,
    const RingConfig *source_rings, int n_source_rings,
    double t, double c, int mode_int,
    double out_force[3]);


/* Устойчивость: полная энергия кольца (кулон + натяжение) и равновесный радиус. */
double ring_total_energy(double a, double q, double r0, double k, double T,
                         double *E_coul_out, double *E_tens_out);
double ring_equilibrium_radius(double q, double r0, double k, double T,
                               double amin, double amax, int n, int *stable);

#endif /* DISLOCATION_RINGS_V4_H */
