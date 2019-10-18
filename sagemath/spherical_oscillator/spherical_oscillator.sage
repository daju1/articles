import sys
reload(sys)
sys.setdefaultencoding('utf8')

import numpy as np

attach("tzap.spyx")
attach("float_formatting.sage")


def plot_spherical_capascitor_phi_R0(q, t, min_R0, max_R0, step_R0):
    phi_lw_data = [ (R0_i, phi_lw(-q, t, R0_i, Rneg, 0, 0, r_min) + phi_lw(q, t, R0_i, Rpos, 0, 0, r_min)) for R0_i  in np.arange(min_R0, max_R0, step_R0) ]
    p = list_plot (phi_lw_data)
    p.save("results/spherical_capascitor_phi" + "_Rneg=" + str(Rneg) + "_Rpos= " + str(Rpos) + ".png")


def plot_spherical_oscillator_phi_R0(q, t, min_R0, max_R0, step_R0, Rneg, v0neg, a0neg, r_min_neg, Rpos, v0pos, a0pos, r_min_pos ):
    phi_lw_data = [ (R0_i, phi_lw(-q, t, R0_i, Rneg, v0neg, a0neg, r_min_neg) + phi_lw(q, t, R0_i, Rpos, v0pos, a0pos, r_min_pos)) for R0_i in np.arange(min_R0, max_R0, step_R0) ]

    p = list_plot (phi_lw_data)
    pname = "results/spherical_oscillator_phi" + "_t=" + float_formatting(t) + "_Rneg=" + float_formatting(Rneg) + "_Rpos= " + float_formatting(Rpos) + "_a0neg=" + float_formatting(a0neg) + "_a0pos=" + float_formatting(a0pos) + ".png"
    print pname
    p.save(pname)


def plots_spherical_oscillator_phi_R0(q, t1, t2, dt, min_R0, max_R0, step_R0, Rneg, v0neg, a0neg, r_min_neg, Rpos, v0pos, a0pos, r_min_pos):
    phi_lw_data = [ (R0_i, phi_lw(-q, t_i, R0_i, Rneg, v0neg, a0neg, r_min) + phi_lw(q, t_i, R0_i, Rpos, v0pos, a0pos, r_min)) for R0_i in np.arange(min_R0, max_R0, step_R0) for t_i in np.arange(t1, t2, dt)]

    p = list_plot (phi_lw_data)
    pname = "results/spherical_oscillator_phi" + "_t=" + float_formatting(t1) + ".." + float_formatting(t2) + "_Rneg=" + float_formatting(Rneg) + "_Rpos= " + float_formatting(Rpos) + "_a0neg=" + float_formatting(a0neg) + "_a0pos=" + float_formatting(a0pos) + ".png"
    print pname
    p.save(pname)

q = 1
r_min = 0.1
t = 0.2

# Data for plotting of phi_lw of unmoved spherical capacitor
min_R0 = -20.0
max_R0 = 20.0
step_R0 = 0.1
Rneg=2
Rpos=1

plot_spherical_capascitor_phi_R0(q, t, min_R0, max_R0, step_R0)


# Data for plotting of phi_lw of spherical oscillator with expanding negative sphere
v0pos = 0
v0neg = 0

a0pos = 0.05
a0neg = -0.15
plot_spherical_oscillator_phi_R0(q, t, min_R0, max_R0, step_R0, Rneg, v0neg, a0neg, r_min, Rpos, v0pos, a0pos, r_min )

t1 = 0.0
dt = 0.5
t2 = 20.0
plots_spherical_oscillator_phi_R0(q, t1, t2, dt, min_R0, max_R0, step_R0, Rneg, v0neg, a0neg, r_min, Rpos, v0pos, a0pos, r_min)

a0pos = 0
a0neg = 0
t = 0

# Результат расчёта по формуле 116 скалярного потенциала Лиенара Вихерта c учётом запаздывания сферического конденсатора, в котором обе обкладки - как внешняя отрицательная  так и внутренняя положительная разлетаются наружу. Скорости обкладок  и. Начальная фаза центрально-симметричного взрыва. Для сравнения на том же графике приведён результат расчёта по формуле 19 без учёта запаздывания
v0pos = 1/6
v0neg = 1/3
plot_spherical_oscillator_phi_R0(q, t, min_R0, max_R0, step_R0, Rneg, v0neg, a0neg, r_min, Rpos, v0pos, a0pos, r_min )

# Скалярный потенциал Лиенара Вихерта c учётом и без учёта запаздывания сферического конденсатора, в котором положительная обкладка разлетается наружу, а отрицательная обкладка схлопывается внутрь. Обе обкладки движутся навстречу друг другу. Скорости обкладок  и. Потенциальная яма внутри положительной обкладки при учёте запаздывания оказывается глубже, чем без учёта. Интересный момент, что при учёте запаздывания излом потенциальной кривой в области внутренней положительно заряженной вкладки практически полностью исчёз.
v0pos = 1/3
v0neg = -2/3
plot_spherical_oscillator_phi_R0(q, t, min_R0, max_R0, step_R0, Rneg, v0neg, a0neg, r_min, Rpos, v0pos, a0pos, r_min )
