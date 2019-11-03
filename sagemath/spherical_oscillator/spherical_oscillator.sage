import sys
reload(sys)
sys.setdefaultencoding('utf8')

import numpy as np

attach("tzap.spyx")
attach("float_formatting.sage")
c = get_light_veloncity()


# Data for plotting of t_zap
def plot_spherical_capascitor_tzap_theta(q, t, min_R0, max_R0, step_R0, r0, v0, a0):
    npoints = 180
    thetas = np.arange(0*pi/npoints, npoints*pi/npoints + 1*pi/npoints, 1*pi/npoints)
    t_zap_data = [ (theta_i, tzap(q, t, R0_i, r0, v0, a0, theta_i, r_min)) for theta_i in thetas for R0_i in np.arange(min_R0, max_R0, step_R0)]

    p = list_plot (t_zap_data)
    p.save("results/spherical_oscillator_tzap_theta" + suffix(t, r0, a0) + ".png")

# Data for plotting of t_zap
def plot_spherical_capascitor_tzap_R0(q, t, min_R0, max_R0, step_R0, r0, v0, a0):
    npoints = 2
    thetas = np.arange(0*pi/npoints, npoints*pi/npoints + 1*pi/npoints, 1*pi/npoints)
    t_zap_data = [ (R0_i, tzap(q, t, R0_i, r0, v0, a0, theta_i, r_min)) for theta_i in thetas for R0_i in np.arange(min_R0, max_R0, step_R0)]

    p = list_plot (t_zap_data)
    p.save("results/spherical_oscillator_tzap_R0" + suffix(t, r0, a0) + ".png")


def plot_spherical_capascitor_phi_R0(q, t, min_R0, max_R0, step_R0):
    phi_lw_data = [ (R0_i, phi_lw(-q, t, R0_i, Rneg, 0, 0, r_min) + phi_lw(q, t, R0_i, Rpos, 0, 0, r_min)) for R0_i  in np.arange(min_R0, max_R0, step_R0) ]
    p = list_plot (phi_lw_data)
    p.save("results/spherical_capascitor_phi" + "_Rneg=" + str(Rneg) + "_Rpos= " + str(Rpos) + ".png")


def plot_spherical_oscillator_phi_R0(q, t, min_R0, max_R0, step_R0, Rneg, v0neg, a0neg, r_min_neg, Rpos, v0pos, a0pos, r_min_pos ):
    phi_lw_data = [ (R0_i, phi_lw(-q, t, R0_i, Rneg, v0neg, a0neg, r_min_neg) + phi_lw(q, t, R0_i, Rpos, v0pos, a0pos, r_min_pos)) for R0_i in np.arange(min_R0, max_R0, step_R0) ]

    p = list_plot (phi_lw_data)
    pname = "results/spherical_oscillator_phi" + \
        "_t=" + float_formatting(t) + \
        "_Rneg=" + float_formatting(Rneg) + "_Rpos= " + float_formatting(Rpos) + \
        "_v0neg=" + float_formatting(v0neg) + "_v0pos=" + float_formatting(v0pos) + \
        "_a0neg=" + float_formatting(a0neg) + "_a0pos=" + float_formatting(a0pos) + \
        ".png"
    print pname
    p.save(pname)


def plots_spherical_oscillator_phi_R0(q, t1, t2, dt, min_R0, max_R0, step_R0, Rneg, v0neg, a0neg, r_min_neg, Rpos, v0pos, a0pos, r_min_pos):
    phi_lw_data = [ (R0_i, phi_lw(-q, t_i, R0_i, Rneg, v0neg, a0neg, r_min) + phi_lw(q, t_i, R0_i, Rpos, v0pos, a0pos, r_min)) for R0_i in np.arange(min_R0, max_R0, step_R0) for t_i in np.arange(t1, t2, dt)]

    p = list_plot (phi_lw_data)
    pname = "results/spherical_oscillators_phi" + \
        "_t=" + float_formatting(t1) + ".." + float_formatting(t2) + \
        "_Rneg=" + float_formatting(Rneg) + "_Rpos= " + float_formatting(Rpos) + \
        "_v0neg=" + float_formatting(v0neg) + "_v0pos=" + float_formatting(v0pos) + \
        "_a0neg=" + float_formatting(a0neg) + "_a0pos=" + float_formatting(a0pos) + \
        ".png"
    print pname
    p.save(pname)

q = 1
r_min = 0.1
t = 0

# Data for plotting of phi_lw of unmoved spherical capacitor
min_R0 = -5.0/3.0*c
max_R0 = 5.0/3.0*c
step_R0 = 0.25/3.0*c
Rneg=2.0/3.0*c
Rpos=1.0/3.0*c

#plot_spherical_capascitor_phi_R0(q, t, min_R0, max_R0, step_R0)


# Data for plotting of p0.5hi_lw of spherical oscillator with expanding negative sphere
v0pos = 0
v0neg = 0

a0pos = 0.05*c
a0neg = -0.15*c
plot_spherical_oscillator_phi_R0(q, t, min_R0, max_R0, step_R0, Rneg, v0neg, a0neg, r_min, Rpos, v0pos, a0pos, r_min )
plot_spherical_capascitor_tzap_theta(q, t, min_R0, max_R0, step_R0, Rpos, v0pos, a0pos)
step_R0 = 0.1/3.0*c
plot_spherical_capascitor_tzap_R0(q, t, min_R0, max_R0, step_R0, Rpos, v0pos, a0pos)

t1 = 0.0
dt = 1.0
t2 = 5.0
#plots_spherical_oscillator_phi_R0(q, t1, t2, dt, min_R0, max_R0, step_R0, Rneg, v0neg, a0neg, r_min, Rpos, v0pos, a0pos, r_min)

a0pos = 0
a0neg = 0
t = 0

# Расчёт по формуле 19 скалярного потенциала Лиенара Вихерта (пока что без учёта запаздывания) для двух обкладок сферического конденсатора, в котором внешняя отрицательная обкладка разлетается наружу со скоростью 1.0/3.0, а внутренняя обкладка покоится, имеет вид. 
v0pos = 0.0
v0neg = 1.0/3.0*c
#plot_spherical_oscillator_phi_R0(q, t, min_R0, max_R0, step_R0, Rneg, v0neg, a0neg, r_min, Rpos, v0pos, a0pos, r_min )

# Скалярный потенциал Лиенара Вихерта (без учёта запаздывания) сферического конденсатора того же размера, в котором внешняя отрицательная обкладка схлопывается внутрь с той же скоростью  -1.0/3.0. Внутри внутренней положительной сферы появляется потенциальная яма для положительных зарядов.
v0pos = 0.0
v0neg = -1.0/3.0*c
#plot_spherical_oscillator_phi_R0(q, t, min_R0, max_R0, step_R0, Rneg, v0neg, a0neg, r_min, Rpos, v0pos, a0pos, r_min )

# Результат расчёта по формуле 116 скалярного потенциала Лиенара Вихерта c учётом запаздывания сферического конденсатора, в котором обе обкладки - как внешняя отрицательная  так и внутренняя положительная разлетаются наружу. Скорости обкладок  и. Начальная фаза центрально-симметричного взрыва. Для сравнения на том же графике приведён результат расчёта по формуле 19 без учёта запаздывания
v0pos = 1.0/6.0*c
v0neg = 1.0/3.0*c
#plot_spherical_oscillator_phi_R0(q, t, min_R0, max_R0, step_R0, Rneg, v0neg, a0neg, r_min, Rpos, v0pos, a0pos, r_min )

# Скалярный потенциал Лиенара Вихерта c учётом и без учёта запаздывания сферического конденсатора, в котором положительная обкладка разлетается наружу, а отрицательная обкладка схлопывается внутрь. Обе обкладки движутся навстречу друг другу. Скорости обкладок  и. Потенциальная яма внутри положительной обкладки при учёте запаздывания оказывается глубже, чем без учёта. Интересный момент, что при учёте запаздывания излом потенциальной кривой в области внутренней положительно заряженной вкладки практически полностью исчёз.
v0pos = 1.0/3.0*c
v0neg = -1.0/3.0*c
#plot_spherical_oscillator_phi_R0(q, t, min_R0, max_R0, step_R0, Rneg, v0neg, a0neg, r_min, Rpos, v0pos, a0pos, r_min )

t = 0
Rneg=1.0/3.0*c
Rpos=1.0/3.0*c
v0pos = 0
v0neg = 0
a0pos = 0.001*c
a0neg = 0.1*c
#plot_spherical_oscillator_phi_R0(q, t, min_R0, max_R0, step_R0, Rneg, v0neg, a0neg, r_min, Rpos, v0pos, a0pos, r_min )
#plots_spherical_oscillator_phi_R0(q, t1, t2, dt, min_R0, max_R0, step_R0, Rneg, v0neg, a0neg, r_min, Rpos, v0pos, a0pos, r_min)

