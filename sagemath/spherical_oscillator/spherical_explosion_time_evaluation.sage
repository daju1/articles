import sys
reload(sys)
sys.setdefaultencoding('utf8')

import numpy as np

attach("tzap.spyx")
attach("float_formatting.sage")

c = get_light_veloncity()
q = 1
t = 5
R0 = 2
r0 = 1
v0 = 0
a0 = 0.1
r_min = 0.1

v0_p = 0
v0_n = 0

a0_p = 0.0001
a0_n = 0.01


step_R0 = 1.0
min_R0 = r0 + step_R0
max_R0 = 7.0

t1 = 0
t2 = 10
dt = 0.001

all_plot_data_E = []
all_plot_data_E_p = []
all_plot_data_E_n = []
all_plot_data_EE = []
all_plot_data_E1 = []
all_plot_data_E1_p = []
all_plot_data_E1_n = []
all_plot_data_E2 = []
all_plot_data_E2_p = []
all_plot_data_E2_n = []
for R0_i in np.arange(min_R0, max_R0, step_R0):
    plot_data_E = []
    plot_data_E_p = []
    plot_data_E_n = []
    plot_data_EE = []
    plot_data_E1 = []
    plot_data_E1_p = []
    plot_data_E1_n = []
    plot_data_E2 = []
    plot_data_E2_p = []
    plot_data_E2_n = []
    for t_i in np.arange(t1, t2, dt):
        (phi_p, E1_p, E2_p, error_p, r_p) = phi_and_E_lw(+q, t_i, R0_i, r0, v0_p, a0_p, r_min)
        (phi_n, E1_n, E2_n, error_n, r_n) = phi_and_E_lw(-q, t_i, R0_i, r0, v0_n, a0_n, r_min)
        print (phi_p, E1_p, E2_p, error_p, r_p)
        print (phi_n, E1_n, E2_n, error_n, r_n)
        if r_n < R0_i - step_R0:
            plot_data_E += [(t_i, E1_p + E2_p + E1_n + E2_n)]
            plot_data_E_p += [(t_i, E1_p + E2_p)]
            plot_data_E_n += [(t_i, E1_n + E2_n)]
            plot_data_EE += [(t_i, E1_p - E2_p + E1_n - E2_n)]
            plot_data_E1 += [(t_i, E1_p + E1_n)]
            plot_data_E1_p += [(t_i, E1_p)]
            plot_data_E1_n += [(t_i, E1_n)]
            plot_data_E2 += [(t_i, E2_p + E2_n)]
            plot_data_E2_p += [(t_i, E2_p)]
            plot_data_E2_n += [(t_i, E2_n)]
    all_plot_data_E += plot_data_E
    all_plot_data_E_p += plot_data_E_p
    all_plot_data_E_n += plot_data_E_n
    all_plot_data_EE += plot_data_EE
    all_plot_data_E1 += plot_data_E1
    all_plot_data_E1_p += plot_data_E1_p
    all_plot_data_E1_n += plot_data_E1_n
    all_plot_data_E2 += plot_data_E2
    all_plot_data_E2_p += plot_data_E2_p
    all_plot_data_E2_n += plot_data_E2_n

    if len(plot_data_E) > 0:
        p = list_plot(plot_data_E)
        pname = "results/spherical_explosion_E_t" + "_R0=" + float_formatting(R0_i) + ".png"
        print pname
        p.save(pname)

    if len(plot_data_E_p) > 0:
        p = list_plot(plot_data_E_p)
        pname = "results/spherical_explosion_E_p_t" + "_R0=" + float_formatting(R0_i) + ".png"
        print pname
        p.save(pname)

    if len(plot_data_E_n) > 0:
        p = list_plot(plot_data_E_n)
        pname = "results/spherical_explosion_E_n_t" + "_R0=" + float_formatting(R0_i) + ".png"
        print pname
        p.save(pname)

    if len(plot_data_EE) > 0:
        p = list_plot(plot_data_EE)
        pname = "results/spherical_explosion_EE_t" + "_R0=" + float_formatting(R0_i) + ".png"
        print pname
        p.save(pname)

    if len(plot_data_E1) > 0:
        p = list_plot(plot_data_E1)
        pname = "results/spherical_explosion_E1_t" + "_R0=" + float_formatting(R0_i) + ".png"
        print pname
        p.save(pname)

    if len(plot_data_E1_p) > 0:
        p = list_plot(plot_data_E1_p)
        pname = "results/spherical_explosion_E1_p_t" + "_R0=" + float_formatting(R0_i) + ".png"
        print pname
        p.save(pname)

    if len(plot_data_E1_n) > 0:
        p = list_plot(plot_data_E1_n)
        pname = "results/spherical_explosion_E1_n_t" + "_R0=" + float_formatting(R0_i) + ".png"
        print pname
        p.save(pname)

    if len(plot_data_E2) > 0:
        p = list_plot(plot_data_E2)
        pname = "results/spherical_explosion_E2_t" + "_R0=" + float_formatting(R0_i) + ".png"
        print pname
        p.save(pname)

    if len(plot_data_E2_p) > 0:
        p = list_plot(plot_data_E2_p)
        pname = "results/spherical_explosion_E2_p_t" + "_R0=" + float_formatting(R0_i) + ".png"
        print pname
        p.save(pname)

    if len(plot_data_E2_n) > 0:
        p = list_plot(plot_data_E2_n)
        pname = "results/spherical_explosion_E2_n_t" + "_R0=" + float_formatting(R0_i) + ".png"
        print pname
        p.save(pname)

p = list_plot(all_plot_data_E)
pname = "results/spherical_explosion_all_E_t.png"
print pname
p.save(pname)

p = list_plot(all_plot_data_E_p)
pname = "results/spherical_explosion_all_E_p_t.png"
print pname
p.save(pname)

p = list_plot(all_plot_data_E_n)
pname = "results/spherical_explosion_all_E_n_t.png"
print pname
p.save(pname)

p = list_plot(all_plot_data_EE)
pname = "results/spherical_explosion_all_EE_t.png"
print pname
p.save(pname)

p = list_plot(all_plot_data_E1)
pname = "results/spherical_explosion_all_E1_t.png"
print pname
p.save(pname)

p = list_plot(all_plot_data_E1_p)
pname = "results/spherical_explosion_all_E1_p_t.png"
print pname
p.save(pname)

p = list_plot(all_plot_data_E1_n)
pname = "results/spherical_explosion_all_E1_n_t.png"
print pname
p.save(pname)

p = list_plot(all_plot_data_E2)
pname = "results/spherical_explosion_all_E2_t.png"
print pname
p.save(pname)

p = list_plot(all_plot_data_E2_p)
pname = "results/spherical_explosion_all_E2_p_t.png"
print pname
p.save(pname)

p = list_plot(all_plot_data_E2_n)
pname = "results/spherical_explosion_all_E2_n_t.png"
print pname
p.save(pname)
