import sys
reload(sys)
sys.setdefaultencoding('utf8')

import numpy as np

attach("tzap.spyx")
attach("float_formatting.sage")

q = 1
t = 5
R0 = 2
r0 = 1
v0 = 0
a0 = 0.1
r_min = 0.1


(phi, E1, E2, error, r) = phi_and_E_lw(q, t, R0, r0, v0, a0, r_min)
print phi, E1, E2, error, r

step_R0 = 1.0
min_R0 = r0 + step_R0
max_R0 = 10.0

t1 = 0
t2 = 100
dt = 0.1

v0_p = 0
v0_n = 0


a0_p = 0.0001
a0_n = 0.01


all_plot_data_E = []
all_plot_data_EE = []
all_plot_data_E1 = []
all_plot_data_E2 = []
for R0_i in np.arange(min_R0, max_R0, step_R0):
    plot_data_E = []
    plot_data_EE = []
    # plot_data_rn = []
    # plot_data_rp = []
    plot_data_E1 = []
    plot_data_E2 = []
    for t_i in np.arange(t1, t2, dt):
        (phi_p, E1_p, E2_p, error_p, r_p) = phi_and_E_lw(+q, t_i, R0_i, r0, v0_p, a0_p, r_min)
        (phi_n, E1_n, E2_n, error_n, r_n) = phi_and_E_lw(-q, t_i, R0_i, r0, v0_n, a0_n, r_min)
        print (phi_p, E1_p, E2_p, error_p, r_p)
        print (phi_n, E1_n, E2_n, error_n, r_n)
        if r_n < R0_i - step_R0:
            plot_data_E += [(t_i, E1_p + E2_p + E1_n + E2_n)]
            plot_data_EE += [(t_i, E1_p - E2_p + E1_n - E2_n)]
            plot_data_E1 += [(t_i, E1_p + E1_n)]
            plot_data_E2 += [(t_i, E2_p + E2_n)]
    all_plot_data_E += plot_data_E
    all_plot_data_EE += plot_data_EE
    all_plot_data_E1 += plot_data_E1
    all_plot_data_E2 += plot_data_E2

    p = list_plot(plot_data_E)
    pname = "results/spherical_explosion_E" + "_R0=" + float_formatting(R0_i) + ".png"
    print pname
    p.save(pname)

    p = list_plot(plot_data_EE)
    pname = "results/spherical_explosion_EE" + "_R0=" + float_formatting(R0_i) + ".png"
    print pname
    p.save(pname)

    p = list_plot(plot_data_E1)
    pname = "results/spherical_explosion_E1" + "_R0=" + float_formatting(R0_i) + ".png"
    print pname
    p.save(pname)

    p = list_plot(plot_data_E2)
    pname = "results/spherical_explosion_E2" + "_R0=" + float_formatting(R0_i) + ".png"
    print pname
    p.save(pname)

p = list_plot(all_plot_data_E)
pname = "results/spherical_explosion_all_E_t.png"
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

p = list_plot(all_plot_data_E2)
pname = "results/spherical_explosion_all_E2_t.png"
print pname
p.save(pname)



all_plot_data_E = []
all_plot_data_EE = []
all_plot_data_E1 = []
all_plot_data_E2 = []
for t_i in np.arange(t1, t2, dt):
    plot_data_E = []
    plot_data_EE = []
    plot_data_E1 = []
    plot_data_E2 = []
    for R0_i in np.arange(min_R0, max_R0, step_R0):
        (phi_p, E1_p, E2_p, error_p, r_p) = phi_and_E_lw(+q, t_i, R0_i, r0, v0_p, a0_p, r_min)
        (phi_n, E1_n, E2_n, error_n, r_n) = phi_and_E_lw(-q, t_i, R0_i, r0, v0_n, a0_n, r_min)
        print (phi_p, E1_p, E2_p, error_p, r_p)
        print (phi_n, E1_n, E2_n, error_n, r_n)
        if r_n < R0_i - step_R0:
            plot_data_E += [(R0_i, E1_p + E2_p + E1_n + E2_n)]
            plot_data_EE += [(R0_i, E1_p - E2_p + E1_n - E2_n)]
            plot_data_E1 += [(R0_i, E1_p + E1_n)]
            plot_data_E2 += [(R0_i, E2_p + E2_n)]
    all_plot_data_E += plot_data_E
    all_plot_data_EE += plot_data_EE
    all_plot_data_E1 += plot_data_E1
    all_plot_data_E2 += plot_data_E2


p = list_plot(all_plot_data_E)
pname = "results/spherical_explosion_all_E_R0.png"
print pname
p.save(pname)

p = list_plot(all_plot_data_EE)
pname = "results/spherical_explosion_all_EE_R0.png"
print pname
p.save(pname)

p = list_plot(all_plot_data_E1)
pname = "results/spherical_explosion_all_E1_R0.png"
print pname
p.save(pname)

p = list_plot(all_plot_data_E2)
pname = "results/spherical_explosion_all_E2_R0.png"
print pname
p.save(pname)