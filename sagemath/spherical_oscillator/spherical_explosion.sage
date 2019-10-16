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


(phi, E1, E2, error) = phi_and_E_lw(q, t, R0, r0, v0, a0, r_min)
print phi, E1, E2, error

step_R0 = 0.5
min_R0 = r0 + step_R0
max_R0 = 5.0

t1 = 0
t2 = 5
dt = 0.5

v0_p = 0
v0_n = 0


a0_p = 0.0001
a0_n = 0.01


all_plot_data = []
for R0_i in np.arange(min_R0, max_R0, step_R0):
    plot_data = []
    for t_i in np.arange(t1, t2, dt):
        (phi_p, E1_p, E2_p, error_p) = phi_and_E_lw(+q, t_i, R0_i, r0, v0_p, a0_p, r_min)
        (phi_n, E1_n, E2_n, error_n) = phi_and_E_lw(-q, t_i, R0_i, r0, v0_n, a0_n, r_min)
        print (phi_p, E1_p, E2_p, error_p)
        print (phi_n, E1_n, E2_n, error_n)
        plot_data += [(t_i, E1_p + E2_p + E1_n + E2_n)]
    all_plot_data += plot_data
    p = list_plot(plot_data)
    pname = "results/spherical_explosion_" + "_R0=" + float_formatting(R0_i) + ".png"
    print pname
    p.save(pname)

p = list_plot(all_plot_data)
pname = "results/spherical_explosion_all_R0.png"
print pname
p.save(pname)
