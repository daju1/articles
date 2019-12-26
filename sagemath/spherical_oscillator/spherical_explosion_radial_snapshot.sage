import sys
reload(sys)
sys.setdefaultencoding('utf8')

import numpy as np

attach("tzap.spyx")
attach("float_formatting.sage")
attach("get_min_max_of_data.sage")

from sage.plot.line import Line

c = get_light_veloncity()
'''
q = 1

r0 = 1

r_min = 0.1

v0_p = 0
v0_n = 0

a0_p = 0.0001
a0_n = 0.01


step_R0 = 0.1
min_R0 = r0 + step_R0
max_R0 = 20.0

t1 = 0
t2 = 100
dt = 0.5
dt_all = 2.0
'''

def spherical_explosion_radial_snapshot(q, t1, t2, dt, r0, v0_p, v0_n, a0_p, a0_n, step_R0, min_R0, max_R0, r_min):
    all_plot_data_phi = []
    all_plot_data_phi_p = []
    all_plot_data_phi_n = []

    all_plot_data_A = []
    all_plot_data_A_p = []
    all_plot_data_A_n = []

    all_plot_data_E = []
    all_plot_data_E_p = []
    all_plot_data_E_n = []

    all_plot_data_E1 = []
    all_plot_data_E1_p = []
    all_plot_data_E1_n = []
    all_plot_data_E2 = []
    all_plot_data_E2_p = []
    all_plot_data_E2_n = []
    for t_i in np.arange(t1, t2, dt):
        plot_data_phi = []
        plot_data_phi_p = []
        plot_data_phi_n = []

        plot_data_A = []
        plot_data_A_p = []
        plot_data_A_n = []

        plot_data_E = []
        plot_data_E_p = []
        plot_data_E_n = []

        plot_data_E1 = []
        plot_data_E1_p = []
        plot_data_E1_n = []
        plot_data_E2 = []
        plot_data_E2_p = []
        plot_data_E2_n = []

        r_p = get_r_of_sphere(+q, t_i, r0, v0_p, a0_p, r_min)
        r_n = get_r_of_sphere(-q, t_i, r0, v0_n, a0_n, r_min)

        for R0_i in np.arange(min_R0, max_R0, step_R0):
            (phi_p, A_p, E1_p, E2_p, error_p) = phi_and_E_lw(+q, t_i, R0_i, r0, v0_p, a0_p, r_min)
            (phi_n, A_n, E1_n, E2_n, error_n) = phi_and_E_lw(-q, t_i, R0_i, r0, v0_n, a0_n, r_min)
            print (phi_p, A_p, E1_p, E2_p, error_p, r_p)
            print (phi_n, A_n, E1_n, E2_n, error_n, r_n)

            plot_data_phi += [(R0_i, phi_p + phi_n)]
            plot_data_phi_p += [(R0_i, phi_p)]
            plot_data_phi_n += [(R0_i, phi_n)]

            plot_data_A += [(R0_i, A_p + A_n)]
            plot_data_A_p += [(R0_i, A_p)]
            plot_data_A_n += [(R0_i, A_n)]

            if r_n < R0_i - step_R0:
                plot_data_E += [(R0_i, E1_p + E2_p + E1_n + E2_n)]
                plot_data_E_p += [(R0_i, E1_p + E2_p)]
                plot_data_E_n += [(R0_i, E1_n + E2_n)]

                plot_data_E1 += [(R0_i, E1_p + E1_n)]
                plot_data_E1_p += [(R0_i, E1_p)]
                plot_data_E1_n += [(R0_i, E1_n)]
                plot_data_E2 += [(R0_i, E2_p + E2_n)]
                plot_data_E2_p += [(R0_i, E2_p)]
                plot_data_E2_n += [(R0_i, E2_n)]
        if abs(t_i/dt_all) > dt/2:
            all_plot_data_phi += plot_data_phi
            all_plot_data_phi_p += plot_data_phi_p
            all_plot_data_phi_n += plot_data_phi_n

            all_plot_data_A += plot_data_A
            all_plot_data_A_p += plot_data_A_p
            all_plot_data_A_n += plot_data_A_n

            all_plot_data_E += plot_data_E
            all_plot_data_E_p += plot_data_E_p
            all_plot_data_E_n += plot_data_E_n

            all_plot_data_E1 += plot_data_E1
            all_plot_data_E1_p += plot_data_E1_p
            all_plot_data_E1_n += plot_data_E1_n
            all_plot_data_E2 += plot_data_E2
            all_plot_data_E2_p += plot_data_E2_p
            all_plot_data_E2_n += plot_data_E2_n

        if len(plot_data_phi) > 0:
            p = list_plot(plot_data_phi)
            pname = "results/spherical_explosion_phi_R0" + "_t=" + float_formatting(t_i) + ".png"
            print pname
            (min_phi, max_phi) = get_min_max_of_data(plot_data_phi)
            L = Line([r_p,r_p], [min_phi, max_phi],{'alpha':1,'thickness':1,'rgbcolor':(1,0,0),'legend_label':''})
            p.add_primitive(L)
            L = Line([r_n,r_n], [min_phi, max_phi],{'alpha':1,'thickness':1,'rgbcolor':(0,0,1),'legend_label':''})
            p.add_primitive(L)
            p.save(pname)

        if len(plot_data_phi_p) > 0:
            p = list_plot(plot_data_phi_p)
            pname = "results/spherical_explosion_phi_p_R0" + "_t=" + float_formatting(t_i) + ".png"
            print pname
            p.save(pname)

        if len(plot_data_phi_n) > 0:
            p = list_plot(plot_data_phi_n)
            pname = "results/spherical_explosion_phi_n_R0" + "_t=" + float_formatting(t_i) + ".png"
            print pname
            p.save(pname)

        if len(plot_data_A) > 0:
            p = list_plot(plot_data_A)
            pname = "results/spherical_explosion_A_R0" + "_t=" + float_formatting(t_i) + ".png"
            print pname
            p.save(pname)

        if len(plot_data_A_p) > 0:
            p = list_plot(plot_data_A_p)
            pname = "results/spherical_explosion_A_p_R0" + "_t=" + float_formatting(t_i) + ".png"
            print pname
            p.save(pname)

        if len(plot_data_A_n) > 0:
            p = list_plot(plot_data_A_n)
            pname = "results/spherical_explosion_A_n_R0" + "_t=" + float_formatting(t_i) + ".png"
            print pname
            p.save(pname)

        if len(plot_data_E) > 0:
            p = list_plot(plot_data_E)
            pname = "results/spherical_explosion_E_R0" + "_t=" + float_formatting(t_i) + ".png"
            print pname
            p.save(pname)

        if len(plot_data_E_p) > 0:
            p = list_plot(plot_data_E_p)
            pname = "results/spherical_explosion_E_p_R0" + "_t=" + float_formatting(t_i) + ".png"
            print pname
            p.save(pname)

        if len(plot_data_E_n) > 0:
            p = list_plot(plot_data_E_n)
            pname = "results/spherical_explosion_E_n_R0" + "_t=" + float_formatting(t_i) + ".png"
            print pname
            p.save(pname)

        if len(plot_data_E1) > 0:
            p = list_plot(plot_data_E1)
            pname = "results/spherical_explosion_E1_R0" + "_t=" + float_formatting(t_i) + ".png"
            print pname
            p.save(pname)

        if len(plot_data_E1_p) > 0:
            p = list_plot(plot_data_E1_p)
            pname = "results/spherical_explosion_E1_p_R0" + "_t=" + float_formatting(t_i) + ".png"
            print pname
            p.save(pname)

        if len(plot_data_E1_n) > 0:
            p = list_plot(plot_data_E1_n)
            pname = "results/spherical_explosion_E1_n_R0" + "_t=" + float_formatting(t_i) + ".png"
            print pname
            p.save(pname)

        if len(plot_data_E2) > 0:
            p = list_plot(plot_data_E2)
            pname = "results/spherical_explosion_E2_R0" + "_t=" + float_formatting(t_i) + ".png"
            print pname
            p.save(pname)

        if len(plot_data_E2_p) > 0:
            p = list_plot(plot_data_E2_p)
            pname = "results/spherical_explosion_E2_p_R0" + "_t=" + float_formatting(t_i) + ".png"
            print pname
            p.save(pname)

        if len(plot_data_E2_n) > 0:
            p = list_plot(plot_data_E2_n)
            pname = "results/spherical_explosion_E2_n_R0" + "_t=" + float_formatting(t_i) + ".png"
            print pname
            p.save(pname)


    p = list_plot(all_plot_data_phi)
    pname = "results/spherical_explosion_all_phi_R0.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_phi_p)
    pname = "results/spherical_explosion_all_phi_p_R0.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_phi_n)
    pname = "results/spherical_explosion_all_phi_n_R0.png"
    print pname
    p.save(pname)


    p = list_plot(all_plot_data_A)
    pname = "results/spherical_explosion_all_A_R0.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_A_p)
    pname = "results/spherical_explosion_all_A_p_R0.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_A_n)
    pname = "results/spherical_explosion_all_A_n_R0.png"
    print pname
    p.save(pname)


    p = list_plot(all_plot_data_E)
    pname = "results/spherical_explosion_all_E_R0.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_E_p)
    pname = "results/spherical_explosion_all_E_p_R0.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_E_n)
    pname = "results/spherical_explosion_all_E_n_R0.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_E1)
    pname = "results/spherical_explosion_all_E1_R0.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_E1_p)
    pname = "results/spherical_explosion_all_E1_p_R0.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_E1_n)
    pname = "results/spherical_explosion_all_E1_n_R0.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_E2)
    pname = "results/spherical_explosion_all_E2_R0.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_E2_p)
    pname = "results/spherical_explosion_all_E2_p_R0.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_E2_n)
    pname = "results/spherical_explosion_all_E2_n_R0.png"
    print pname
    p.save(pname)

# spherical_explosion_radial_snapshot(q, t1, t2, dt, r0, v0_p, v0_n, a0_p, a0_n, step_R0, min_R0, max_R0, r_min)