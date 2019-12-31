import sys
reload(sys)
sys.setdefaultencoding('utf8')

import numpy as np

attach("tzap.spyx")
attach("float_formatting.sage")
attach("get_min_max_of_data.sage")
attach("plot_r_of_sphere.sage")


c = get_light_veloncity()
'''
q = 1


r0 = 1

r_min = 0.1

v0_p = 0
v0_n = 0

a0_p = 0.0001
a0_n = 0.01


step_R0 = 1.0
min_R0 = r0 + step_R0
max_R0 = 10.0

t1 = 0
t2 = 10
dt = 0.1
'''

def spherical_explosion_time_evaluation(q, t1, t2, dt, r0, v0_p, v0_n, a0_p, a0_n, step_R0, min_R0, max_R0, r_min):
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
    for R0_i in np.arange(min_R0, max_R0, step_R0):

        r_p_met = False
        r_n_met = False

        plot_data_r_p = []
        plot_data_r_n = []

        plot_data_v_p = []
        plot_data_v_n = []

        t_r_p = NaN
        t_r_n = NaN

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
        plot_data_E_E1_E2 = []
        for t_i in np.arange(t1, t2, dt):
            r_p = get_r_of_sphere(+q, t_i, r0, v0_p, a0_p, r_min)
            r_n = get_r_of_sphere(-q, t_i, r0, v0_n, a0_n, r_min)

            if r_p_met is False and r_p >= R0_i:
                r_p_met = True
                t_r_p = t_i - dt + dt * (R0_i - r_p_pre) / (r_p - r_p_pre)

            if r_n_met is False and r_n >= R0_i:
                r_n_met = True
                t_r_n = t_i - dt + dt * (R0_i - r_n_pre) / (r_n - r_n_pre)

            r_p_pre = r_p
            r_n_pre = r_n

            v_p = get_v_of_sphere(+q, t_i, v0_p, a0_p)
            v_n = get_v_of_sphere(-q, t_i, v0_n, a0_n)

            plot_data_r_p += [(t_i, r_p)]
            plot_data_r_n += [(t_i, r_n)]

            plot_data_v_p += [(t_i, v_p/c)]
            plot_data_v_n += [(t_i, v_n/c)]

            (phi_p, A_p, E1_p, E2_p, error_p) = phi_and_E_lw(+q, t_i, R0_i, r0, v0_p, a0_p, r_min)
            (phi_n, A_n, E1_n, E2_n, error_n) = phi_and_E_lw(-q, t_i, R0_i, r0, v0_n, a0_n, r_min)
            print (phi_p, A_p, E1_p, E2_p, error_p, r_p)
            print (phi_n, A_n, E1_n, E2_n, error_n, r_n)

            if r_n < R0_i - step_R0:
                plot_data_phi += [(t_i, phi_p + phi_n)]
                plot_data_phi_p += [(t_i, phi_p)]
                plot_data_phi_n += [(t_i, phi_n)]

            plot_data_A += [(t_i, A_p + A_n)]
            plot_data_A_p += [(t_i, A_p)]
            plot_data_A_n += [(t_i, A_n)]

            if r_n < R0_i - step_R0:
                plot_data_E += [(t_i, E1_p + E2_p + E1_n + E2_n)]
                plot_data_E_p += [(t_i, E1_p + E2_p)]
                plot_data_E_n += [(t_i, E1_n + E2_n)]

                plot_data_E1 += [(t_i, E1_p + E1_n)]
                plot_data_E1_p += [(t_i, E1_p)]
                plot_data_E1_n += [(t_i, E1_n)]
                plot_data_E2 += [(t_i, E2_p + E2_n)]
                plot_data_E2_p += [(t_i, E2_p)]
                plot_data_E2_n += [(t_i, E2_n)]
                plot_data_E_E1_E2 += [(t_i, E1_p + E2_p + E1_n + E2_n), (t_i, E1_p + E1_n), (t_i, E2_p + E2_n)]

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


        dir = os.getcwd()  + "/results/R0=" + float_formatting(R0_i)
        print "dir = ", dir

        try:
            os.mkdir(dir)
        except:
            pass

        folder = dir + "/"

        if len(plot_data_r_p) > 0:
            p = list_plot(plot_data_r_p)
            pname = folder + "spherical_explosion_r_p_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print pname
            plot_r_of_sphere(plot_data_r_p, p, r_p_met, r_n_met, t_r_p, t_r_n)
            p.save(pname)

        if len(plot_data_r_n) > 0:
            p = list_plot(plot_data_r_n)
            pname = folder + "spherical_explosion_r_n_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print pname
            plot_r_of_sphere(plot_data_r_n, p, r_p_met, r_n_met, t_r_p, t_r_n)
            p.save(pname)

        if len(plot_data_v_p) > 0:
            p = list_plot(plot_data_v_p)
            pname = folder + "spherical_explosion_v_p_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print pname
            plot_r_of_sphere(plot_data_v_p, p, r_p_met, r_n_met, t_r_p, t_r_n)
            p.save(pname)

        if len(plot_data_v_n) > 0:
            p = list_plot(plot_data_v_n)
            pname = folder + "spherical_explosion_v_n_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print pname
            plot_r_of_sphere(plot_data_v_n, p, r_p_met, r_n_met, t_r_p, t_r_n)
            p.save(pname)

        if len(plot_data_phi) > 0:
            p = list_plot(plot_data_phi)
            pname = folder + "spherical_explosion_phi_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print pname
            plot_r_of_sphere(plot_data_phi, p, r_p_met, r_n_met, t_r_p, t_r_n)
            p.save(pname)

        if len(plot_data_phi_p) > 0:
            p = list_plot(plot_data_phi_p)
            pname = folder + "spherical_explosion_phi_p_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print pname
            plot_r_of_sphere(plot_data_phi_p, p, r_p_met, r_n_met, t_r_p, t_r_n)
            p.save(pname)

        if len(plot_data_phi_n) > 0:
            p = list_plot(plot_data_phi_n)
            pname = folder + "spherical_explosion_phi_n_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print pname
            plot_r_of_sphere(plot_data_phi_n, p, r_p_met, r_n_met, t_r_p, t_r_n)
            p.save(pname)

        if len(plot_data_A) > 0:
            p = list_plot(plot_data_A)
            pname = folder + "spherical_explosion_A_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print pname
            plot_r_of_sphere(plot_data_A, p, r_p_met, r_n_met, t_r_p, t_r_n)
            p.save(pname)

        if len(plot_data_A_p) > 0:
            p = list_plot(plot_data_A_p)
            pname = folder + "spherical_explosion_A_p_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print pname
            plot_r_of_sphere(plot_data_A_p, p, r_p_met, r_n_met, t_r_p, t_r_n)
            p.save(pname)

        if len(plot_data_A_n) > 0:
            p = list_plot(plot_data_A_n)
            pname = folder + "spherical_explosion_A_n_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print pname
            plot_r_of_sphere(plot_data_A_n, p, r_p_met, r_n_met, t_r_p, t_r_n)
            p.save(pname)


        if len(plot_data_E) > 0:
            p = list_plot(plot_data_E)
            pname = folder + "spherical_explosion_E_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print pname
            plot_r_of_sphere(plot_data_E, p, r_p_met, r_n_met, t_r_p, t_r_n)
            p.save(pname)

        if len(plot_data_E_p) > 0:
            p = list_plot(plot_data_E_p)
            pname = folder + "spherical_explosion_E_p_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print pname
            plot_r_of_sphere(plot_data_E_p, p, r_p_met, r_n_met, t_r_p, t_r_n)
            p.save(pname)

        if len(plot_data_E_n) > 0:
            p = list_plot(plot_data_E_n)
            pname = folder + "spherical_explosion_E_n_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print pname
            plot_r_of_sphere(plot_data_E_n, p, r_p_met, r_n_met, t_r_p, t_r_n)
            p.save(pname)

        if len(plot_data_E1) > 0:
            p = list_plot(plot_data_E1)
            pname = folder + "spherical_explosion_E1_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print pname
            plot_r_of_sphere(plot_data_E1, p, r_p_met, r_n_met, t_r_p, t_r_n)
            p.save(pname)

        if len(plot_data_E1_p) > 0:
            p = list_plot(plot_data_E1_p)
            pname = folder + "spherical_explosion_E1_p_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print pname
            plot_r_of_sphere(plot_data_E1_p, p, r_p_met, r_n_met, t_r_p, t_r_n)
            p.save(pname)

        if len(plot_data_E1_n) > 0:
            p = list_plot(plot_data_E1_n)
            pname = folder + "spherical_explosion_E1_n_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print pname
            plot_r_of_sphere(plot_data_E1_n, p, r_p_met, r_n_met, t_r_p, t_r_n)
            p.save(pname)

        if len(plot_data_E2) > 0:
            p = list_plot(plot_data_E2)
            pname = folder + "spherical_explosion_E2_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print pname
            plot_r_of_sphere(plot_data_E2, p, r_p_met, r_n_met, t_r_p, t_r_n)
            p.save(pname)

        if len(plot_data_E2_p) > 0:
            p = list_plot(plot_data_E2_p)
            pname = folder + "spherical_explosion_E2_p_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print pname
            plot_r_of_sphere(plot_data_E2_p, p, r_p_met, r_n_met, t_r_p, t_r_n)
            p.save(pname)

        if len(plot_data_E2_n) > 0:
            p = list_plot(plot_data_E2_n)
            pname = folder + "spherical_explosion_E2_n_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print pname
            plot_r_of_sphere(plot_data_E2_n, p, r_p_met, r_n_met, t_r_p, t_r_n)
            p.save(pname)

        if len(plot_data_E_E1_E2) > 0:
            p = list_plot(plot_data_E_E1_E2)
            pname = folder + "spherical_explosion_E_E1_E2_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print pname
            plot_r_of_sphere(plot_data_E_E1_E2, p, r_p_met, r_n_met, t_r_p, t_r_n)
            p.save(pname)

    dir = os.getcwd()  + "/results/all_R0"
    print "dir = ", dir

    try:
        os.mkdir(dir)
    except:
        pass

    folder = dir + "/"

    p = list_plot(all_plot_data_phi)
    pname = folder + "spherical_explosion_all_phi_t.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_phi_p)
    pname = folder + "spherical_explosion_all_phi_p_t.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_phi_n)
    pname = folder + "spherical_explosion_all_phi_n_t.png"
    print pname
    p.save(pname)


    p = list_plot(all_plot_data_A)
    pname = folder + "spherical_explosion_all_A_t.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_A_p)
    pname = folder + "spherical_explosion_all_A_p_t.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_A_n)
    pname = folder + "spherical_explosion_all_A_n_t.png"
    print pname
    p.save(pname)


    p = list_plot(all_plot_data_E)
    pname = folder + "spherical_explosion_all_E_t.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_E_p)
    pname = folder + "spherical_explosion_all_E_p_t.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_E_n)
    pname = folder + "spherical_explosion_all_E_n_t.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_E1)
    pname = folder + "spherical_explosion_all_E1_t.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_E1_p)
    pname = folder + "spherical_explosion_all_E1_p_t.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_E1_n)
    pname = folder + "spherical_explosion_all_E1_n_t.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_E2)
    pname = folder + "spherical_explosion_all_E2_t.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_E2_p)
    pname = folder + "spherical_explosion_all_E2_p_t.png"
    print pname
    p.save(pname)

    p = list_plot(all_plot_data_E2_n)
    pname = folder + "spherical_explosion_all_E2_n_t.png"
    print pname
    p.save(pname)

# spherical_explosion_time_evaluation(q, t1, t2, dt, r0, v0_p, v0_n, a0_p, a0_n, step_R0, min_R0, max_R0, r_min)