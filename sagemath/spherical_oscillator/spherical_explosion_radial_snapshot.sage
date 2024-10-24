import sys
try:
    reload(sys)
    sys.setdefaultencoding('utf8')
except:
    pass

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


step_R0 = 0.1
min_R0 = r0 + step_R0
max_R0 = 20.0

t1 = 0
t2 = 100
dt = 0.5
dt_all = 2.0
'''

def spherical_explosion_radial_snapshot(q, t1, t2, dt, r0_p_min, v0_p_min, a0_p_min, r0_p_max, v0_p_max, a0_p_max, r0_n_min, v0_n_min, a0_n_min, r0_n_max, v0_n_max, a0_n_max, step_R0, min_R0, max_R0, r_min, use_dbl_integration = False):
    all_plot_data_phi = []
    all_plot_data_phi_p = []
    all_plot_data_phi_n = []

    all_plot_data_A = []
    all_plot_data_A_p = []
    all_plot_data_A_n = []

    all_plot_data_E = []
    all_plot_data_E_ = []
    all_plot_data_E__p = []
    all_plot_data_E__n = []
    all_plot_data_E_p = []
    all_plot_data_E_n = []

    all_plot_data_E1 = []
    all_plot_data_E1_p = []
    all_plot_data_E1_n = []
    all_plot_data_E2 = []
    all_plot_data_E2_p = []
    all_plot_data_E2_n = []

    all_plot_data_k_p = []
    all_plot_data_k_n = []

    all_plot_data_k1_p = []
    all_plot_data_k1_n = []

    all_plot_data_k2_p = []
    all_plot_data_k2_n = []

    for t_i in np.arange(t1, t2, dt):
        plot_data_phi = []
        plot_data_phi_p = []
        plot_data_phi_n = []

        plot_data_A = []
        plot_data_A_p = []
        plot_data_A_n = []

        plot_data_E_ = []
        plot_data_E = []
        plot_data_E__p = []
        plot_data_E__n = []

        plot_data_E_p = []
        plot_data_E_n = []

        plot_data_E1 = []
        plot_data_E1_p = []
        plot_data_E1_n = []
        plot_data_E2 = []
        plot_data_E2_p = []
        plot_data_E2_n = []
        plot_data_E_E1_E2 = []

        plot_data_k_p = []
        plot_data_k_n = []

        plot_data_k1_p = []
        plot_data_k1_n = []

        plot_data_k2_p = []
        plot_data_k2_n = []

        r_p_min = get_r_of_sphere(+q, t_i, r0_p_min, v0_p_min, a0_p_min, r_min)
        r_p_max = get_r_of_sphere(+q, t_i, r0_p_max, v0_p_max, a0_p_max, r_min)
        r_n_min = get_r_of_sphere(-q, t_i, r0_n_min, v0_n_min, a0_n_min, r_min)
        r_n_max = get_r_of_sphere(-q, t_i, r0_n_max, v0_n_max, a0_n_max, r_min)

        if False:
            (r_p, r_n, capacity, energy) = get_capacity_of_spherical_capacitor(q, t_i, r0, v0_p, a0_p, r0, v0_n, a0_n, r_min)

        r_p_min_met = r_p_min < max_R0 and r_p_min > min_R0
        r_p_max_met = r_p_max < max_R0 and r_p_max > min_R0
        r_n_min_met = r_n_min < max_R0 and r_n_min > min_R0
        r_n_max_met = r_n_max < max_R0 and r_n_max > min_R0

        for R0_i in np.arange(min_R0, max_R0, step_R0):
            if True == use_dbl_integration:
                (phi_p, A_p, E1_p, E2_p, E_p, error_p) = dbl_phi_and_E_lw(+q, t_i, R0_i, r0_p_min, v0_p_min, a0_p_min, r0_p_max, v0_p_max, a0_p_max, r_min)
                (phi_n, A_n, E1_n, E2_n, E_n, error_n) = dbl_phi_and_E_lw(-q, t_i, R0_i, r0_n_min, v0_n_min, a0_n_min, r0_n_max, v0_n_max, a0_n_max, r_min)
            else:
                r0_p = 0.5 * (r0_p_min + r0_p_max)
                r0_n = 0.5 * (r0_n_min + r0_n_max)

                v0_p = 0.5 * (v0_p_min + v0_p_max)
                v0_n = 0.5 * (v0_n_min + v0_n_max)

                a0_p = 0.5 * (a0_p_min + a0_p_max)
                a0_n = 0.5 * (a0_n_min + a0_n_max)

                # (phi_p, A_p, E1_p, E2_p, E_p, error_p) = phi_and_E_lw(+q, t_i, R0_i, r0_p, v0_p, a0_p, r_min)
                # (phi_n, A_n, E1_n, E2_n, E_n, error_n) = phi_and_E_lw(-q, t_i, R0_i, r0_n, v0_n, a0_n, r_min)

                attach("calc_E_with_num_int.sage")
                print ("q_p = ", +q, "t = ", t_i, "R0 = ", R0_i, "r0_p = ", r0_p, "v0_p = ", v0_p, "a0_p = ", a0_p, "r_min =", r_min)
                print ("q_n = ", -q, "t = ", t_i, "R0 = ", R0_i, "r0_n = ", r0_n, "v0_n = ", v0_n, "a0_n = ", a0_n, "r_min =", r_min)
 
                k_p   = k_E     (+q, t_i, R0_i, r0_p, v0_p, a0_p, r_min)
                k_n   = k_E     (-q, t_i, R0_i, r0_n, v0_n, a0_n, r_min)
                k1_p  = k_E1    (+q, t_i, R0_i, r0_p, v0_p, a0_p, r_min)
                k1_n  = k_E1    (-q, t_i, R0_i, r0_n, v0_n, a0_n, r_min)
                k2_p  = k_E2    (+q, t_i, R0_i, r0_p, v0_p, a0_p, r_min)
                k2_n  = k_E2    (-q, t_i, R0_i, r0_n, v0_n, a0_n, r_min)
                phi_p = int_phi (+q, t_i, R0_i, r0_p, v0_p, a0_p, r_min)
                phi_n = int_phi (-q, t_i, R0_i, r0_n, v0_n, a0_n, r_min)
                A_p   = int_A   (+q, t_i, R0_i, r0_p, v0_p, a0_p, r_min)
                A_n   = int_A   (-q, t_i, R0_i, r0_n, v0_n, a0_n, r_min)
                E_p   = int_E   (+q, t_i, R0_i, r0_p, v0_p, a0_p, r_min)
                E_n   = int_E   (-q, t_i, R0_i, r0_n, v0_n, a0_n, r_min)
                E1_p  = int_E1  (+q, t_i, R0_i, r0_p, v0_p, a0_p, r_min)
                E1_n  = int_E1  (-q, t_i, R0_i, r0_n, v0_n, a0_n, r_min)
                E2_p  = int_E2  (+q, t_i, R0_i, r0_p, v0_p, a0_p, r_min)
                E2_n  = int_E2  (-q, t_i, R0_i, r0_n, v0_n, a0_n, r_min)


            print (phi_p, A_p, E1_p, E2_p, E_p, r_p_min, r_p_max)
            print (phi_n, A_n, E1_n, E2_n, E_n, r_n_min, r_n_max)

            plot_data_phi += [(R0_i, phi_p + phi_n)]
            plot_data_phi_p += [(R0_i, phi_p)]
            plot_data_phi_n += [(R0_i, phi_n)]

            plot_data_A += [(R0_i, A_p + A_n)]
            plot_data_A_p += [(R0_i, A_p)]
            plot_data_A_n += [(R0_i, A_n)]

            if r_n_max < R0_i - 50*step_R0:
                plot_data_E_ += [(R0_i, E1_p + E2_p + E1_n + E2_n)]
                plot_data_E += [(R0_i, E_p + E_n)]
                plot_data_E__p += [(R0_i, E1_p + E2_p)]
                plot_data_E__n += [(R0_i, E1_n + E2_n)]
                plot_data_E_p += [(R0_i, E_p)]
                plot_data_E_n += [(R0_i, E_n)]

                plot_data_E1 += [(R0_i, E1_p + E1_n)]
                plot_data_E1_p += [(R0_i, E1_p)]
                plot_data_E1_n += [(R0_i, E1_n)]
                plot_data_E2 += [(R0_i, E2_p + E2_n)]
                plot_data_E2_p += [(R0_i, E2_p)]
                plot_data_E2_n += [(R0_i, E2_n)]
                plot_data_E_E1_E2 += [(R0_i, E1_p + E2_p + E1_n + E2_n), (R0_i, E1_p + E1_n), (R0_i, E2_p + E2_n)]

                plot_data_k_p += [(R0_i, k_p)]
                plot_data_k_n += [(R0_i, k_n)]

                plot_data_k1_p += [(R0_i, k1_p)]
                plot_data_k1_n += [(R0_i, k1_n)]

                plot_data_k2_p += [(R0_i, k2_p)]
                plot_data_k2_n += [(R0_i, k2_n)]


        if True: #abs(t_i/dt_all) > dt/2:
            all_plot_data_phi += plot_data_phi
            all_plot_data_phi_p += plot_data_phi_p
            all_plot_data_phi_n += plot_data_phi_n

            all_plot_data_A += plot_data_A
            all_plot_data_A_p += plot_data_A_p
            all_plot_data_A_n += plot_data_A_n

            all_plot_data_E_ += plot_data_E_
            all_plot_data_E += plot_data_E
            all_plot_data_E__p += plot_data_E__p
            all_plot_data_E__n += plot_data_E__n
            all_plot_data_E_p += plot_data_E_p
            all_plot_data_E_n += plot_data_E_n

            all_plot_data_E1 += plot_data_E1
            all_plot_data_E1_p += plot_data_E1_p
            all_plot_data_E1_n += plot_data_E1_n
            all_plot_data_E2 += plot_data_E2
            all_plot_data_E2_p += plot_data_E2_p
            all_plot_data_E2_n += plot_data_E2_n

            all_plot_data_k_p += plot_data_k_p
            all_plot_data_k_n += plot_data_k_n

            all_plot_data_k1_p += plot_data_k1_p
            all_plot_data_k1_n += plot_data_k1_n

            all_plot_data_k2_p += plot_data_k2_p
            all_plot_data_k2_n += plot_data_k2_n

        dir = os.getcwd()  + "/results_num_int/spherical_explosion_radial_snapshot/"
        if True == use_dbl_integration:
            dir += "dbl_integration_on_volume/"
        dir += "t=" + float_formatting(t_i)
        print ("dir = ", dir)

        try:
            os.makedirs(dir)
        except:
            pass

        folder = dir + "/"

        if len(plot_data_phi) > 0:
            p = list_plot(plot_data_phi)
            pname = folder + "spherical_explosion_phi_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_phi, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

        if len(plot_data_phi_p) > 0:
            p = list_plot(plot_data_phi_p)
            pname = folder + "spherical_explosion_phi_p_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_phi_p, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

        if len(plot_data_phi_n) > 0:
            p = list_plot(plot_data_phi_n)
            pname = folder + "spherical_explosion_phi_n_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_phi_n, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

        if len(plot_data_A) > 0:
            p = list_plot(plot_data_A)
            pname = folder + "spherical_explosion_A_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_A, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

        if len(plot_data_A_p) > 0:
            p = list_plot(plot_data_A_p)
            pname = folder + "spherical_explosion_A_p_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_A_p, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

        if len(plot_data_A_n) > 0:
            p = list_plot(plot_data_A_n)
            pname = folder + "spherical_explosion_A_n_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_A_n, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

        if len(plot_data_E) > 0:
            p = list_plot(plot_data_E)
            pname = folder + "spherical_explosion_E_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

        if len(plot_data_E_p) > 0:
            p = list_plot(plot_data_E_p)
            pname = folder + "spherical_explosion_E_p_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E_p, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

        if len(plot_data_E_n) > 0:
            p = list_plot(plot_data_E_n)
            pname = folder + "spherical_explosion_E_n_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E_n, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

        if len(plot_data_k_p) > 0:
            p = list_plot(plot_data_k_p)
            pname = folder + "spherical_explosion_k_p_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_k_p, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

        if len(plot_data_k_n) > 0:
            p = list_plot(plot_data_k_n)
            pname = folder + "spherical_explosion_k_n_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_k_n, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

        if len(plot_data_k1_p) > 0:
            p = list_plot(plot_data_k1_p)
            pname = folder + "spherical_explosion_k1_p_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_k1_p, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

        if len(plot_data_k1_n) > 0:
            p = list_plot(plot_data_k1_n)
            pname = folder + "spherical_explosion_k1_n_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_k1_n, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

        if len(plot_data_k2_p) > 0:
            p = list_plot(plot_data_k2_p)
            pname = folder + "spherical_explosion_k2_p_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_k2_p, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

        if len(plot_data_k2_n) > 0:
            p = list_plot(plot_data_k2_n)
            pname = folder + "spherical_explosion_k2_n_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_k2_n, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

        if len(plot_data_E1) > 0:
            p = list_plot(plot_data_E1)
            pname = folder + "spherical_explosion_E1_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E1, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

        if len(plot_data_E1_p) > 0:
            p = list_plot(plot_data_E1_p)
            pname = folder + "spherical_explosion_E1_p_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E1_p, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

        if len(plot_data_E1_n) > 0:
            p = list_plot(plot_data_E1_n)
            pname = folder + "spherical_explosion_E1_n_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E1_n, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

        if len(plot_data_E2) > 0:
            p = list_plot(plot_data_E2)
            pname = folder + "spherical_explosion_E2_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E2, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

        if len(plot_data_E2_p) > 0:
            p = list_plot(plot_data_E2_p)
            pname = folder + "spherical_explosion_E2_p_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E2_p, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

        if len(plot_data_E2_n) > 0:
            p = list_plot(plot_data_E2_n)
            pname = folder + "spherical_explosion_E2_n_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E2_n, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

        if len(plot_data_E_E1_E2) > 0:
            p = list_plot(plot_data_E_E1_E2)
            pname = folder + "spherical_explosion_E_E1_E2_R0" + "_t=" + float_formatting(t_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E_E1_E2, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, r_p_min, r_p_max, r_n_min, r_n_max)
            p.save(pname)

    dir = os.getcwd()  + "/results_num_int/spherical_explosion_radial_snapshot/"
    if True == use_dbl_integration:
        dir += "dbl_integration_on_volume/"
    dir += "all_t"
    print ("dir = ", dir)

    try:
        os.makedirs(dir)
    except:
        pass

    folder = dir + "/"

    p = list_plot(all_plot_data_phi)
    pname = folder + "spherical_explosion_all_phi_R0.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_phi_p)
    pname = folder + "spherical_explosion_all_phi_p_R0.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_phi_n)
    pname = folder + "spherical_explosion_all_phi_n_R0.png"
    print (pname)
    p.save(pname)


    p = list_plot(all_plot_data_A)
    pname = folder + "spherical_explosion_all_A_R0.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_A_p)
    pname = folder + "spherical_explosion_all_A_p_R0.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_A_n)
    pname = folder + "spherical_explosion_all_A_n_R0.png"
    print (pname)
    p.save(pname)


    p = list_plot(all_plot_data_E)
    pname = folder + "spherical_explosion_all_E_R0.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_E_p)
    pname = folder + "spherical_explosion_all_E_p_R0.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_E_n)
    pname = folder + "spherical_explosion_all_E_n_R0.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_E1)
    pname = folder + "spherical_explosion_all_E1_R0.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_E1_p)
    pname = folder + "spherical_explosion_all_E1_p_R0.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_E1_n)
    pname = folder + "spherical_explosion_all_E1_n_R0.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_E2)
    pname = folder + "spherical_explosion_all_E2_R0.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_E2_p)
    pname = folder + "spherical_explosion_all_E2_p_R0.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_E2_n)
    pname = "spherical_explosion_all_E2_n_R0.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_k_p)
    pname = folder + "spherical_explosion_all_k_p_R0.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_k_n)
    pname = folder + "spherical_explosion_all_k_n_R0.png"
    print (pname)
    p.save(pname)


    p = list_plot(all_plot_data_k1_p)
    pname = folder + "spherical_explosion_all_k1_p_R0.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_k1_n)
    pname = folder + "spherical_explosion_all_k1_n_R0.png"
    print (pname)
    p.save(pname)


    p = list_plot(all_plot_data_k2_p)
    pname = folder + "spherical_explosion_all_k2_p_R0.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_k2_n)
    pname = folder + "spherical_explosion_all_k2_n_R0.png"
    print (pname)
    p.save(pname)

# spherical_explosion_radial_snapshot(q, t1, t2, dt, r0, v0_p, v0_n, a0_p, a0_n, step_R0, min_R0, max_R0, r_min)