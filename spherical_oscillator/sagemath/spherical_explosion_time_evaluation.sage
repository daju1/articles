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


step_R0 = 1.0
min_R0 = r0 + step_R0
max_R0 = 10.0

t1 = 0
t2 = 10
dt = 0.1
'''

def spherical_explosion_time_evaluation(q, t1, t2, dt, r0_p_min, v0_p_min, a0_p_min, r0_p_max, v0_p_max, a0_p_max, r0_n_min, v0_n_min, a0_n_min, r0_n_max, v0_n_max, a0_n_max, step_R0, min_R0, max_R0, r_min, use_dbl_integration = False):
    all_plot_data_phi = []
    all_plot_data_phi_p = []
    all_plot_data_phi_n = []

    all_plot_data_A = []
    all_plot_data_A_p = []
    all_plot_data_A_n = []

    all_plot_data_E_ = []
    all_plot_data_E = []
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
    all_plot_data_dk = []

    all_plot_data_k1_p = []
    all_plot_data_k1_n = []
    all_plot_data_dk1 = []

    all_plot_data_k2_p = []
    all_plot_data_k2_n = []
    all_plot_data_dk2 = []

    for R0_i in np.arange(min_R0, max_R0, step_R0):

        r_p_min_met = False
        r_n_min_met = False

        r_p_max_met = False
        r_n_max_met = False

        plot_data_r_p_min = []
        plot_data_r_n_min = []

        plot_data_r_p_max = []
        plot_data_r_n_max = []

        plot_data_v_p_min = []
        plot_data_v_n_min = []

        plot_data_v_p_max = []
        plot_data_v_n_max = []

        plot_data_capacity = []
        plot_data_energy = []

        t_r_p_min = NaN
        t_r_p_max = NaN
        t_r_n_min = NaN
        t_r_n_max = NaN

        plot_data_phi = []
        plot_data_phi_p = []
        plot_data_phi_n = []

        plot_data_A = []
        plot_data_A_p = []
        plot_data_A_n = []

        plot_data_E12 = []
        plot_data_E = []
        plot_data_E12_p = []
        plot_data_E12_n = []
        plot_data_E_p = []
        plot_data_E_n = []

        plot_data_E1 = []
        plot_data_E1_p = []
        plot_data_E1_n = []
        plot_data_E2 = []
        plot_data_E2_p = []
        plot_data_E2_n = []
        plot_data_E_E1_E2 = []
        plot_data_E12_E1_E2 = []

        plot_data_k_p = []
        plot_data_k_n = []
        plot_data_k12_p = []
        plot_data_k12_n = []
        plot_data_dk = []

        plot_data_k1_p = []
        plot_data_k1_n = []
        plot_data_dk1 = []

        plot_data_k2_p = []
        plot_data_k2_n = []
        plot_data_dk2 = []

        for t_i in np.arange(t1, t2, dt):
            r_p_min = get_r_of_sphere(+q, t_i, r0_p_min, v0_p_min, a0_p_min, r_min)
            r_p_max = get_r_of_sphere(+q, t_i, r0_p_max, v0_p_max, a0_p_max, r_min)
            r_n_min = get_r_of_sphere(-q, t_i, r0_n_min, v0_n_min, a0_n_min, r_min)
            r_n_max = get_r_of_sphere(-q, t_i, r0_n_max, v0_n_max, a0_n_max, r_min)

            if False:
                (r_p, r_n, capacity, energy) = get_capacity_of_spherical_capacitor(q, t_i, r0, v0_p, a0_p, r0, v0_n, a0_n, r_min)

            if r_p_min_met is False and r_p_min >= R0_i:
                r_p_min_met = True
                t_r_p_min = t_i - dt + dt * (R0_i - r_p_min_pre) / (r_p_min - r_p_min_pre)

            if r_p_max_met is False and r_p_max >= R0_i:
                r_p_max_met = True
                t_r_p_max = t_i - dt + dt * (R0_i - r_p_max_pre) / (r_p_max - r_p_max_pre)

            if r_n_min_met is False and r_n_min >= R0_i:
                r_n_min_met = True
                t_r_n_min = t_i - dt + dt * (R0_i - r_n_min_pre) / (r_n_min - r_n_min_pre)

            if r_n_max_met is False and r_n_max >= R0_i:
                r_n_max_met = True
                t_r_n_max = t_i - dt + dt * (R0_i - r_n_max_pre) / (r_n_max - r_n_max_pre)

            r_p_min_pre = r_p_min
            r_n_min_pre = r_n_min

            r_p_max_pre = r_p_max
            r_n_max_pre = r_n_max

            v_p_min = get_v_of_sphere(+q, t_i, v0_p_min, a0_p_min)
            v_n_min = get_v_of_sphere(-q, t_i, v0_n_min, a0_n_min)

            v_p_max = get_v_of_sphere(+q, t_i, v0_p_max, a0_p_max)
            v_n_max = get_v_of_sphere(-q, t_i, v0_n_max, a0_n_max)

            plot_data_r_p_min += [(t_i, r_p_min)]
            plot_data_r_n_min += [(t_i, r_n_min)]

            plot_data_r_p_max += [(t_i, r_p_max)]
            plot_data_r_n_max += [(t_i, r_n_max)]

            #if r_p != r_n:
            #    plot_data_capacity += [(t_i, capacity)]
            #    plot_data_energy += [(t_i, energy)]

            plot_data_v_p_min += [(t_i, v_p_min/c)]
            plot_data_v_n_min += [(t_i, v_n_min/c)]

            plot_data_v_p_max += [(t_i, v_p_max/c)]
            plot_data_v_n_max += [(t_i, v_n_max/c)]

            if True == use_dbl_integration:
                (phi_p, A_p, E1_p, E2_p, E_p, error_p) = dbl_phi_and_E_lw(+q, t_i, R0_i, r0_p_min, v0_p_min, a0_p_min, r0_p_max, v0_p_max, a0_p_max, r_min)
                (phi_n, A_n, E1_n, E2_n, E_n, error_n) = dbl_phi_and_E_lw(-q, t_i, R0_i, r0_n_min, v0_n_min, a0_n_min, r0_n_max, v0_n_max, a0_n_max, r_min)
            else:
                # (phi_p, A_p, E1_p, E2_p, error_p) = phi_and_E_lw(+q, t_i, R0_i, 0.5 * (r0_p_min + r0_p_max), 0.5 * (v0_p_min + v0_p_max), 0.5 * (a0_p_min + a0_p_max), r_min)
                # (phi_n, A_n, E1_n, E2_n, error_n) = phi_and_E_lw(-q, t_i, R0_i, 0.5 * (r0_n_min + r0_n_max), 0.5 * (v0_n_min + v0_n_max), 0.5 * (a0_n_min + a0_n_max), r_min)
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

            if r_n_max < R0_i - step_R0:
                plot_data_phi += [(t_i, phi_p + phi_n)]
                plot_data_phi_p += [(t_i, phi_p)]
                plot_data_phi_n += [(t_i, phi_n)]

            plot_data_A += [(t_i, A_p + A_n)]
            plot_data_A_p += [(t_i, A_p)]
            plot_data_A_n += [(t_i, A_n)]

            if r_n_max < R0_i - step_R0:
                plot_data_E12 += [(t_i, E1_p + E2_p + E1_n + E2_n)]
                plot_data_E += [(t_i, E_p + E_n)]
                plot_data_E12_p += [(t_i, E1_p + E2_p)]
                plot_data_E12_n += [(t_i, E1_n + E2_n)]
                plot_data_E_p += [(t_i, E_p)]
                plot_data_E_n += [(t_i, E_n)]

                plot_data_E1 += [(t_i, E1_p + E1_n)]
                plot_data_E1_p += [(t_i, E1_p)]
                plot_data_E1_n += [(t_i, E1_n)]
                plot_data_E2 += [(t_i, E2_p + E2_n)]
                plot_data_E2_p += [(t_i, E2_p)]
                plot_data_E2_n += [(t_i, E2_n)]
                plot_data_E_E1_E2 += [(t_i, E_p + E_n), (t_i, E1_p + E1_n), (t_i, E2_p + E2_n)]
                plot_data_E12_E1_E2 += [(t_i, E1_p + E2_p + E1_n + E2_n), (t_i, E1_p + E1_n), (t_i, E2_p + E2_n)]

                plot_data_k_p += [(t_i, k_p)]
                plot_data_k_n += [(t_i, k_n)]
                plot_data_k12_p += [(t_i, k1_p + k2_p)]
                plot_data_k12_n += [(t_i, k1_n + k2_n)]

                plot_data_dk += [(t_i, k_p - k_n)]
                plot_data_k1_p += [(t_i, k1_p)]
                plot_data_k1_n += [(t_i, k1_n)]
                plot_data_dk1 += [(t_i, k1_p - k1_n)]
                plot_data_k2_p += [(t_i, k2_p)]
                plot_data_k2_n += [(t_i, k2_n)]
                plot_data_dk2 += [(t_i, k2_p - k2_n)]

        all_plot_data_phi += plot_data_phi
        all_plot_data_phi_p += plot_data_phi_p
        all_plot_data_phi_n += plot_data_phi_n

        all_plot_data_A += plot_data_A
        all_plot_data_A_p += plot_data_A_p
        all_plot_data_A_n += plot_data_A_n

        all_plot_data_E_ += plot_data_E12
        all_plot_data_E += plot_data_E
        all_plot_data_E__p += plot_data_E12_p
        all_plot_data_E__n += plot_data_E12_n
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
        all_plot_data_dk += plot_data_dk

        all_plot_data_k1_p += plot_data_k1_p
        all_plot_data_k1_n += plot_data_k1_n
        all_plot_data_dk1 += plot_data_dk1

        all_plot_data_k2_p += plot_data_k2_p
        all_plot_data_k2_n += plot_data_k2_n
        all_plot_data_dk2 += plot_data_dk2


        print("r_p_min_met = ", r_p_min_met, "r_p_max_met = ", r_p_max_met, "r_n_min_met = ", r_n_min_met, "r_n_max_met = ", r_n_max_met, "t_r_p_min = ", t_r_p_min, "t_r_p_max = ", t_r_p_max, "t_r_n_min = ", t_r_n_min, "t_r_n_max = ", t_r_n_max)

        dir = os.getcwd()  + "/results_num_int/spherical_explosion_time_evaluation/"
        if True == use_dbl_integration:
            dir += "dbl_integration_on_volume/"
        dir += "R0=" + float_formatting(R0_i)
        print ("dir = ", dir)

        try:
            os.makedirs(dir)
        except:
            pass

        folder = dir + "/"

        if len(plot_data_capacity) > 0:
            p = list_plot(plot_data_capacity)
            pname = folder + "spherical_explosion_capacity_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_capacity, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_energy) > 0:
            p = list_plot(plot_data_energy)
            pname = folder + "spherical_explosion_energy_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_energy, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_r_p_min) > 0:
            p = list_plot(plot_data_r_p_min)
            pname = folder + "spherical_explosion_r_p_min_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_r_p_min, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_r_p_max) > 0:
            p = list_plot(plot_data_r_p_max)
            pname = folder + "spherical_explosion_r_p_max_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_r_p_max, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_r_n_min) > 0:
            p = list_plot(plot_data_r_n_min)
            pname = folder + "spherical_explosion_r_n_min_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_r_n_min, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_r_n_max) > 0:
            p = list_plot(plot_data_r_n_max)
            pname = folder + "spherical_explosion_r_n_max_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_r_n_max, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_v_p_min) > 0:
            p = list_plot(plot_data_v_p_min)
            pname = folder + "spherical_explosion_v_p_min_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_v_p_min, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_v_p_max) > 0:
            p = list_plot(plot_data_v_p_max)
            pname = folder + "spherical_explosion_v_p_max_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_v_p_max, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_v_n_min) > 0:
            p = list_plot(plot_data_v_n_min)
            pname = folder + "spherical_explosion_v_n_min_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_v_n_min, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_v_n_max) > 0:
            p = list_plot(plot_data_v_n_max)
            pname = folder + "spherical_explosion_v_n_max_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_v_n_max, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_phi) > 0:
            p = list_plot(plot_data_phi)
            pname = folder + "spherical_explosion_phi_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_phi, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_phi_p) > 0:
            p = list_plot(plot_data_phi_p)
            pname = folder + "spherical_explosion_phi_p_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_phi_p, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_phi_n) > 0:
            p = list_plot(plot_data_phi_n)
            pname = folder + "spherical_explosion_phi_n_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_phi_n, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_A) > 0:
            p = list_plot(plot_data_A)
            pname = folder + "spherical_explosion_A_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_A, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_A_p) > 0:
            p = list_plot(plot_data_A_p)
            pname = folder + "spherical_explosion_A_p_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_A_p, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_A_n) > 0:
            p = list_plot(plot_data_A_n)
            pname = folder + "spherical_explosion_A_n_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_A_n, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_E12) > 0:
            p = list_plot(plot_data_E12)
            pname = folder + "spherical_explosion_E12_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E12, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname, title="$E_1 + E_2 ( R_0 = " + float_formatting(R0_i) + ")$", axes_labels=['$t$', '$E_1 + E_2$'])

        if len(plot_data_E12_p) > 0:
            p = list_plot(plot_data_E12_p)
            pname = folder + "spherical_explosion_E12_p_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E12_p, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_E12_n) > 0:
            p = list_plot(plot_data_E12_n)
            pname = folder + "spherical_explosion_E12_n_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E12_n, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_E) > 0:
            p = list_plot(plot_data_E)
            pname = folder + "spherical_explosion_E_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname, title="$E (R_0 = " + float_formatting(R0_i) + ")$", axes_labels=['$t$', '$E$'])

        if len(plot_data_E_p) > 0:
            p = list_plot(plot_data_E_p)
            pname = folder + "spherical_explosion_E_p_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E_p, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_E_n) > 0:
            p = list_plot(plot_data_E_n)
            pname = folder + "spherical_explosion_E_n_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E_n, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_k12_p) > 0:
            p = list_plot(plot_data_k12_p)
            pname = folder + "spherical_explosion_k12_p_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_k12_p, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname, title="$Gauss theorem multiplier for positive charge (" + str(min(plot_data_k12_p)) + "..." + str(max(plot_data_k12_p)) + ")$", axes_labels=['$t$', '$k_p$'])

        if len(plot_data_k12_n) > 0:
            p = list_plot(plot_data_k12_n)
            pname = folder + "spherical_explosion_k12_n_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_k12_n, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname, title="$Gauss theorem multiplier for negative charge (" + str(min(plot_data_k12_n)) + "..." + str(max(plot_data_k12_n)) + ")$", axes_labels=['$t$', '$k_n$'])


        if len(plot_data_k_p) > 0:
            p = list_plot(plot_data_k_p)
            pname = folder + "spherical_explosion_k_p_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_k_p, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname, title="$Gauss theorem multiplier for positive charge (" + str(min(plot_data_k_p)[1]) + "..." + str(max(plot_data_k_p)[1]) + ")$", axes_labels=['$t$', '$k_p$'])

        if len(plot_data_k_n) > 0:
            p = list_plot(plot_data_k_n)
            pname = folder + "spherical_explosion_k_n_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_k_n, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname, title="$Gauss theorem multiplier for negative charge (" + str(min(plot_data_k_n)[1]) + "..." + str(max(plot_data_k_n)[1]) + ")$", axes_labels=['$t$', '$k_n$'])

        if len(plot_data_dk) > 0:
            p = list_plot(plot_data_dk)
            pname = folder + "spherical_explosion_dk_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_dk, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_k1_p) > 0:
            p = list_plot(plot_data_k1_p)
            pname = folder + "spherical_explosion_k1_p_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_k1_p, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_k1_n) > 0:
            p = list_plot(plot_data_k1_n)
            pname = folder + "spherical_explosion_k1_n_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_k1_n, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_dk1) > 0:
            p = list_plot(plot_data_dk1)
            pname = folder + "spherical_explosion_dk1_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_dk1, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_k2_p) > 0:
            p = list_plot(plot_data_k2_p)
            pname = folder + "spherical_explosion_k2_p_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_k2_p, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_k2_n) > 0:
            p = list_plot(plot_data_k2_n)
            pname = folder + "spherical_explosion_k2_n_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_k2_n, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_dk2) > 0:
            p = list_plot(plot_data_dk2)
            pname = folder + "spherical_explosion_dk2_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_dk2, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_E1) > 0:
            p = list_plot(plot_data_E1)
            pname = folder + "spherical_explosion_E1_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E1, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname, title="$E_1 (R_0 = " + float_formatting(R0_i) + ")$", axes_labels=['$t$', '$E_1$'])

        if len(plot_data_E1_p) > 0:
            p = list_plot(plot_data_E1_p)
            pname = folder + "spherical_explosion_E1_p_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E1_p, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_E1_n) > 0:
            p = list_plot(plot_data_E1_n)
            pname = folder + "spherical_explosion_E1_n_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E1_n, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_E2) > 0:
            p = list_plot(plot_data_E2)
            pname = folder + "spherical_explosion_E2_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E2, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname, title="$E_2 (R_0 = " + float_formatting(R0_i) + ")$", axes_labels=['$t$', '$E_2$'])

        if len(plot_data_E2_p) > 0:
            p = list_plot(plot_data_E2_p)
            pname = folder + "spherical_explosion_E2_p_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E2_p, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_E2_n) > 0:
            p = list_plot(plot_data_E2_n)
            pname = folder + "spherical_explosion_E2_n_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E2_n, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_E12_E1_E2) > 0:
            p = list_plot(plot_data_E12_E1_E2)
            pname = folder + "spherical_explosion_E12_E1_E2_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E12_E1_E2, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname)

        if len(plot_data_E_E1_E2) > 0:
            p = list_plot(plot_data_E_E1_E2)
            pname = folder + "spherical_explosion_E_E1_E2_t" + "_R0=" + float_formatting(R0_i) + ".png"
            print (pname)
            plot_r_of_sphere(plot_data_E_E1_E2, p, r_p_min_met, r_p_max_met, r_n_min_met, r_n_max_met, t_r_p_min, t_r_p_max, t_r_n_min, t_r_n_max)
            p.save(pname, title="$E, E_1, E_2 (R_0 = " + float_formatting(R0_i) + ")$", axes_labels=['$t$', '$E, E_1, E_2$'])

    dir = os.getcwd()  + "/results_num_int/spherical_explosion_time_evaluation/"
    if True == use_dbl_integration:
        dir += "dbl_integration_on_volume/"
    dir +="all_R0"
    print ("dir = ", dir)

    try:
        os.makedirs(dir)
    except:
        pass

    folder = dir + "/"

    p = list_plot(all_plot_data_phi)
    pname = folder + "spherical_explosion_all_phi_t.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_phi_p)
    pname = folder + "spherical_explosion_all_phi_p_t.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_phi_n)
    pname = folder + "spherical_explosion_all_phi_n_t.png"
    print (pname)
    p.save(pname)


    p = list_plot(all_plot_data_A)
    pname = folder + "spherical_explosion_all_A_t.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_A_p)
    pname = folder + "spherical_explosion_all_A_p_t.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_A_n)
    pname = folder + "spherical_explosion_all_A_n_t.png"
    print (pname)
    p.save(pname)


    p = list_plot(all_plot_data_E_)
    pname = folder + "spherical_explosion_all_E__t.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_E)
    pname = folder + "spherical_explosion_all_E_t.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_E_p)
    pname = folder + "spherical_explosion_all_E_p_t.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_E_n)
    pname = folder + "spherical_explosion_all_E_n_t.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_E1)
    pname = folder + "spherical_explosion_all_E1_t.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_E1_p)
    pname = folder + "spherical_explosion_all_E1_p_t.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_E1_n)
    pname = folder + "spherical_explosion_all_E1_n_t.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_E2)
    pname = folder + "spherical_explosion_all_E2_t.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_E2_p)
    pname = folder + "spherical_explosion_all_E2_p_t.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_E2_n)
    pname = folder + "spherical_explosion_all_E2_n_t.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_k_p)
    pname = folder + "spherical_explosion_all_k_p_t.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_k_n)
    pname = folder + "spherical_explosion_all_k_n_t.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_dk)
    pname = folder + "spherical_explosion_all_dk_t.png"
    print (pname)
    p.save(pname)


    p = list_plot(all_plot_data_k1_p)
    pname = folder + "spherical_explosion_all_k1_p_t.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_k1_n)
    pname = folder + "spherical_explosion_all_k1_n_t.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_dk1)
    pname = folder + "spherical_explosion_all_dk1_t.png"
    print (pname)
    p.save(pname)


    p = list_plot(all_plot_data_k2_p)
    pname = folder + "spherical_explosion_all_k2_p_t.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_k2_n)
    pname = folder + "spherical_explosion_all_k2_n_t.png"
    print (pname)
    p.save(pname)

    p = list_plot(all_plot_data_dk2)
    pname = folder + "spherical_explosion_all_dk2_t.png"
    print (pname)
    p.save(pname)

