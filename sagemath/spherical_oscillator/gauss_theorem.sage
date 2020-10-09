import sys
try:
    reload(sys)
    sys.setdefaultencoding('utf8')
except:
    pass
import numpy as np

attach("tzap.spyx")
attach("float_formatting.sage")
c = get_light_veloncity()
m = get_electric_field_multiplier()


def get_folder(t, r0, v0, a0, R0):
    dir = os.getcwd()  + "/results/gauss_theorem/"
    dir += suffix_trvaR(t, r0, v0, a0, R0)
    print ("dir = ", dir)

    try:
        os.makedirs(dir)
    except:
        pass

    folder = dir + "/"
    return folder


def plot_gauss_theorem_theta_E(q, trange, R0, r0, v0, a0, r_min):
    npoints = 180
    thetas = np.arange(0*pi/npoints, npoints*pi/npoints + 1*pi/npoints, 1*pi/npoints)
    data_tzap = []
    data_rzap = []
    data_E = []
    data_E_sin = []
    for t in trange:
        for theta in thetas:
            (t2, r1) = tzap(q, t, R0, r0, v0, a0, theta, r_min)
            (phi, A, E1, E2, E, error) = phi_and_E_integrand(q, t, R0, r0, v0, a0, theta, r_min)
            data_E += [ (theta, E/m) ]
            data_E_sin += [ (theta, E/m*sin(theta)) ]
            data_tzap += [ (theta, t2) ]
            data_rzap += [ (theta, r1) ]

    p = list_plot (data_E)
    p.save(folder + "gauss_theorem_theta_E" + suffix(t, r0, a0) + ".png")

    p = list_plot (data_E_sin)
    p.save(folder + "gauss_theorem_theta_E_sin" + suffix(t, r0, a0) + ".png")

    p = list_plot (data_tzap)
    p.save(folder + "gauss_theorem_theta_tzap" + suffix(t, r0, a0) + ".png")

    p = list_plot (data_rzap)
    p.save(folder + "gauss_theorem_theta_rzap" + suffix(t, r0, a0) + ".png")

def plot_gauss_theorem_t_E(q, trange, R0, r0, v0, a0, r_min):
    attach("calc_E_with_num_int.sage")
    data_k = []
    data_E = []
    data_int_E = []
    data_int_E1 = []
    data_int_E2 = []
    data_E_int_E = []
    data_E1_int_E1 = []
    data_E2_int_E2 = []
    for t in trange:
        k   = k_E   (q, t, R0, r0, v0, a0, r_min)
        iE  = int_E (q, t, R0, r0, v0, a0, r_min)
        iE1 = int_E1(q, t, R0, r0, v0, a0, r_min)
        iE2 = int_E2(q, t, R0, r0, v0, a0, r_min)
        (phi, A, E1, E2, E, error) = phi_and_E_lw(q, t, R0, r0, v0, a0, r_min)
        print (phi, A, E1, E2, E, error)
        data_k += [(t, k)]
        data_E += [(t, E)]

        data_int_E += [(t, iE)]

        data_E_int_E += [(t, E)]
        data_E_int_E += [(t, iE)]

        data_E1_int_E1 += [(t, E1)]
        data_E1_int_E1 += [(t, iE1)]

        data_E2_int_E2 += [(t, E2)]
        data_E2_int_E2 += [(t, iE2)]

    p = list_plot (data_k)
    p.save(folder + "gauss_theorem_k" + suffix(t2, r0, a0) + ".png")

    p = list_plot (data_E)
    p.save(folder + "gauss_theorem_E" + suffix(t2, r0, a0) + ".png")

    p = list_plot (data_int_E)
    p.save(folder + "gauss_theorem_int_E" + suffix(t2, r0, a0) + ".png")

    p = list_plot (data_E_int_E)
    p.save(folder + "gauss_theorem_E_int_E" + suffix(t2, r0, a0) + ".png")

    p = list_plot (data_E1_int_E1)
    p.save(folder + "gauss_theorem_E1_int_E1" + suffix(t2, r0, a0) + ".png")

    p = list_plot (data_E2_int_E2)
    p.save(folder + "gauss_theorem_E2_int_E2" + suffix(t2, r0, a0) + ".png")

'''
q = 1
R0 = 1 * c
r0 = 0.1 * c
v0 = 0
a0 = 0.05 * c
r_min = 0.01
t1 = 0
t2 = 5
dt = 0.1
trange = np.arange(t1, t2, dt)

folder = get_folder(t2, r0, v0, a0, R0)
plot_gauss_theorem_theta_E(q, trange, R0, r0, v0, a0, r_min)
plot_gauss_theorem_t_E    (q, trange, R0, r0, v0, a0, r_min)

'''


# these are half values from results of copper_explosion_lw.sage
q_p =  2.12744210611535
r0_p =  0.000277815949170315 
v0_p =  0.000000000000000 
a0_p =  2.77810674609282e6 
r_min = 0.0000555631898340631

q_n =  -2.12744210611535 
r0_n =  0.000277815949170315 
v0_n =  0.000000000000000 
a0_n =  9.48896834491850e8 
r_min = 0.0000555631898340631

R0 =  0.15

t1 = 0
Delta_t =  0.000347500236321034
t2 = Delta_t
dt = 0.000001

trange = np.arange(t1, t2, dt)
'''
folder = get_folder(t2, r0_p, v0_p, a0_p, R0)
plot_gauss_theorem_theta_E(q_p, trange, R0, r0_p, v0_p, a0_p, r_min)
plot_gauss_theorem_t_E    (q_p, trange, R0, r0_p, v0_p, a0_p, r_min)

folder = get_folder(t2, r0_n, v0_n, a0_n, R0)
plot_gauss_theorem_theta_E(q_n, trange, R0, r0_n, v0_n, a0_n, r_min)
plot_gauss_theorem_t_E    (q_n, trange, R0, r0_n, v0_n, a0_n, r_min)


r0_pos_min = r0
v0_pos_min = 0
a0_pos_min = a0
r0_pos_max = r0
v0_pos_max = 0
a0_pos_max = a0

r0_neg_min = r0
v0_neg_min = 0
a0_neg_min = 2 * a0
r0_neg_max = r0
v0_neg_max = 0
a0_neg_max = 2 * a0

step_R0 = R0/2
min_R0 = R0
max_R0 = 2 * R0
'''
# these are results of copper_explosion_lw.sage
r0_pos_min =  0.000000000000000
r0_neg_min =  0.000000000000000
r0_pos_max =  0.000555631898340631
r0_neg_max =  0.000555631898340631

v0_pos_min =  0
v0_neg_min =  0

#v0_pos_max =  0
#v0_neg_max =  0
v0_pos_max =  678.777253763262
v0_neg_max =  231844.794418677

a0_pos_min =  0.000000000000000
a0_neg_min =  0.000000000000000
a0_pos_max =  5.55621349218564e6
a0_neg_max =  1.89779366898370e9

q =  2.12744210611535

# step_R0 = 100*r0_pos_max
step_R0 = 0.05

min_R0 = step_R0
max_R0 = 20.0 * step_R0


attach("spherical_explosion_time_evaluation.sage")
use_dbl_integration = False
spherical_explosion_time_evaluation(q, t1, t2, dt, r0_pos_min, v0_pos_min, a0_pos_min, r0_pos_max, v0_pos_max, a0_pos_max, r0_neg_min, v0_neg_min, a0_neg_min, r0_neg_max, v0_neg_max, a0_neg_max, step_R0, min_R0, max_R0, r_min)
