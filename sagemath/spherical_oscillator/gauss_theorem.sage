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


def plot_gauss_theorem_theta_E(q, R0, r0, v0, a0, r_min):
    npoints = 180
    thetas = np.arange(0*pi/npoints, npoints*pi/npoints + 1*pi/npoints, 1*pi/npoints)
    data_tzap = []
    data_rzap = []
    data_E = []
    data_E_sin = []
    for t in (0, 1, 2):
        for theta in thetas:
            (t2, r1) = tzap(q, t, R0, r0, v0, a0, theta, r_min)
            (phi, A, E1, E2, E, error) = phi_and_E_integrand(q, t, R0, r0, v0, a0, theta, r_min)
            data_E += [ (theta, E/m) ]
            data_E_sin += [ (theta, E/m*sin(theta)) ]
            data_tzap += [ (theta, t2) ]
            data_rzap += [ (theta, r1) ]

    p = list_plot (data_E)
    p.save("results/gauss_theorem_theta_E" + suffix(t, r0, a0) + ".png")

    p = list_plot (data_E_sin)
    p.save("results/gauss_theorem_theta_E_sin" + suffix(t, r0, a0) + ".png")

    p = list_plot (data_tzap)
    p.save("results/gauss_theorem_theta_tzap" + suffix(t, r0, a0) + ".png")

    p = list_plot (data_rzap)
    p.save("results/gauss_theorem_theta_rzap" + suffix(t, r0, a0) + ".png")

q = 1
R0 = 1 * c
r0 = 0.25 * c
v0 = 0
a0 = 0.1 * c
r_min = 0.01
# plot_gauss_theorem_theta_E(q, R0, r0, v0, a0, r_min)

attach("../field_of_deyna_cylinder/num_int.sage")

def phi_and_E_integrand_E(q, t, R0, r0, v0, a0, theta, r_min):
    (phi, A, E1, E2, E, error) = phi_and_E_integrand(q, t, R0, r0, v0, a0, theta, r_min)
    return E

calc_E = lambda q, t, R0, r0, v0, a0, theta, r_min : phi_and_E_integrand_E(q, t, R0, r0, v0, a0, theta, r_min) / m
print (calc_E(q, 0, R0, r0, v0, a0, 0, r_min))


k_E = lambda q, t, R0, r0, v0, a0, r_min : num_int(lambda theta : calc_E ( q, t, R0, r0, v0, a0, theta, r_min) * sin(theta) * R0^2, 0, pi )

data_k = []
for t in np.arange(0, 50, 1):
    k = k_E(q, t, R0, r0, v0, a0, r_min)
    data_k += [(t, k)]

p = list_plot (data_k)
p.save("results/gauss_theorem_k" + suffix(t, r0, a0) + ".png")
