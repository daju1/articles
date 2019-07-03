import sys
reload(sys)
sys.setdefaultencoding('utf8')

import numpy as np

c = var("c")
t = var("t")
R0 = var("R0")
r0 = var("r0")
theta=var("theta")
t_zap = var("t_zap")

t = 5
R0 = 2
r0 = 1
a0 = 2
attach("tzap.spyx")

# Data for plotting of t_zap
npoints = 180
thetas = np.arange(0*pi/npoints, npoints*pi/npoints + 1*pi/npoints, 1*pi/npoints)
t_zap_data = [ (thetha_i, tzap(t, R0, r0, a0, thetha_i)) for thetha_i in thetas]

#print "t_zap_data = ", t_zap_data

p = list_plot (t_zap_data)
p.save("results/spherical_oscillator_t_zap" + suffix(t, r0, a0) + ".png")


t = 5
q = 1

# Data for plotting of phi_lw of unmoved spherical capacitor
min_R0 = -5.0
max_R0 = 5.0
step_R0 = 0.01

phi_lw_data = [ (R0_i, phi_lw(-q, t, R0_i, 2, 0) + phi_lw(q, t, R0_i, 1, 0)) for R0_i  in np.arange(min_R0, max_R0, step_R0) ]
#print "phi_lw_data = ", phi_lw_data

p = list_plot (phi_lw_data)
try:
    p.show()
except ValueError as ex:
    print ex

p.save("results/spherical_capascitor_phi" + "_Rneg=2" + "_Rpos=1 " + ".png")


philw = phi_lw(q, t, R0, r0, 0)
print "philw =", philw

philw = phi_lw(q, t, R0, r0, a0)
print "philw =", philw
