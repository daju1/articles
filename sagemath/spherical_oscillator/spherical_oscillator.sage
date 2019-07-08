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

q = 1
t = 5
R0 = 2
r0 = 1
a0 = 0.1

attach("tzap.spyx")
attach("float_formatting.sage")

philw = phi_lw(q, t, R0, r0, 0)
print "philw =", philw

philw = phi_lw(q, t, R0, r0, a0)
print "philw =", philw

# Data for plotting of t_zap
npoints = 180
thetas = np.arange(0*pi/npoints, npoints*pi/npoints + 1*pi/npoints, 1*pi/npoints)
t_zap_data = [ (thetha_i, tzap(t, R0, r0, a0, thetha_i)) for thetha_i in thetas]

p = list_plot (t_zap_data)
p.save("results/spherical_oscillator_t_zap" + suffix(t, r0, a0) + ".png")

t = 0.2

# Data for plotting of phi_lw of unmoved spherical capacitor
min_R0 = -20.0
max_R0 = 20.0
step_R0 = 0.1
Rneg=2
Rpos=1


phi_lw_data = [ (R0_i, phi_lw(-q, t, R0_i, Rneg, 0) + phi_lw(q, t, R0_i, Rpos, 0)) for R0_i  in np.arange(min_R0, max_R0, step_R0) ]

p = list_plot (phi_lw_data)
p.save("results/spherical_capascitor_phi" + "_Rneg=" + str(Rneg) + "_Rpos= " + str(Rpos) + ".png")

# Data for plotting of phi_lw of spherical oscillator with expanding negative sphere
a0pos = 0.05
a0neg = -0.15
phi_lw_data = [ (R0_i, phi_lw(-q, t, R0_i, Rneg, a0neg) + phi_lw(q, t, R0_i, Rpos, a0pos)) for R0_i in np.arange(min_R0, max_R0, step_R0) ]

p = list_plot (phi_lw_data)
pname = "results/spherical_oscillator_phi" + "_t=" + float_formatting(t) + "_Rneg=" + float_formatting(Rneg) + "_Rpos= " + float_formatting(Rpos) + "_a0neg=" + float_formatting(a0neg) + "_a0pos=" + float_formatting(a0pos) + ".png"
print pname
p.save(pname)

t1 = 0.0
dt = 0.5
t2 = 20.0
phi_lw_data = [ (R0_i, phi_lw(-q, t_i, R0_i, Rneg, a0neg) + phi_lw(q, t_i, R0_i, Rpos, a0pos)) for R0_i in np.arange(min_R0, max_R0, step_R0) for t_i in np.arange(t1, t2, dt)]

p = list_plot (phi_lw_data)
pname = "results/spherical_oscillator_phi" + "_t=" + float_formatting(t1) + ".." + float_formatting(t2) + "_Rneg=" + float_formatting(Rneg) + "_Rpos= " + float_formatting(Rpos) + "_a0neg=" + float_formatting(a0neg) + "_a0pos=" + float_formatting(a0pos) + ".png"
print pname
p.save(pname)

