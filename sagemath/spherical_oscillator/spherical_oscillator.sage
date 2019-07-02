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

t_zap = tzap(t, R0, r0, a0, 0*pi/8)
print "tzap (0*pi/8)=", t_zap

t_zap = tzap(t, R0, r0, a0, 1*pi/8)
print "tzap (1*pi/8)=", t_zap

t_zap = tzap(t, R0, r0, a0, 2*pi/8)
print "tzap (2*pi/8)=", t_zap

t_zap = tzap(t, R0, r0, a0, 3*pi/8)
print "tzap (3*pi/8)=", t_zap

t_zap = tzap(t, R0, r0, a0, 4*pi/8)
print "tzap (4*pi/8)=", t_zap

t_zap = tzap(t, R0, r0, a0, 5*pi/8)
print "tzap (5*pi/8)=", t_zap

t_zap = tzap(t, R0, r0, a0, 6*pi/8)
print "tzap (6*pi/8)=", t_zap

t_zap = tzap(t, R0, r0, a0, 7*pi/8)
print "tzap (7*pi/8)=", t_zap

t_zap = tzap(t, R0, r0, a0, 8*pi/8)
print "tzap (8*pi/8)=", t_zap

# Data for plotting
npoints = 180
thetas = np.arange(0*pi/npoints, npoints*pi/npoints + 1*pi/npoints, 1*pi/npoints)
t_zap_data = [ (tzap(t, R0, r0, a0, thetha_i), thetha_i) for thetha_i in thetas]

print t_zap_data

p = list_plot (t_zap_data)
print p
p.show()
