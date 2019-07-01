import sys
reload(sys)
sys.setdefaultencoding('utf8')


c = var("c")
t = var("t")
R0 = var("R0")
r0 = var("r0")
theta=var("theta")
t_zap = var("t_zap")

t = 5
R0 = 2
r0 = 1
a0 = 0.1
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

#ftzap(t, R, r, a, th) = tzap(t, R, r, a, th)

#plot(ftzap(t, R0, r0, a0, theta), (theta, 0, pi))
