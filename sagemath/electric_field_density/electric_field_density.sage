import sys
reload(sys)
sys.setdefaultencoding('utf8')
# sage: %display ascii_art

theta_0 = var("theta_0")
v = var("v")
c = var("c")
z = var("z")
r0 = var("r0")
dR = var("dR")
dt_zap = var("dt_zap")

x_1 = solve(x^2 + (x*cot(theta_0) + v*dt_zap)^2 - (r0 + c*dt_zap)^2, x)
z_1 = solve((z*tan(theta_0))^2  + (z + v*dt_zap)^2 - (r0 + c*dt_zap)^2, z)

print "\nx_1 =", x_1
print "\nz_1 =", z_1


# x_1 = [
# x == -(dt_zap*v*cot(theta_0) + sqrt(c^2*dt_zap^2 - dt_zap^2*v^2 + 2*c*dt_zap*r0 + (c^2*dt_zap^2 + 2*c*dt_zap*r0 + r0^2)*cot(theta_0)^2 + r0^2))/(cot(theta_0)^2 + 1),
# x == -(dt_zap*v*cot(theta_0) - sqrt(c^2*dt_zap^2 - dt_zap^2*v^2 + 2*c*dt_zap*r0 + (c^2*dt_zap^2 + 2*c*dt_zap*r0 + r0^2)*cot(theta_0)^2 + r0^2))/(cot(theta_0)^2 + 1)
# ]

#  /                            __________________________________________________________________________________________________\  #
#  |                           /  2       2                         2  2     2   / 2       2                     2\    2          |  #
# -\dt_zap*v*cot(theta_0) +- \/  c *dt_zap  + 2*c*dt_zap*r0 - dt_zap *v  + r0  + \c *dt_zap  + 2*c*dt_zap*r0 + r0 /*cot (theta_0) /  #
# ---------------------------------------------------------------------------------------------------------------------------------- #
#                                                             2                                                                      #
#                                                          cot (theta_0) + 1                                                         #

# z_1 = [
# z == -(dt_zap*v + sqrt(-dt_zap^2*v^2*tan(theta_0)^2 + c^2*dt_zap^2 + 2*c*dt_zap*r0 + (c^2*dt_zap^2 + 2*c*dt_zap*r0 + r0^2)*tan(theta_0)^2 + r0^2))/(tan(theta_0)^2 + 1),
# z == -(dt_zap*v - sqrt(-dt_zap^2*v^2*tan(theta_0)^2 + c^2*dt_zap^2 + 2*c*dt_zap*r0 + (c^2*dt_zap^2 + 2*c*dt_zap*r0 + r0^2)*tan(theta_0)^2 + r0^2))/(tan(theta_0)^2 + 1)
# ]

def Z_1(theta_0, r0, dt_zap, v, c):
    if (theta_0 > pi/2):
        return  -(dt_zap*v + sqrt(-dt_zap^2*v^2*tan(theta_0)^2 + c^2*dt_zap^2 + 2*c*dt_zap*r0 + (c^2*dt_zap^2 + 2*c*dt_zap*r0 + r0^2)*tan(theta_0)^2 + r0^2))/(tan(theta_0)^2 + 1)
    return      -(dt_zap*v - sqrt(-dt_zap^2*v^2*tan(theta_0)^2 + c^2*dt_zap^2 + 2*c*dt_zap*r0 + (c^2*dt_zap^2 + 2*c*dt_zap*r0 + r0^2)*tan(theta_0)^2 + r0^2))/(tan(theta_0)^2 + 1)


#  /               ________________________________________________________________________________________________________________\  #
#  |              /  2       2                         2  2    2              2   / 2       2                     2\    2          |  #
# -\dt_zap*v +- \/  c *dt_zap  + 2*c*dt_zap*r0 - dt_zap *v *tan (theta_0) + r0  + \c *dt_zap  + 2*c*dt_zap*r0 + r0 /*tan (theta_0) /  #
# ----------------------------------------------------------------------------------------------------------------------------------- #
#                                                             2                                                                       #
#                                                          tan (theta_0) + 1                                                          #

dR  = solve((r0 + dR)^2 + 2*(r0 + dR) * v * dt_zap * cos(theta_0) + (v * dt_zap)^2 - (r0 + c*dt_zap)^2, dR)
print "\ndR =", dR

# dR = [
# dR == -dt_zap*v*cos(theta_0) - r0 - sqrt(c^2*dt_zap^2 + 2*c*dt_zap*r0 + (dt_zap^2*cos(theta_0)^2 - dt_zap^2)*v^2 + r0^2),
# dR == -dt_zap*v*cos(theta_0) - r0 + sqrt(c^2*dt_zap^2 + 2*c*dt_zap*r0 + (dt_zap^2*cos(theta_0)^2 - dt_zap^2)*v^2 + r0^2)
# ]

dR1 = -dt_zap*v*cos(theta_0) - r0 - sqrt(c^2*dt_zap^2 + 2*c*dt_zap*r0 + (dt_zap^2*cos(theta_0)^2 - dt_zap^2)*v^2 + r0^2)
dR2 = -dt_zap*v*cos(theta_0) - r0 + sqrt(c^2*dt_zap^2 + 2*c*dt_zap*r0 + (dt_zap^2*cos(theta_0)^2 - dt_zap^2)*v^2 + r0^2)

dR1_dt = derivative (dR1, dt_zap)
dR2_dt = derivative (dR2, dt_zap)

print "\ndR1_dt =", dR1_dt
print "\ndR2_dt =", dR2_dt

# dR1_dt = -v*cos(theta_0) - (c^2*dt_zap + (dt_zap*cos(theta_0)^2 - dt_zap)*v^2 + c*r0)/sqrt(c^2*dt_zap^2 + 2*c*dt_zap*r0 + (dt_zap^2*cos(theta_0)^2 - dt_zap^2)*v^2 + r0^2)
# dR2_dt = -v*cos(theta_0) + (c^2*dt_zap + (dt_zap*cos(theta_0)^2 - dt_zap)*v^2 + c*r0)/sqrt(c^2*dt_zap^2 + 2*c*dt_zap*r0 + (dt_zap^2*cos(theta_0)^2 - dt_zap^2)*v^2 + r0^2)

print simplify(limit (dR1_dt, dt_zap = 0))
print simplify(limit (dR2_dt, dt_zap = 0))

# -(sqrt(r0^2)*v*cos(theta_0) + c*r0)/sqrt(r0^2)
# -(sqrt(r0^2)*v*cos(theta_0) - c*r0)/sqrt(r0^2)


def Z_1_0(theta_0, r0, dt_zap, v, c):
    return  -(dt_zap*v + sqrt(-dt_zap^2*v^2*tan(theta_0)^2 + c^2*dt_zap^2 + 2*c*dt_zap*r0 + (c^2*dt_zap^2 + 2*c*dt_zap*r0 + r0^2)*tan(theta_0)^2 + r0^2))/(tan(theta_0)^2 + 1)

def Z_1_1(theta_0, r0, dt_zap, v, c):
    return  -(dt_zap*v - sqrt(-dt_zap^2*v^2*tan(theta_0)^2 + c^2*dt_zap^2 + 2*c*dt_zap*r0 + (c^2*dt_zap^2 + 2*c*dt_zap*r0 + r0^2)*tan(theta_0)^2 + r0^2))/(tan(theta_0)^2 + 1)


def X_1(theta_0, r0, dt_zap, v, c):
    return tan(theta_0) * Z_1(theta_0, r0, dt_zap, v, c)

p = plot( Z_1_0(theta_0, 0.1, 0.1, 0.5, 1), (theta_0, 0.001, pi-0.001) )
pname = "results/plot_of_z_1_0_depending_on_theta_0" + "_r0=0.1_dt_zap=0.5_v=0.5" + ".png"
p.save(pname)

p = plot( Z_1_1(theta_0, 0.1, 0.1, 0.5, 1), (theta_0, 0.001, pi-0.001) )
pname = "results/plot_of_z_1_1_depending_on_theta_0" + "_r0=0.1_dt_zap=0.5_v=0.5" + ".png"
p.save(pname)

p = plot( Z_1(theta_0, 0.1, 0.1, 0.5, 1), (theta_0, 0.001, pi-0.001) )
pname = "results/plot_of_z_1_depending_on_theta_0" + "_r0=0.1_dt_zap=0.5_v=0.5" + ".png"
p.save(pname)

# DOTO: these need to be reworked and tested
x1 = piecewise([((0,pi/2), x_1[0]), ([pi/2,pi], x_1[1])], var=theta_0);  x1
z1 = piecewise([((0,pi/2), z_1[0]), ([pi/2,pi], z_1[1])], var=theta_0);  z1

print "\nx1 =", x1
print "\nz1 =", z1

#x1 = -(dt_zap*v*cot(theta_0) + sqrt(c^2*dt_zap^2 - dt_zap^2*v^2 + 2*c*dt_zap*r0 + (c^2*dt_zap^2 + 2*c*dt_zap*r0 + r0^2)*cot(theta_0)^2 + r0^2))/(cot(theta_0)^2 + 1)
#z1 = -(dt_zap*v + sqrt(-dt_zap^2*v^2*tan(theta_0)^2 + c^2*dt_zap^2 + 2*c*dt_zap*r0 + (c^2*dt_zap^2 + 2*c*dt_zap*r0 + r0^2)*tan(theta_0)^2 + r0^2))/(tan(theta_0)^2 + 1)

z1_ = x1 * cot(theta_0)
x1_ = z1 * tan(theta_0)

tan_theta_1 = x1/(z1+v*dt_zap)
# cot_theta_1 = (z1 + v*dt_zap)/x1
cot_theta_1 = cot(theta_0) + (v*dt_zap)/x1
Cot_theta_1 = cot(theta_0) + (v*dt_zap) / X_1(theta_0, r0, dt_zap, v, c)

def theta_1(theta_0, r0, dt_zap, v, c):
    #return X_1(theta_0, r0, dt_zap, v, c) / Z_1(theta_0, r0, dt_zap, v, c)
    _x1 = X_1(theta_0, r0, dt_zap, v, c)
    _z1 = Z_1(theta_0, r0, dt_zap, v, c)
    #print "\ntheta_0, _x1, _z1", theta_0, _x1.n(), _z1.n()
    #return arctan2(_x1, _z1)
    return arctan2(_x1, _z1 + v*dt_zap)

#def theta_1(theta_0, r0, dt_zap, v, c):
#    return arccot(cot(theta_0) + (v*dt_zap) / X_1(theta_0, r0, dt_zap, v, c))


def theta_5(theta_0, v):
    return theta_1(theta_0, 0.1, 0.5, v, 1)#/theta_0

print [theta_5(pi/4, 0.5), theta_5(3*pi/4, 0.5)]

#plot( [theta_5(pi/4, v), theta_5(3*pi/4, v)], (v,0,0.999) ).show()

print( [
    #theta_5(0*pi/16, 0.5), #  indeterminate expression: 0 * infinity encountered.
    theta_5(1*pi/16, 0.5),
    theta_5(2*pi/16, 0.5),
    theta_5(3*pi/16, 0.5),
    theta_5(4*pi/16, 0.5),
    theta_5(5*pi/16, 0.5),
    theta_5(6*pi/16, 0.5),
    theta_5(7*pi/16, 0.5),
    #theta_5(8*pi/16, 0.5), #  indeterminate expression: 0 * infinity encountered. #_z1 = Z_1(theta_0, r0, dt_zap, v, c) #  return      -(dt_zap*v - sqrt(-dt_zap**_sage_const_2 *v**_sage_const_2 *tan(theta_0)**_sage_const_2  + c**_sage_const_2 *dt_zap**_sage_const_2  + _sage_const_2 *c*dt)
    theta_5(9*pi/16, 0.5),
    theta_5(10*pi/16, 0.5),
    theta_5(11*pi/16, 0.5),
    theta_5(12*pi/16, 0.5),
    theta_5(14*pi/16, 0.5),
    theta_5(15*pi/16, 0.5),
    ])

p = plot( [
    #theta_5(0*pi/16, v), #  indeterminate expression: 0 * infinity encountered.
    theta_5(1*pi/16, v),
    theta_5(2*pi/16, v),
    theta_5(3*pi/16, v),
    theta_5(4*pi/16, v),
    theta_5(5*pi/16, v),
    theta_5(6*pi/16, v),
    theta_5(7*pi/16, v),
    #theta_5(8*pi/16, v), #  indeterminate expression: 0 * infinity encountered. #_z1 = Z_1(theta_0, r0, dt_zap, v, c) #  return      -(dt_zap*v - sqrt(-dt_zap**_sage_const_2 *v**_sage_const_2 *tan(theta_0)**_sage_const_2  + c**_sage_const_2 *dt_zap**_sage_const_2  + _sage_const_2 *c*dt)
    theta_5(9*pi/16, v),
    theta_5(10*pi/16, v),
    theta_5(11*pi/16, v),
    theta_5(12*pi/16, v),
    theta_5(14*pi/16, v),
    theta_5(15*pi/16, v),
    ], (v,0,0.999))

pname = "results/theta_1_depending_on_v_" + "_r0=0.1_dt_zap=0.5" + ".png"
p.save(pname)

p = contour_plot( theta_1(theta_0, 0.1, 0.5, v, 1), (theta_0, 0.001, pi-0.001), (v,0,0.999) )
pname = "results/contour_plot_of_theta_1_depending_on_v_and_theta_0" + "_r0=0.1_dt_zap=0.5" + ".png"
p.save(pname)

#n = simplify(x1/(v*dt_zap))
#print n

# -(dt_zap*v*cot(theta_0) + sqrt(c^2*dt_zap^2 - dt_zap^2*v^2 + 2*c*dt_zap*r0 + (c^2*dt_zap^2 + 2*c*dt_zap*r0 + r0^2)*cot(theta_0)^2 + r0^2))/((cot(theta_0)^2 + 1)*dt_zap*v)

#  /                           __________________________________________________________________________________________________\
#  |                          /  2       2                         2  2     2   / 2       2                     2\    2          |
# -\dt_zap*v*cot(theta_0) + \/  c *dt_zap  + 2*c*dt_zap*r0 - dt_zap *v  + r0  + \c *dt_zap  + 2*c*dt_zap*r0 + r0 /*cot (theta_0) /
# ---------------------------------------------------------------------------------------------------------------------------------
#                                                             /   2             \
#                                                    dt_zap*v*\cot (theta_0) + 1/

import numpy as np
from sage.plot.circle import Circle
g = Graphics()
min_theta_0 = 5*pi/180
max_theta_0=176*pi/180
step_theta_0=10*pi/180

R0 = 100
dtzap = 0.1

z_plot_data = []
x_plot_data = []
for theta0 in np.arange (min_theta_0, max_theta_0, step_theta_0):
    z1=Z_1(theta0, R0, dtzap, 0.5, 1)
    x1=X_1(theta0, R0, dtzap, 0.5, 1)
    z_plot_data += [(theta0, z1)]
    x_plot_data += [(theta0, x1)]

p = list_plot( z_plot_data )
pname = "results/list_plot_of_z_1_depending_on_theta_0" + "_r0=100_dt_zap=0.5_v=0.5" + ".png"
p.save(pname)

p = list_plot( x_plot_data )
pname = "results/list_plot_of_x_1_depending_on_theta_0" + "_r0=100_dt_zap=0.5_v=0.5" + ".png"
p.save(pname)

for theta0 in np.arange (min_theta_0, max_theta_0, step_theta_0):
    zq = 0
    za = 0
    xa = 0
    vv = 0.8
    cc = 1
    theta = theta0.n()
    r = R0
    tt = 0
    plot_data = []
    for i in range (0, 1001):
        tt += dtzap
        _x1 = X_1(theta, r, dtzap, vv, cc)
        _z1 = Z_1(theta, r, dtzap, vv, cc)
        theta1 = arctan2(_x1, _z1 + vv*dtzap)
        theta = theta1
        za = _z1 - zq
        xa = _x1
        zq += dtzap*vv
        r = R0 + cc * tt
        # print "\ntt theta r rr _x1 _z1 xa za zq", tt, theta, r, rr, _x1, _z1, xa, za, zq
        rr = sqrt(_x1^2 + _z1^2)

        if (0 == i % 100 and theta0 == min_theta_0):
	    g += circle((0, -zq), rr, rgbcolor=hue(r/500))

        plot_data += [(xa, za)]

    print "theta0 = ", theta0

    g += list_plot(plot_data)

pname = "results/normals" + ".png"
print pname
g.save(pname)
# g.show()
