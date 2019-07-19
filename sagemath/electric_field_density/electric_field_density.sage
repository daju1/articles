import sys
reload(sys)
sys.setdefaultencoding('utf8')
# sage: %display ascii_art

theta_0 = var("theta_0")
v = var("v")
c = var("c")
z = var("z")
r0 = var("r0")
dt_zap = var("dt_zap")

x_1 = solve(x^2 + (x*cot(theta_0) + v*dt_zap)^2 - (r0 + c*dt_zap)^2, x)
z_1 = solve((z*tan(theta_0))^2  + (z + v*dt_zap)^2 - (r0 + c*dt_zap)^2, z)

print "\nx_1 =", x_1
print "\nz_1 =", z_1

# x_1 = [
# x == -(dt_zap*v*cot(theta_0) + sqrt(c^2*dt_zap^2 - dt_zap^2*v^2 + 2*c*dt_zap*r0 + (c^2*dt_zap^2 + 2*c*dt_zap*r0 + r0^2)*cot(theta_0)^2 + r0^2))/(cot(theta_0)^2 + 1),
# x == -(dt_zap*v*cot(theta_0) - sqrt(c^2*dt_zap^2 - dt_zap^2*v^2 + 2*c*dt_zap*r0 + (c^2*dt_zap^2 + 2*c*dt_zap*r0 + r0^2)*cot(theta_0)^2 + r0^2))/(cot(theta_0)^2 + 1)
# ]

def X_1(theta_0, r0, dt_zap, v, c):
    if (theta_0 > pi/2):
        return (dt_zap*v*cot(theta_0) + sqrt(c^2*dt_zap^2 - dt_zap^2*v^2 + 2*c*dt_zap*r0 + (c^2*dt_zap^2 + 2*c*dt_zap*r0 + r0^2)*cot(theta_0)^2 + r0^2))/(cot(theta_0)^2 + 1)
    return    -(dt_zap*v*cot(theta_0) - sqrt(c^2*dt_zap^2 - dt_zap^2*v^2 + 2*c*dt_zap*r0 + (c^2*dt_zap^2 + 2*c*dt_zap*r0 + r0^2)*cot(theta_0)^2 + r0^2))/(cot(theta_0)^2 + 1)


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
         return  (dt_zap*v - sqrt(-dt_zap^2*v^2*tan(theta_0)^2 + c^2*dt_zap^2 + 2*c*dt_zap*r0 + (c^2*dt_zap^2 + 2*c*dt_zap*r0 + r0^2)*tan(theta_0)^2 + r0^2))/(tan(theta_0)^2 + 1)
    #    return -(dt_zap*v + sqrt(-dt_zap^2*v^2*tan(theta_0)^2 + c^2*dt_zap^2 + 2*c*dt_zap*r0 + (c^2*dt_zap^2 + 2*c*dt_zap*r0 + r0^2)*tan(theta_0)^2 + r0^2))/(tan(theta_0)^2 + 1)
    return      -(dt_zap*v - sqrt(-dt_zap^2*v^2*tan(theta_0)^2 + c^2*dt_zap^2 + 2*c*dt_zap*r0 + (c^2*dt_zap^2 + 2*c*dt_zap*r0 + r0^2)*tan(theta_0)^2 + r0^2))/(tan(theta_0)^2 + 1)


#  /               ________________________________________________________________________________________________________________\  #
#  |              /  2       2                         2  2    2              2   / 2       2                     2\    2          |  #
# -\dt_zap*v +- \/  c *dt_zap  + 2*c*dt_zap*r0 - dt_zap *v *tan (theta_0) + r0  + \c *dt_zap  + 2*c*dt_zap*r0 + r0 /*tan (theta_0) /  #
# ----------------------------------------------------------------------------------------------------------------------------------- #
#                                                             2                                                                       #
#                                                          tan (theta_0) + 1                                                          #

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

def theta_2(theta_0, v):
    return theta_1(theta_0, 0.1, 0.01, v, 1)

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


#plot( theta_2(theta_0, 0.5), (theta_0, 0.001, pi-0.001) ).show()

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

