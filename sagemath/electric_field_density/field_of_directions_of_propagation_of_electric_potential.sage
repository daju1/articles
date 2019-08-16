import sys
reload(sys)
sys.setdefaultencoding('utf8')

v = var("v")
c = var("c")
x = var("x")
z = var("z")
r0 = var("r0")
t = var('t')
t_z = var('t_z')
ttz = var('ttz')
tt = var('tt')
C_0 = var("C_0")
theta_0 = var("theta_0")
theta_z = var("theta_z")
theta_1 = var("theta_1")

assume(v > 0, c > 0, v < c)
#Is ttz positive, negative or zero?
assume(ttz>0)
# Is c*ttz+r0 positive, negative or zero?
#assume(c*ttz+r0>0)
# Is r0 positive, negative or zero?
assume(r0>0)
#Is theta_1-theta_0 positive, negative or zero?
#assume(theta_1<theta_0)
#Is theta_0 positive, negative or zero?
assume(theta_0>0)
assume(theta_0<pi)
#Is theta_1 positive, negative or zero?
assume(theta_1>0)
assume(theta_1<pi)

#theta_z = function('theta_z')(tt)
#de = diff(theta_z, tt) - (v / (r0 + c*tt)) / (sin(theta_z) + ((cos(theta_z))^2) / sin(theta_z))
#desolve(de, theta_z)

eq = integral( (v / (r0 + c*tt)), tt, 0, ttz) + integral( (sin(theta_z) + ((cos(theta_z))^2) / sin(theta_z)), theta_z, theta_0, theta_1) - C_0
print "\neq  =", eq 

#sage: 
#eq  = v*(log(c*ttz + r0)/c - log(r0)/c) - C_0 + 1/2*log(abs(cos(theta_0) + 1)) - 1/2*log(abs(cos(theta_0) - 1)) - 1/2*log(abs(cos(theta_1) + 1)) + 1/2*log(abs(cos(theta_1) - 1))

#sage: sage: %display ascii_art
#sage: eq
#
#                /  log(r0)   log(c*ttz + r0)\   log(|cos(theta_0) - 1|)   log(|cos(theta_0) + 1|)   log(|cos(theta_1) - 1|)   log(|cos(theta_1) + 1|)
#-C_0 + v*|- ------- + -------------------| - -------------------------- + ---------------------------- + ---------------------------- - -----------------------------
#                \     c             c                /                  2                                2                                    2                                2           
#

eq1 = (((c*ttz + r0) / r0) ^ (2*v/c)) * (1+cos(theta_0)) / (1-cos(theta_0)) * (1-cos(theta_1)) / (1+cos(theta_1))-1
print "\neq1  =", eq1

th1 = solve(eq1, theta_1)
print "\nth1  =", th1

#th1  = [
#theta_1 == arccos((c*ttz/r0 + 1)^(2*v/c)*cos(theta_0)/((c*ttz/r0 + 1)^(2*v/c)*cos(theta_0) + (c*ttz/r0 + 1)^(2*v/c) - cos(theta_0) + 1) + (c*ttz/r0 + 1)^(2*v/c)/((c*ttz/r0 + 1)^(2*v/c)*cos(theta_0) + (c*ttz/r0 + 1)^(2*v/c) - cos(theta_0) + 1) + cos(theta_0)/((c*ttz/r0 + 1)^(2*v/c)*cos(theta_0) + (c*ttz/r0 + 1)^(2*v/c) - cos(theta_0) + 1) - 1/((c*ttz/r0 + 1)^(2*v/c)*cos(theta_0) + (c*ttz/r0 + 1)^(2*v/c) - cos(theta_0) + 1))
#]
