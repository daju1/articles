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
tt1 = var('tt1')
tt2 = var('tt2')

C_0 = var("C_0")
theta_0 = var("theta_0")
theta_z = var("theta_z")
theta_1 = var("theta_1")
theta_2 = var("theta_2")

assume(v > 0, c > 0, v < c)
#Is ttz positive, negative or zero?
assume(ttz>0)
# Is c*ttz+r0 positive, negative or zero?
#assume(c*ttz+r0>0)
# Is r0 positive, negative or zero?
assume(r0>0)
#Is theta_1-theta_0 positive, negative or zero?
assume(theta_1<theta_0)
#Is theta_0 positive, negative or zero?
assume(theta_0>0)
assume(theta_0<pi)
#Is theta_1 positive, negative or zero?
assume(theta_1>0)
#assume(theta_1<pi)

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


tan_theta_2_half = tan(theta_1 / 2) * ( (r0+c*tt1) / (r0+c*tt2) ) ^ (v / c)

sin_theta_2  = 2 * tan_theta_2_half / (1 + tan_theta_2_half ^2)

cos_theta_2  = (1 - tan_theta_2_half ^2) / (1 + tan_theta_2_half ^2)

x2 = (r0 + c * tt2) * sin_theta_2
z2 = (r0 + c * tt2) * cos_theta_2 - v * tt2

print "\ntan_theta_2_half = ", tan_theta_2_half
print "\nsin_theta_2 = ", sin_theta_2
print "\ncos_theta_2 = ", cos_theta_2
print "\nx2 = ", x2
print "\nz2 = ", z2

# tan_theta_2_half =  ((c*tt1 + r0)/(c*tt2 + r0))^(v/c)*tan(1/2*theta_1)
# sin_theta_2 =  2*((c*tt1 + r0)/(c*tt2 + r0))^(v/c)*tan(1/2*theta_1)/(((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)
# cos_theta_2 =  -(((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 - 1)/(((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)
# x2 =  2*(c*tt2 + r0)*((c*tt1 + r0)/(c*tt2 + r0))^(v/c)*tan(1/2*theta_1)/(((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)
# z2 =  -tt2*v - (((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 - 1)*(c*tt2 + r0)/(((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)

dx2_dtt2 = derivative(x2, tt2)
dz2_dtt2 = derivative(z2, tt2)

dx2_dtheta_1 = derivative(x2, theta_1)
dz2_dtheta_1 = derivative(z2, theta_1)

print "\ndx2_dtt2 = ", dx2_dtt2
print "\ndz2_dtt2 = ", dz2_dtt2

print "\ndx2_dtheta_1 = ", dx2_dtheta_1
print "\ndz2_dtheta_1 = ", dz2_dtheta_1

# dx2_dtt2 =  4*(c*tt1 + r0)*v*((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c - 1)*((c*tt1 + r0)/(c*tt2 + r0))^(v/c)*tan(1/2*theta_1)^3/((((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)^2*(c*tt2 + r0)) - 2*(c*tt1 + r0)*v*((c*tt1 + r0)/(c*tt2 + r0))^(v/c - 1)*tan(1/2*theta_1)/((((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)*(c*tt2 + r0)) + 2*c*((c*tt1 + r0)/(c*tt2 + r0))^(v/c)*tan(1/2*theta_1)/(((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)
# dz2_dtt2 =  2*(c*tt1 + r0)*v*((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c - 1)*tan(1/2*theta_1)^2/((((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)*(c*tt2 + r0)) - 2*(((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 - 1)*(c*tt1 + r0)*v*((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c - 1)*tan(1/2*theta_1)^2/((((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)^2*(c*tt2 + r0)) - (((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 - 1)*c/(((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1) - v
# dx2_dtheta_1 =  -2*(c*tt2 + r0)*(tan(1/2*theta_1)^2 + 1)*((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*((c*tt1 + r0)/(c*tt2 + r0))^(v/c)*tan(1/2*theta_1)^2/(((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)^2 + (c*tt2 + r0)*(tan(1/2*theta_1)^2 + 1)*((c*tt1 + r0)/(c*tt2 + r0))^(v/c)/(((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)
# dz2_dtheta_1 =  -(c*tt2 + r0)*(tan(1/2*theta_1)^2 + 1)*((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)/(((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1) + (((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 - 1)*(c*tt2 + r0)*(tan(1/2*theta_1)^2 + 1)*((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)/(((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)^2

J = dx2_dtt2 * dz2_dtheta_1 - dx2_dtheta_1 * dz2_dtt2

print "\nJ = ", J

# J =  -2*(2*(c*tt1 + r0)*v*((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c - 1)*((c*tt1 + r0)/(c*tt2 + r0))^(v/c)*tan(1/2*theta_1)^3/((((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)^2*(c*tt2 + r0)) - (c*tt1 + r0)*v*((c*tt1 + r0)/(c*tt2 + r0))^(v/c - 1)*tan(1/2*theta_1)/((((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)*(c*tt2 + r0)) + c*((c*tt1 + r0)/(c*tt2 + r0))^(v/c)*tan(1/2*theta_1)/(((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1))*((c*tt2 + r0)*(tan(1/2*theta_1)^2 + 1)*((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)/(((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1) - (((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 - 1)*(c*tt2 + r0)*(tan(1/2*theta_1)^2 + 1)*((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)/(((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)^2) + (2*(c*tt2 + r0)*(tan(1/2*theta_1)^2 + 1)*((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*((c*tt1 + r0)/(c*tt2 + r0))^(v/c)*tan(1/2*theta_1)^2/(((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)^2 - (c*tt2 + r0)*(tan(1/2*theta_1)^2 + 1)*((c*tt1 + r0)/(c*tt2 + r0))^(v/c)/(((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1))*(2*(c*tt1 + r0)*v*((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c - 1)*tan(1/2*theta_1)^2/((((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)*(c*tt2 + r0)) - 2*(((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 - 1)*(c*tt1 + r0)*v*((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c - 1)*tan(1/2*theta_1)^2/((((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)^2*(c*tt2 + r0)) - (((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 - 1)*c/(((c*tt1 + r0)/(c*tt2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1) - v
