import sys
reload(sys)
sys.setdefaultencoding('utf8')

import numpy as np
attach("../spherical_oscillator/float_formatting.sage")

v = var("v")
c = var("c")
x = var("x")
z = var("z")
r0 = var("r0")
t = var('t')
t_z = var('t_z')
tauz = var('tauz')
tau = var('tau')
tau1 = var('tau1')
tau2 = var('tau2')

C_0 = var("C_0")
theta_0 = var("theta_0")
theta_z = var("theta_z")
theta_1 = var("theta_1")
theta_2 = var("theta_2")

assume(v > 0, c > 0, v < c)
#Is tauz positive, negative or zero?
assume(tauz>0)
# Is c*tauz+r0 positive, negative or zero?
#assume(c*tauz+r0>0)
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

assume((c*tau1+r0)>0)
assume((c*tau2+r0)>0)

#theta_z = function('theta_z')(tau)
#de = diff(theta_z, tau) - (v / (r0 + c*tau)) / (sin(theta_z) + ((cos(theta_z))^2) / sin(theta_z))
#desolve(de, theta_z)

eq = integral( (v / (r0 + c*tau)), tau, 0, tauz) + integral( (sin(theta_z) + ((cos(theta_z))^2) / sin(theta_z)), theta_z, theta_0, theta_1) - C_0
print "\neq  =", eq 

#sage: 
#eq  = v*(log(c*tauz + r0)/c - log(r0)/c) - C_0 + 1/2*log(abs(cos(theta_0) + 1)) - 1/2*log(abs(cos(theta_0) - 1)) - 1/2*log(abs(cos(theta_1) + 1)) + 1/2*log(abs(cos(theta_1) - 1))

#sage: sage: %display ascii_art
#sage: eq
#
#                /  log(r0)   log(c*tauz + r0)\   log(|cos(theta_0) - 1|)   log(|cos(theta_0) + 1|)   log(|cos(theta_1) - 1|)   log(|cos(theta_1) + 1|)
#-C_0 + v*|- ------- + -------------------| - -------------------------- + ---------------------------- + ---------------------------- - -----------------------------
#                \     c             c                /                  2                                2                                    2                                2           
#

eq1 = (((c*tauz + r0) / r0) ^ (2*v/c)) * (1+cos(theta_0)) / (1-cos(theta_0)) * (1-cos(theta_1)) / (1+cos(theta_1))-1
print "\neq1  =", eq1

th1 = solve(eq1, theta_1)
print "\nth1  =", th1

#th1  = [
#theta_1 == arccos((c*tauz/r0 + 1)^(2*v/c)*cos(theta_0)/((c*tauz/r0 + 1)^(2*v/c)*cos(theta_0) + (c*tauz/r0 + 1)^(2*v/c) - cos(theta_0) + 1) + (c*tauz/r0 + 1)^(2*v/c)/((c*tauz/r0 + 1)^(2*v/c)*cos(theta_0) + (c*tauz/r0 + 1)^(2*v/c) - cos(theta_0) + 1) + cos(theta_0)/((c*tauz/r0 + 1)^(2*v/c)*cos(theta_0) + (c*tauz/r0 + 1)^(2*v/c) - cos(theta_0) + 1) - 1/((c*tauz/r0 + 1)^(2*v/c)*cos(theta_0) + (c*tauz/r0 + 1)^(2*v/c) - cos(theta_0) + 1))
#]

# theta_1 - координатный угол криволинейной сиситемы координат, определяемый углом наклона кривой направления распространения электрического потенциала в момент времени tau1
# theta_2 - текущий угол наклона кривой направления распространения электрического потенциала в момент времени tau2

# интеграл дифференциального уравнения кривой направления распространения электрического потенциала
tan_theta_2_half = tan(theta_1 / 2) * ( (r0+c*tau1) / (r0+c*tau2) ) ^ (v / c)

# синус и косинус текущего угла наклона  кривой направления распространения электрического потенциала
sin_theta_2  = 2 * tan_theta_2_half / (1 + tan_theta_2_half ^2)
cos_theta_2  = (1 - tan_theta_2_half ^2) / (1 + tan_theta_2_half ^2)

# декартовы координаты кривой направления распространения электрического потенциала
x2 = (r0 + c * tau2) * sin_theta_2
z2 = (r0 + c * tau2) * cos_theta_2 - v * tau2

print "\ntan_theta_2_half = ", tan_theta_2_half
print "\nsin_theta_2 = ", sin_theta_2
print "\ncos_theta_2 = ", cos_theta_2
print "\nx2 = ", x2
print "\nz2 = ", z2

# tan_theta_2_half =  ((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*tan(1/2*theta_1)
# sin_theta_2 =  2*((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*tan(1/2*theta_1)/(((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)
# cos_theta_2 =  -(((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 - 1)/(((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)
# x2 =  2*(c*tau2 + r0)*((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*tan(1/2*theta_1)/(((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)
# z2 =  -tau2*v - (((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 - 1)*(c*tau2 + r0)/(((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)

# составляем якобиан криволинейной системы координат
# как определитель матрицы производных декартовых координат x2 и z2 по криволинейным координатам  tau2 и theta_1
dx2_dtau2 = derivative(x2, tau2)
dz2_dtau2 = derivative(z2, tau2)

dx2_dtheta_1 = derivative(x2, theta_1)
dz2_dtheta_1 = derivative(z2, theta_1)

print "\ndx2_dtau2 = ", dx2_dtau2
print "\ndz2_dtau2 = ", dz2_dtau2

print "\ndx2_dtheta_1 = ", dx2_dtheta_1
print "\ndz2_dtheta_1 = ", dz2_dtheta_1

# dx2_dtau2 =  4*(c*tau1 + r0)*v*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c - 1)*((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*tan(1/2*theta_1)^3/((((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)^2*(c*tau2 + r0)) - 2*(c*tau1 + r0)*v*((c*tau1 + r0)/(c*tau2 + r0))^(v/c - 1)*tan(1/2*theta_1)/((((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)*(c*tau2 + r0)) + 2*c*((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*tan(1/2*theta_1)/(((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)
# dz2_dtau2 =  2*(c*tau1 + r0)*v*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c - 1)*tan(1/2*theta_1)^2/((((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)*(c*tau2 + r0)) - 2*(((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 - 1)*(c*tau1 + r0)*v*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c - 1)*tan(1/2*theta_1)^2/((((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)^2*(c*tau2 + r0)) - (((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 - 1)*c/(((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1) - v
# dx2_dtheta_1 =  -2*(c*tau2 + r0)*(tan(1/2*theta_1)^2 + 1)*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*tan(1/2*theta_1)^2/(((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)^2 + (c*tau2 + r0)*(tan(1/2*theta_1)^2 + 1)*((c*tau1 + r0)/(c*tau2 + r0))^(v/c)/(((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)
# dz2_dtheta_1 =  -(c*tau2 + r0)*(tan(1/2*theta_1)^2 + 1)*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)/(((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1) + (((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 - 1)*(c*tau2 + r0)*(tan(1/2*theta_1)^2 + 1)*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)/(((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)^2

# якобиан криволинейной системы координат
J = dx2_dtau2 * dz2_dtheta_1 - dx2_dtheta_1 * dz2_dtau2
print "\nJ = ", J

# J =  -2*(2*(c*tau1 + r0)*v*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c - 1)*((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*tan(1/2*theta_1)^3/((((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)^2*(c*tau2 + r0)) - (c*tau1 + r0)*v*((c*tau1 + r0)/(c*tau2 + r0))^(v/c - 1)*tan(1/2*theta_1)/((((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)*(c*tau2 + r0)) + c*((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*tan(1/2*theta_1)/(((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1))*((c*tau2 + r0)*(tan(1/2*theta_1)^2 + 1)*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)/(((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1) - (((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 - 1)*(c*tau2 + r0)*(tan(1/2*theta_1)^2 + 1)*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)/(((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)^2) + (2*(c*tau2 + r0)*(tan(1/2*theta_1)^2 + 1)*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*tan(1/2*theta_1)^2/(((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)^2 - (c*tau2 + r0)*(tan(1/2*theta_1)^2 + 1)*((c*tau1 + r0)/(c*tau2 + r0))^(v/c)/(((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1))*(2*(c*tau1 + r0)*v*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c - 1)*tan(1/2*theta_1)^2/((((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)*(c*tau2 + r0)) - 2*(((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 - 1)*(c*tau1 + r0)*v*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c - 1)*tan(1/2*theta_1)^2/((((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)^2*(c*tau2 + r0)) - (((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 - 1)*c/(((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1) - v

J = J.substitute(tan(1/2*theta_1) == sqrt((1-cos(theta_1))/(1+cos(theta_1))) ).full_simplify()
print "\nJ = ", J

# J =  -2*(((c*tau2 - (c*tau2 + r0)*cos(theta_1) + r0)*((c*tau1 + r0)/(c*tau2 + r0))^(3*v/c) - (c*tau2 + (c*tau2 + r0)*cos(theta_1) + r0)*((c*tau1 + r0)/(c*tau2 + r0))^(v/c))*v + (c^2*tau2 + c*r0 - (c^2*tau2 + c*r0)*cos(theta_1))*((c*tau1 + r0)/(c*tau2 + r0))^(3*v/c) + (c^2*tau2 + c*r0 + (c^2*tau2 + c*r0)*cos(theta_1))*((c*tau1 + r0)/(c*tau2 + r0))^(v/c))/(2*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*sin(theta_1)^2 - (sin(theta_1)^2 + 2*cos(theta_1) - 2)*((c*tau1 + r0)/(c*tau2 + r0))^(4*v/c) - sin(theta_1)^2 + 2*cos(theta_1) + 2)




# предел якобиана при устремлении разности текущего момента и запаздывающего момента в точке наблюдения (tau2 = t - t_zap_2)
# к разности текущего момента и момента начала движения заряда (tau1 = t - t_zap_1)
J0 = J.limit(tau2 = tau1)
print "\nJ0 = ", J0

# J0 =  -(c*r0*tan(1/2*theta_1)^2 + c*r0 + (c^2*tan(1/2*theta_1)^2 + c^2)*tau1 + (r0*tan(1/2*theta_1)^2 + (c*tan(1/2*theta_1)^2 - c)*tau1 - r0)*v)/(tan(1/2*theta_1)^2 + 1)

# sage: %display ascii_art
# sage: J0
#
# /        2/theta_1\              / 2    2/theta_1\    2\     /      2/theta_1\            /     2/theta_1\    \\\ 
#-|c*r0*tan |-------| + c*r0 + tau1*|c *tan |-------| + c | + v*|r0*tan |-------| - r0 + tau1*|c*tan |-------| - c|||
# \         \   2   /              \       \   2   /     /     \       \   2   /            \      \   2   /    /// 
#-------------------------------------------------------------------------------------------------------------------
#                                                    2/theta_1\                                                     
#                                                 tan |-------| + 1                                                 
#                                                     \   2   /              

J0 = J0.substitute(tan(1/2*theta_1) == sqrt((1-cos(theta_1))/(1+cos(theta_1))) ).full_simplify()
print "\nJ0 = ", J0
# J0 =  -c^2*tau1 + (c*tau1 + r0)*v*cos(theta_1) - c*r0

J0 = J0.collect(r0+c*tau1)
print "\nJ0 = ", J0

R_zap_1 = var("R_zap_1")
J0 = J0.substitute((c*tau1 + r0)==R_zap_1)
print "\nJ0 = ", J0

# J0 = - (r0+c*tau1) * (c - v * cos(theta_1)


# метрический тензор
# x_1 = tau
# x_2 = theta_1

g11 = dx2_dtau2    * dx2_dtau2    + dz2_dtau2    * dz2_dtau2
g22 = dx2_dtheta_1 * dx2_dtheta_1 + dz2_dtheta_1 * dz2_dtheta_1

g12 = dx2_dtau2    * dx2_dtheta_1 + dz2_dtau2    * dz2_dtheta_1
g21 = dx2_dtheta_1 * dx2_dtau2    + dz2_dtheta_1 * dz2_dtau2

print "\ng11 = ", g11
print "\ng22 = ", g22
print "\ng12 = ", g12
print "\ng21 = ", g21

g11 = g11.substitute(tan(1/2*theta_1) == sqrt((1-cos(theta_1))/(1+cos(theta_1))) ).full_simplify()
g22 = g22.substitute(tan(1/2*theta_1) == sqrt((1-cos(theta_1))/(1+cos(theta_1))) ).full_simplify()
g12 = g12.substitute(tan(1/2*theta_1) == sqrt((1-cos(theta_1))/(1+cos(theta_1))) ).full_simplify()
g21 = g21.substitute(tan(1/2*theta_1) == sqrt((1-cos(theta_1))/(1+cos(theta_1))) ).full_simplify()

print "\ng11 = ", g11
print "\ng22 = ", g22
print "\ng12 = ", g12
print "\ng21 = ", g21

# g11 =  (2*c^2*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*sin(theta_1)^2 - c^2*sin(theta_1)^2 - (2*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*sin(theta_1)^2 + (sin(theta_1)^2 + 2*cos(theta_1) - 2)*((c*tau1 + r0)/(c*tau2 + r0))^(4*v/c) + sin(theta_1)^2 - 2*cos(theta_1) - 2)*v^2 + 2*c^2*cos(theta_1) + 2*c^2 + 2*(c*sin(theta_1)^2 - (c*sin(theta_1)^2 + 2*c*cos(theta_1) - 2*c)*((c*tau1 + r0)/(c*tau2 + r0))^(4*v/c) - 2*c*cos(theta_1) - 2*c)*v - (c^2*sin(theta_1)^2 + 2*c^2*cos(theta_1) - 2*c^2)*((c*tau1 + r0)/(c*tau2 + r0))^(4*v/c))/(2*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*sin(theta_1)^2 - (sin(theta_1)^2 + 2*cos(theta_1) - 2)*((c*tau1 + r0)/(c*tau2 + r0))^(4*v/c) - sin(theta_1)^2 + 2*cos(theta_1) + 2)

# g22 =  4*(c^2*tau2^2 + 2*c*r0*tau2 + r0^2)*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)/(2*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*sin(theta_1)^2 - (sin(theta_1)^2 + 2*cos(theta_1) - 2)*((c*tau1 + r0)/(c*tau2 + r0))^(4*v/c) - sin(theta_1)^2 + 2*cos(theta_1) + 2)

# g12 =  -8*(((c*tau2 + (c*tau2 + r0)*cos(theta_1) + r0)*(-(cos(theta_1) - 1)/(cos(theta_1) + 1))^(3/2) - (c*tau2 - (c*tau2 + r0)*cos(theta_1) + r0)*sqrt(-(cos(theta_1) - 1)/(cos(theta_1) + 1)))*((c*tau1 + r0)/(c*tau2 + r0))^(6*v/c)*sin(theta_1)^2 + (((3*c*tau2 + (c*tau2 + r0)*cos(theta_1) + 3*r0)*(-(cos(theta_1) - 1)/(cos(theta_1) + 1))^(3/2) + (c*tau2 + (c*tau2 + r0)*cos(theta_1) + r0)*sqrt(-(cos(theta_1) - 1)/(cos(theta_1) + 1)))*sin(theta_1)^2 - 4*(c*tau2 + (c*tau2 + r0)*cos(theta_1) + r0)*(-(cos(theta_1) - 1)/(cos(theta_1) + 1))^(3/2))*((c*tau1 + r0)/(c*tau2 + r0))^(4*v/c))*v/(6*((c*tau1 + r0)/(c*tau2 + r0))^(4*v/c)*sin(theta_1)^4 + sin(theta_1)^4 - 4*(cos(theta_1) + 2)*sin(theta_1)^2 + (sin(theta_1)^4 + 4*(cos(theta_1) - 2)*sin(theta_1)^2 - 8*cos(theta_1) + 8)*((c*tau1 + r0)/(c*tau2 + r0))^(8*v/c) - 4*(sin(theta_1)^4 + 2*(cos(theta_1) - 1)*sin(theta_1)^2)*((c*tau1 + r0)/(c*tau2 + r0))^(6*v/c) - 4*(sin(theta_1)^4 - 2*(cos(theta_1) + 1)*sin(theta_1)^2)*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c) + 8*cos(theta_1) + 8)

# g21 =  -8*(((c*tau2 + (c*tau2 + r0)*cos(theta_1) + r0)*(-(cos(theta_1) - 1)/(cos(theta_1) + 1))^(3/2) - (c*tau2 - (c*tau2 + r0)*cos(theta_1) + r0)*sqrt(-(cos(theta_1) - 1)/(cos(theta_1) + 1)))*((c*tau1 + r0)/(c*tau2 + r0))^(6*v/c)*sin(theta_1)^2 + (((3*c*tau2 + (c*tau2 + r0)*cos(theta_1) + 3*r0)*(-(cos(theta_1) - 1)/(cos(theta_1) + 1))^(3/2) + (c*tau2 + (c*tau2 + r0)*cos(theta_1) + r0)*sqrt(-(cos(theta_1) - 1)/(cos(theta_1) + 1)))*sin(theta_1)^2 - 4*(c*tau2 + (c*tau2 + r0)*cos(theta_1) + r0)*(-(cos(theta_1) - 1)/(cos(theta_1) + 1))^(3/2))*((c*tau1 + r0)/(c*tau2 + r0))^(4*v/c))*v/(6*((c*tau1 + r0)/(c*tau2 + r0))^(4*v/c)*sin(theta_1)^4 + sin(theta_1)^4 - 4*(cos(theta_1) + 2)*sin(theta_1)^2 + (sin(theta_1)^4 + 4*(cos(theta_1) - 2)*sin(theta_1)^2 - 8*cos(theta_1) + 8)*((c*tau1 + r0)/(c*tau2 + r0))^(8*v/c) - 4*(sin(theta_1)^4 + 2*(cos(theta_1) - 1)*sin(theta_1)^2)*((c*tau1 + r0)/(c*tau2 + r0))^(6*v/c) - 4*(sin(theta_1)^4 - 2*(cos(theta_1) + 1)*sin(theta_1)^2)*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c) + 8*cos(theta_1) + 8)

# длина дуги координатной кривой tau2 в пределах theta_1 от theta_min до theta_max
theta_min = var("theta_min")
theta_max = var("theta_max")
assume(theta_max-theta_min>0)
assume(pi-theta_max>0)

L = integral( sqrt(g22), (theta_1, theta_min, theta_max) )
print "\nL = ", L

# L =  2*c*tau2*arctan((c*tau1 + r0)^(v/c)*sin(theta_max)/((c*tau2 + r0)^(v/c)*(cos(theta_max) + 1))) - 2*c*tau2*arctan((c*tau1 + r0)^(v/c)*sin(theta_min)/((c*tau2 + r0)^(v/c)*(cos(theta_min) + 1))) + 2*r0*arctan((c*tau1 + r0)^(v/c)*sin(theta_max)/((c*tau2 + r0)^(v/c)*(cos(theta_max) + 1))) - 2*r0*arctan((c*tau1 + r0)^(v/c)*sin(theta_min)/((c*tau2 + r0)^(v/c)*(cos(theta_min) + 1)))

# отношение длины дуги координатной кривой tau2 в пределах theta_1 от theta_min до theta_max к запаздывающему радиусу и к разности коррдинатных углов (theta_max-theta_min)
p = L/((c*tau2 + r0)*(theta_max-theta_min))
print "\np = ", p

# p =  2*(c*tau2*arctan((c*tau1 + r0)^(v/c)*sin(theta_max)/((c*tau2 + r0)^(v/c)*(cos(theta_max) + 1))) - c*tau2*arctan((c*tau1 + r0)^(v/c)*sin(theta_min)/((c*tau2 + r0)^(v/c)*(cos(theta_min) + 1))) + r0*arctan((c*tau1 + r0)^(v/c)*sin(theta_max)/((c*tau2 + r0)^(v/c)*(cos(theta_max) + 1))) - r0*arctan((c*tau1 + r0)^(v/c)*sin(theta_min)/((c*tau2 + r0)^(v/c)*(cos(theta_min) + 1))))/((c*tau2 + r0)*(theta_max - theta_min))

p = p.full_simplify()
print "\np = ", p

# p =  2*(arctan((c*tau1 + r0)^(v/c)*sin(theta_max)/((c*tau2 + r0)^(v/c)*(cos(theta_max) + 1))) - arctan((c*tau1 + r0)^(v/c)*sin(theta_min)/((c*tau2 + r0)^(v/c)*(cos(theta_min) + 1))))/(theta_max - theta_min)


# вычисление площади криволинейной клетки образованной координатными кривыми


tau_min = var("tau_min")
tau_max = var("tau_max")

assume(tau_max-tau_min>0)


from sage.symbolic.integration.integral import definite_integral
A = definite_integral(J, tau2, tau_min, tau_max)
print "\nA = ", A

from sage.symbolic.integration.integral import indefinite_integral
A = indefinite_integral(J, tau2)
print "\nA = ", A

A = integral (J, (tau2, tau_min, tau_max) )
print "\nA = ", A

T = integral( J, (theta_1, theta_min, theta_max) )
print "\nT = ", T

# T =  -2*((2*c^2*tau2 - (c^2*tau2 + c*r0)*sin(theta_max)^2 + 2*c*r0 + 2*(c^2*tau2 + c*r0)*cos(theta_max))*arctan(((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*sin(theta_max)/(cos(theta_max) + 1))*cos(theta_min)^2 + 2*(2*c^2*tau2 - (c^2*tau2 + c*r0)*sin(theta_max)^2 + 2*c*r0 + 2*(c^2*tau2 + c*r0)*cos(theta_max))*arctan(((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*sin(theta_max)/(cos(theta_max) + 1))*cos(theta_min) - (((c*tau2 + (c*tau2 + r0)*cos(theta_max) + r0)*sin(theta_max)*sin(theta_min)^2 - ((c*tau2 + r0)*cos(theta_min)*sin(theta_max)^2 + (c*tau2 + r0)*sin(theta_max)^2)*sin(theta_min))*((c*tau1 + r0)/(c*tau2 + r0))^(3*v/c) + ((c*tau2 + (c*tau2 + r0)*cos(theta_max) + r0)*cos(theta_min)^2*sin(theta_max) + 2*(c*tau2 + (c*tau2 + r0)*cos(theta_max) + r0)*cos(theta_min)*sin(theta_max) + (c*tau2 + (c*tau2 + r0)*cos(theta_max) + r0)*sin(theta_max) - ((c*tau2 + r0)*cos(theta_max)^2 + c*tau2 + 2*(c*tau2 + r0)*cos(theta_max) + ((c*tau2 + r0)*cos(theta_max)^2 + c*tau2 + 2*(c*tau2 + r0)*cos(theta_max) + r0)*cos(theta_min) + r0)*sin(theta_min))*((c*tau1 + r0)/(c*tau2 + r0))^(v/c))*v + ((c^2*tau2 + c*r0)*arctan(((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*sin(theta_max)/(cos(theta_max) + 1))*sin(theta_max)^2*sin(theta_min)^2 - (c^2*tau2 + c*r0)*arctan(((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*sin(theta_min)/(cos(theta_min) + 1))*sin(theta_max)^2*sin(theta_min)^2)*((c*tau1 + r0)/(c*tau2 + r0))^(4*v/c) + ((c^2*tau2 + c*r0)*arctan(((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*sin(theta_max)/(cos(theta_max) + 1))*cos(theta_min)^2*sin(theta_max)^2 + 2*(c^2*tau2 + c*r0)*arctan(((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*sin(theta_max)/(cos(theta_max) + 1))*cos(theta_min)*sin(theta_max)^2 + (c^2*tau2 + c*r0)*arctan(((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*sin(theta_max)/(cos(theta_max) + 1))*sin(theta_max)^2 + (2*c^2*tau2 - (c^2*tau2 + c*r0)*sin(theta_max)^2 + 2*c*r0 + 2*(c^2*tau2 + c*r0)*cos(theta_max))*arctan(((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*sin(theta_max)/(cos(theta_max) + 1))*sin(theta_min)^2 - (2*(c^2*tau2 + c*r0)*cos(theta_min)*sin(theta_max)^2 + 2*(c^2*tau2 + c*r0)*sin(theta_max)^2 + (c^2*tau2 + (c^2*tau2 + c*r0)*cos(theta_max)^2 - (c^2*tau2 + c*r0)*sin(theta_max)^2 + c*r0 + 2*(c^2*tau2 + c*r0)*cos(theta_max))*sin(theta_min)^2)*arctan(((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*sin(theta_min)/(cos(theta_min) + 1)))*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c) + (2*c^2*tau2 - (c^2*tau2 + c*r0)*sin(theta_max)^2 + 2*c*r0 + 2*(c^2*tau2 + c*r0)*cos(theta_max))*arctan(((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*sin(theta_max)/(cos(theta_max) + 1)) - (2*c^2*tau2 + 2*(c^2*tau2 + c*r0)*cos(theta_max)^2 - (c^2*tau2 + (c^2*tau2 + c*r0)*cos(theta_max)^2 + c*r0 + 2*(c^2*tau2 + c*r0)*cos(theta_max))*sin(theta_min)^2 + 2*c*r0 + 4*(c^2*tau2 + c*r0)*cos(theta_max) + 2*(c^2*tau2 + (c^2*tau2 + c*r0)*cos(theta_max)^2 + c*r0 + 2*(c^2*tau2 + c*r0)*cos(theta_max))*cos(theta_min))*arctan(((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*sin(theta_min)/(cos(theta_min) + 1)))/(((c*tau1 + r0)/(c*tau2 + r0))^(4*v/c)*sin(theta_max)^2*sin(theta_min)^2 + (cos(theta_max)^2 + 2*cos(theta_max) + 1)*cos(theta_min)^2 + (cos(theta_min)^2*sin(theta_max)^2 + 2*cos(theta_min)*sin(theta_max)^2 + (cos(theta_max)^2 + 2*cos(theta_max) + 1)*sin(theta_min)^2 + sin(theta_max)^2)*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c) + cos(theta_max)^2 + 2*(cos(theta_max)^2 + 2*cos(theta_max) + 1)*cos(theta_min) + 2*cos(theta_max) + 1)

T = T.full_simplify()
print "\nT = ", T

# T =  2*((2*c^2*tau2 - (c^2*tau2 + c*r0)*sin(theta_max)^2 + 2*c*r0 + 2*(c^2*tau2 + c*r0)*cos(theta_max))*arctan(((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*sin(theta_max)/(cos(theta_max) + 1))*sin(theta_min)^2 - 2*(2*c^2*tau2 - (c^2*tau2 + c*r0)*sin(theta_max)^2 + 2*c*r0 + 2*(c^2*tau2 + c*r0)*cos(theta_max))*arctan(((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*sin(theta_max)/(cos(theta_max) + 1))*cos(theta_min) + (((c*tau2 + (c*tau2 + r0)*cos(theta_max) + r0)*sin(theta_max)*sin(theta_min)^2 - ((c*tau2 + r0)*cos(theta_min)*sin(theta_max)^2 + (c*tau2 + r0)*sin(theta_max)^2)*sin(theta_min))*((c*tau1 + r0)/(c*tau2 + r0))^(3*v/c) - ((c*tau2 + (c*tau2 + r0)*cos(theta_max) + r0)*sin(theta_max)*sin(theta_min)^2 - 2*(c*tau2 + (c*tau2 + r0)*cos(theta_max) + r0)*cos(theta_min)*sin(theta_max) - 2*(c*tau2 + (c*tau2 + r0)*cos(theta_max) + r0)*sin(theta_max) - ((c*tau2 + r0)*sin(theta_max)^2 - 2*c*tau2 - 2*(c*tau2 + r0)*cos(theta_max) + ((c*tau2 + r0)*sin(theta_max)^2 - 2*c*tau2 - 2*(c*tau2 + r0)*cos(theta_max) - 2*r0)*cos(theta_min) - 2*r0)*sin(theta_min))*((c*tau1 + r0)/(c*tau2 + r0))^(v/c))*v - ((c^2*tau2 + c*r0)*arctan(((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*sin(theta_max)/(cos(theta_max) + 1))*sin(theta_max)^2*sin(theta_min)^2 - (c^2*tau2 + c*r0)*arctan(((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*sin(theta_min)/(cos(theta_min) + 1))*sin(theta_max)^2*sin(theta_min)^2)*((c*tau1 + r0)/(c*tau2 + r0))^(4*v/c) - 2*((c^2*tau2 + c*r0)*arctan(((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*sin(theta_max)/(cos(theta_max) + 1))*cos(theta_min)*sin(theta_max)^2 + (c^2*tau2 + c*r0)*arctan(((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*sin(theta_max)/(cos(theta_max) + 1))*sin(theta_max)^2 + (c^2*tau2 - (c^2*tau2 + c*r0)*sin(theta_max)^2 + c*r0 + (c^2*tau2 + c*r0)*cos(theta_max))*arctan(((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*sin(theta_max)/(cos(theta_max) + 1))*sin(theta_min)^2 - ((c^2*tau2 + c*r0)*cos(theta_min)*sin(theta_max)^2 + (c^2*tau2 + c*r0)*sin(theta_max)^2 + (c^2*tau2 - (c^2*tau2 + c*r0)*sin(theta_max)^2 + c*r0 + (c^2*tau2 + c*r0)*cos(theta_max))*sin(theta_min)^2)*arctan(((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*sin(theta_min)/(cos(theta_min) + 1)))*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c) - 2*(2*c^2*tau2 - (c^2*tau2 + c*r0)*sin(theta_max)^2 + 2*c*r0 + 2*(c^2*tau2 + c*r0)*cos(theta_max))*arctan(((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*sin(theta_max)/(cos(theta_max) + 1)) + (4*c^2*tau2 - 2*(c^2*tau2 + c*r0)*sin(theta_max)^2 - (2*c^2*tau2 - (c^2*tau2 + c*r0)*sin(theta_max)^2 + 2*c*r0 + 2*(c^2*tau2 + c*r0)*cos(theta_max))*sin(theta_min)^2 + 4*c*r0 + 4*(c^2*tau2 + c*r0)*cos(theta_max) + 2*(2*c^2*tau2 - (c^2*tau2 + c*r0)*sin(theta_max)^2 + 2*c*r0 + 2*(c^2*tau2 + c*r0)*cos(theta_max))*cos(theta_min))*arctan(((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*sin(theta_min)/(cos(theta_min) + 1)))/(((c*tau1 + r0)/(c*tau2 + r0))^(4*v/c)*sin(theta_max)^2*sin(theta_min)^2 + (cos(theta_max)^2 + 2*cos(theta_max) + 1)*cos(theta_min)^2 + (cos(theta_min)^2*sin(theta_max)^2 + 2*cos(theta_min)*sin(theta_max)^2 + (cos(theta_max)^2 + 2*cos(theta_max) + 1)*sin(theta_min)^2 + sin(theta_max)^2)*((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c) + cos(theta_max)^2 + 2*(cos(theta_max)^2 + 2*cos(theta_max) + 1)*cos(theta_min) + 2*cos(theta_max) + 1)

#S = integral (T, tau2 )
#print "\nS = ", S


#S = integral (T, (tau2, tau_min, tau_max) )
#print "\nS = ", S


# график кривых распространения электрического потенциала движущегося заряда и изоповерхностей запаздывающего момента
from sage.plot.circle import Circle
g = Graphics()

dt = 0.2
r0 = 0.2
c = 1
v = 0.9
tau1 = 5

g += circle((0, 0), r0, rgbcolor=hue(r0/500))

for theta_i in np.arange(0, 360, 15):
    theta_1 = theta_i * pi / 180.0
    plot_data = []
    for tau_i in np.arange(0, tau1-r0/c+dt, dt):
        tau2 = tau_i
        x2 =  2*(c*tau2 + r0)*((c*tau1 + r0)/(c*tau2 + r0))^(v/c)*tan(1/2*theta_1)/(((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)
        z2 =  -tau2*v - (((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 - 1)*(c*tau2 + r0)/(((c*tau1 + r0)/(c*tau2 + r0))^(2*v/c)*tan(1/2*theta_1)^2 + 1)
        
        plot_data += [(x2, z2)]
        
        if (dt/2 > (tau2+r0/c) % 1 and theta_i == 0):
            zq = tau2*v
            rr = c*tau2 + r0
            g += circle((0, -zq), rr, rgbcolor=hue(rr/5))
        
    g += list_plot(plot_data,size=2)

pname = "results/field_of_directions_of_propagation_of_electric_potential.png"
print pname
g.save(pname)


# график отношения длины дуги координатной кривой tau2 в пределах theta_1 от theta_min до theta_max к запаздывающему радиусу и к разности коррдинатных углов (theta_max-theta_min)
g2 = Graphics()

theta_min = 0
for theta_i in np.arange(5, 180-5, 5):
    theta_max = theta_i * pi / 180.0
    plot_data = []
    for tau_i in np.arange(0, tau1-r0/c+dt, dt):
        tau2 = tau_i
        print "\ntau2 = ", tau2
        print "\ntheta_min = ", theta_min
        print "\ntheta_max = ", theta_max

        p =  2*(arctan((c*tau1 + r0)^(v/c)*sin(theta_max)/((c*tau2 + r0)^(v/c)*(cos(theta_max) + 1))) - arctan((c*tau1 + r0)^(v/c)*sin(theta_min)/((c*tau2 + r0)^(v/c)*(cos(theta_min) + 1))))/(theta_max - theta_min)
        print "\np = ", p
        p = p.n()
        print "\np = ", p

        plot_data += [(tau2, p)]

    g2 += list_plot(plot_data,size=2)
    theta_min = theta_max

pname = "results/relation_L_per_R_and_delta_theta.png"
print pname
g2.save(pname)


