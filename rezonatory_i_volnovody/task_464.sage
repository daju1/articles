import sys
reload(sys)
sys.setdefaultencoding('utf8')

# rot E = i * omega / c * H
# rot H = - i * epsilon * omega / c * epsilon * E
# mu = 1
# Laplasian E + epsilon * omega^2 / c^2 * E - grad div E = 0
# Laplasian H + epsilon * omega^2 / c^2 * H + 1/epsilon * [ grad Epsilon x rot H] = 0

# уравнение Даламбера
# mu * epsilon / c^2 * diff(E, t, 2) - Laplasian E = 0

# Ищем решение в виде
# E = E(x) * exp(k * z)

# подставляя в уравнение Даламбера приходим к обыкновенному дифференциальному уравнению

# diff (E(x), x, 2) + (mu * epsilon * omega^2 / c^2 - k^2) * E(x) = 0

# для диэлектрика
# kappa^2 = (mu * epsilon * omega^2 / c^2 - k^2)
# diff (E(x), x, 2) + kappa^2 * E(x) = 0

# для вакуума с граничными условиями на бесконечности
# - s^2 = (mu * epsilon * omega^2 / c^2 - k^2)
# diff (E(x), x, 2) - s^2 * E(x) = 0

E = function('E')(x)
kappa = var ('kappa')
assume(kappa>0)
x = var('x')
z = var('z')

de1 = diff (E, x, 2) + kappa^2 * E == 0
E_dielectric = desolve(de1, E, ivar = x)
print E_dielectric
# _K2*cos(kappa*x) + _K1*sin(kappa*x)

k = var ('k')
s = var ('s')
assume(s>0)
de2 = diff (E, x, 2) - s^2 * E == 0
E_vacuum = desolve(de2, E, ivar = x)
print E_vacuum
# _K1*e^(s*x) + _K2*e^(-s*x)



a = var('a')
epsilon = var('epsilon')
_B1_z = var ('_B1_z')
_B2_z = var ('_B2_z')
_A1_z = var ('_A1_z')
_A2_z = var ('_A2_z')

_B1_x = var ('_B1_x')
_B2_x = var ('_B2_x')
_A1_x = var ('_A1_x')
_A2_x = var ('_A2_x')

# E_dielectric_z = E_dielectric.subs(var('_K1'), _B1_z).subs(var('_K2'), _B2_z)

E_dielectric_z (x) = _B2_z*cos(kappa*x) + _B1_z*sin(kappa*x)
E_dielectric_x (x) = _B2_x*cos(kappa*x) + _B1_x*sin(kappa*x)

E_vacuum_z (x) = _A1_z*e^(s*x) + _A2_z*e^(-s*x)
E_vacuum_x (x) = _A1_x*e^(s*x) + _A2_x*e^(-s*x)


eqEz_a = E_dielectric_z (a) == E_vacuum_z (a)
print "eqEz_a = ",  eqEz_a
eqEz_ma = E_dielectric_z (-a) == E_vacuum_z (-a)
print "eqEz_ma = ",  eqEz_ma

# eqEz_a  =  _B2_z*cos(a*kappa)  + _B1_z*sin(a*kappa)  == _A1_z*e^(a*s) + _A2_z*e^(-a*s)
# eqEz_ma =  _B2_z*cos(-a*kappa) + _B1_z*sin(-a*kappa) == _A2_z*e^(a*s) + _A1_z*e^(-a*s)

eqEx_a = E_dielectric_x (a) * epsilon == E_vacuum_x (a)
print "eqEx_a =", eqEx_a
eqEx_ma = E_dielectric_x (-a) * epsilon == E_vacuum_x (-a)
print "eqEx_ma =", eqEx_ma

# eq3 = epsilon * (_B2_x*cos(kappa*a)    + _B1_x*sin(kappa*a))    == _A1_x*e^(s*a)    + _A2_x*e^(-s*a)
# eq4 = epsilon * (_B2_x*cos(kappa*(-a)) + _B1_x*sin(kappa*(-a))) == _A1_x*e^(s*(-a)) + _A2_x*e^(-s*(-a))
# eqEx_a  = (_B2_x*cos(a*kappa)  + _B1_x*sin(a*kappa))  * epsilon == _A1_x*e^(a*s) + _A2_x*e^(-a*s)
# eqEx_ma = (_B2_x*cos(-a*kappa) + _B1_x*sin(-a*kappa)) * epsilon == _A2_x*e^(a*s) + _A1_x*e^(-a*s)

c = var('c')
omega = var('omega')
# H_y = c / (I * omega) * (diff(E_x, z) - diff(E_z, x))

H_dielectric_y (x) = c / (I * omega) * (I*k*E_dielectric_x (x) - diff(E_dielectric_z (x), x))
print "H_dielectric_y (x) ", H_dielectric_y (x)

H_vacuum_y (x) = c / (I * omega) * (I*k*E_vacuum_x (x) - diff(E_vacuum_z (x), x))
print "H_vacuum_y (x) =", H_vacuum_y (x)

eqHy_a = H_vacuum_y (a) == H_dielectric_y (a)
print "eqHy_a =", eqHy_a

eqHy_ma = H_vacuum_y (-a) == H_dielectric_y (-a)
print "eqHy_ma =", eqHy_ma

# H_x = c / (I * omega) * (diff(E_z, y) - diff(E_y, z))
# H_z = c / (I * omega) * (diff(E_y, x) - diff(E_x, y))
# E_y = 0

# H_dielectric_x (x) = c / (I * omega) * (diff(E_dielectric_z (x), y))

res = solve([eqEz_a, eqEx_a, eqHy_a, _B2_z==0, _B1_x==0, _A1_z==0, _A1_x==0],
    #_B1_z,
    _B2_z,
    _A1_z,
    _A2_z,

    _B1_x,
    _B2_x,
    _A1_x,
    _A2_x,
    # kappa, s
    )
print "res =", res

'''
res = [
[_B2_z == 0,
 _A1_z == 0,
 _A2_z == _B1_z*e^(a*s)*sin(a*kappa),
 _B1_x == 0,
 _B2_x == (I*_B1_z*kappa*cos(a*kappa) + I*_B1_z*s*sin(a*kappa))/((epsilon*k - k)*cos(a*kappa)),
 _A1_x == 0,
 _A2_x == (I*_B1_z*epsilon*kappa*cos(a*kappa) + I*_B1_z*epsilon*s*sin(a*kappa))*e^(a*s)/(epsilon*k - k)]
]
'''

