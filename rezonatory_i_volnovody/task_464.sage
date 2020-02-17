import sys
reload(sys)
sys.setdefaultencoding('utf8')

# rot E = i * omega / c * H
# rot H = - i * epsilon * omega / c * epsilon * E
# mu = 1
# Laplasian E + epsilon * omega^2 / c^2 * E - grad div E = 0
# Laplasian H 

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

de1 = diff (E, x, 2) + kappa^2 * E == 0
res1 = desolve(de1, E, ivar = x)
print res1
# _K2*cos(kappa*x) + _K1*sin(kappa*x)

s = var ('s')
assume(s>0)
de2 = diff (E, x, 2) - s^2 * E == 0
res2 = desolve(de2, E, ivar = x)
print res2
# _K1*e^(s*x) + _K2*e^(-s*x)



a = var('a')
a = var('epsilon')
_B1_z = var ('_B1_z')
_B2_z = var ('_B2_z')
_A1_z = var ('_A1_z')
_A2_z = var ('_A2_z')

_B1_x = var ('_B1_x')
_B2_x = var ('_B2_x')
_A1_x = var ('_A1_x')
_A2_x = var ('_A2_x')

eq1 = _B2_z*cos(kappa*a)    + _B1_z*sin(kappa*a)    == _A1_z*e^(s*a)    + _A2_z*e^(-s*a)
eq2 = _B2_z*cos(kappa*(-a)) + _B1_z*sin(kappa*(-a)) == _A1_z*e^(s*(-a)) + _A2_z*e^(-s*(-a))

eq3 = epsilon * (_B2_x*cos(kappa*a)    + _B1_x*sin(kappa*a))    == _A1_x*e^(s*a)    + _A2_x*e^(-s*a)
eq4 = epsilon * (_B2_x*cos(kappa*(-a)) + _B1_x*sin(kappa*(-a))) == _A1_x*e^(s*(-a)) + _A2_x*e^(-s*(-a))



