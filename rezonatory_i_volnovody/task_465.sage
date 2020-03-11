import sys
reload(sys)
sys.setdefaultencoding('utf8')


mu = var('mu')
assume(mu>0)
epsilon = var('epsilon')
assume(epsilon>0)

s = var ('s')
assume(s>0)

kappa = var ('kappa')
assume(kappa>0)

k = var ('k')

a = var('a')
assume(a>0)

x = var('x')
y = var('y')
z = var('z')

c = var('c')
assume(c>0)
omega = var('omega')
sigma = var('sigma')


# rot E = i * mu * omega / c * H
# rot H = 4*pi*j - i * epsilon * omega / c * epsilon * E

# Laplasian E - grad div E + mu * epsilon * omega^2 / c^2 * E  = - mu * i * omega / c^2 * 4*pi*j
# в металле
# j = E * sigma
# Laplasian E - grad div E + mu * epsilon * omega^2 / c^2 * E  + mu * i * omega / c^2 * 4*pi * sigma * E = 0


# уравнение Даламбера
# mu * epsilon / c^2 * diff(E, t, 2) - Laplasian E = 0
# Laplasian E + (mu * epsilon * omega^2 / c^2  + mu * i * omega / c^2 * 4 * pi * sigma ) * E = 0

# Ищем решение в виде
# E = E(x) * exp(I * k * z) * exp(-I * omega * t)

# подставляя в уравнение Даламбера приходим к обыкновенному дифференциальному уравнению
# diff (E(x), x, 2) + (mu * epsilon * omega^2 / c^2 + mu * i * omega / c^2 * 4 * pi * sigma - k^2) * E(x) = 0

# для диэлектрика
# kappa^2 = (mu * epsilon * omega^2 / c^2 - k^2)
# diff (E(x), x, 2) + kappa^2 * E(x) = 0

# исходя из вида дифференциальных уравнений для диэлектрика и металла
# записываем соотнощения для волнового вектора в виде уравнений
eq_kappa = mu * epsilon * omega^2 / c^2 - k^2 == kappa^2
# eq_s     = omega^2 / c^2 - k^2 == - s^2
eq_s     = omega^2 / c^2 + s^2 == k^2

# решаем полученные обыкновенные дифференциальные уравнения
E = function('E')(x)

de_dielectric = diff (E, x, 2) + kappa^2 * E == 0
E_dielectric = desolve(de_dielectric, E, ivar = x)

de_metall  = diff (E, x, 2) + (mu * epsilon * omega^2 / c^2 + mu * I * omega / c^2 * 4 * pi * sigma - k^2) * E == 0
E_metall = desolve(de_metall, E, ivar = x)

# получая таким образом выражения для электрического поля в диэлектрике и в металле
print E_dielectric
# _K2*cos(kappa*x) + _K1*sin(kappa*x)
print E_metall
# _K2*cos(1/2*sqrt(-4*c^2*k^2 + 4*epsilon*mu*omega^2 + 16*I*pi*mu*omega*sigma)*x/c) + _K1*sin(1/2*sqrt(-4*c^2*k^2 + 4*epsilon*mu*omega^2 + 16*I*pi*mu*omega*sigma)*x/c)