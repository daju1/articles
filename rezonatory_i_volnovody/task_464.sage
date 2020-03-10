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


# rot E = i * omega / c * H
# rot H = - i * epsilon * omega / c * epsilon * E
# mu = 1
# Laplasian E + epsilon * omega^2 / c^2 * E - grad div E = 0
# Laplasian H + epsilon * omega^2 / c^2 * H + 1/epsilon * [ grad Epsilon x rot H] = 0

# уравнение Даламбера
# mu * epsilon / c^2 * diff(E, t, 2) - Laplasian E = 0

# Ищем решение в виде
# E = E(x) * exp(I * k * z) * exp(-I * omega * t)

# подставляя в уравнение Даламбера приходим к обыкновенному дифференциальному уравнению
# diff (E(x), x, 2) + (mu * epsilon * omega^2 / c^2 - k^2) * E(x) = 0

# для диэлектрика
# kappa^2 = (mu * epsilon * omega^2 / c^2 - k^2)
# diff (E(x), x, 2) + kappa^2 * E(x) = 0

# для вакуума с граничными условиями на бесконечности
# - s^2 = (mu * epsilon * omega^2 / c^2 - k^2)
# diff (E(x), x, 2) - s^2 * E(x) = 0

# исходя из вида дифференциальных уравнений для диэлектрика и вакуума
# записываем соотнощения для волнового вектора в виде уравнений
eq_kappa = mu * epsilon * omega^2 / c^2 - k^2 == kappa^2
# eq_s     = omega^2 / c^2 - k^2 == - s^2
eq_s     = omega^2 / c^2 + s^2 == k^2

# решаем полученные обыкновенные дифференциальные уравнения
E = function('E')(x)

de1 = diff (E, x, 2) + kappa^2 * E == 0
E_dielectric = desolve(de1, E, ivar = x)
de2 = diff (E, x, 2) - s^2 * E == 0
E_vacuum = desolve(de2, E, ivar = x)

# получая таким образом выражения для электрического поля в диэлектрике и в вакууме
print E_dielectric
# _K2*cos(kappa*x) + _K1*sin(kappa*x)
print E_vacuum
# _K1*e^(s*x) + _K2*e^(-s*x)

_B1_z = var ('_B1_z')
_B2_z = var ('_B2_z')
_A_z = var ('_A_z')

_B1_x = var ('_B1_x')
_B2_x = var ('_B2_x')
_A_x = var ('_A_x')

# E_dielectric_z = E_dielectric.subs(var('_K1'), _B1_z).subs(var('_K2'), _B2_z)

# выражения для электрического поля в диэлектрике и в вакууме
E_dielectric_z (x) = _B1_z*cos(kappa*x) + _B2_z*sin(kappa*x)
E_dielectric_x (x) = _B1_x*cos(kappa*x) + _B2_x*sin(kappa*x)

# x > a
E_vacuum2_z (x) = _A_z*e^(-s*x)
E_vacuum2_x (x) = _A_x*e^(-s*x)

# x < -a
E_vacuum1_z (x) = -_A_z*e^(s*x)
E_vacuum1_x (x) = _A_x*e^(s*x)

# выражения для магнитного поля в диэлектрике и в вакууме
# H_y = c / (I * omega) * (diff(E_x, z) - diff(E_z, x))

H_dielectric_y (x) = c / (I * omega) * (I*k*E_dielectric_x (x) - diff(E_dielectric_z (x), x))
H_vacuum1_y (x) = c / (I * omega) * (I*k*E_vacuum1_x (x) - diff(E_vacuum1_z (x), x))
H_vacuum2_y (x) = c / (I * omega) * (I*k*E_vacuum2_x (x) - diff(E_vacuum2_z (x), x))

print "H_dielectric_y (x) ", H_dielectric_y (x)
print "H_vacuum1_y (x) =", H_vacuum1_y (x)
print "H_vacuum2_y (x) =", H_vacuum2_y (x)

# уравнение Максвелла для тока смещения
# rot_H = -I * epsilon * omega / c * E
# rot_H_x = - diff(H_y, z)
# rot_H_z = diff(H_y, x)

rot_H_dielectric_x = - I * k * H_dielectric_y (x)
rot_H_vacuum1_x    = - I * k * H_vacuum1_y (x)
rot_H_vacuum2_x    = - I * k * H_vacuum2_y (x)
rot_H_dielectric_z = diff(H_dielectric_y (x), x)
rot_H_vacuum1_z    = diff(H_vacuum1_y (x), x)
rot_H_vacuum2_z    = diff(H_vacuum2_y (x), x)

print "rot_H_dielectric_x =", rot_H_dielectric_x
print "rot_H_vacuum1_x =", rot_H_vacuum1_x
print "rot_H_vacuum2_x =", rot_H_vacuum2_x
print "rot_H_dielectric_z =", rot_H_dielectric_z
print "rot_H_vacuum1_z =", rot_H_vacuum1_z
print "rot_H_vacuum2_z =", rot_H_vacuum2_z

eq_rot_H_dielectric_x = rot_H_dielectric_x == -I * epsilon * omega / c * E_dielectric_x (x)
eq_rot_H_dielectric_z = rot_H_dielectric_z == -I * epsilon * omega / c * E_dielectric_z (x)
eq_rot_H_vacuum1_x     = rot_H_vacuum1_x     == -I * epsilon * omega / c * E_vacuum1_x (x)
eq_rot_H_vacuum2_x     = rot_H_vacuum2_x     == -I * epsilon * omega / c * E_vacuum2_x (x)
eq_rot_H_vacuum1_z     = rot_H_vacuum1_z     == -I * epsilon * omega / c * E_vacuum1_z (x)
eq_rot_H_vacuum2_z     = rot_H_vacuum2_z     == -I * epsilon * omega / c * E_vacuum2_z (x)

# граничные условия на гринице вакуума и диэлектрика
# тангенциальная компонента электрического поля
eqEz_a  = E_dielectric_z (a)  == E_vacuum2_z (a)
eqEz_ma = E_dielectric_z (-a) == E_vacuum1_z (-a)

print "eqEz_a = ",  eqEz_a
print "eqEz_ma = ",  eqEz_ma

# нормальная компонента электрического смещения
eqEx_a  = E_dielectric_x (a) * epsilon  == E_vacuum2_x (a)
eqEx_ma = E_dielectric_x (-a) * epsilon == E_vacuum1_x (-a)

print "eqEx_a =", eqEx_a
print "eqEx_ma =", eqEx_ma

# eqEz_a =  _B1_z*cos(a*kappa) + _B2_z*sin(a*kappa) == _A1_z*e^(a*s) + _A2_z*e^(-a*s)
# eqEz_ma =  _B1_z*cos(a*kappa) - _B2_z*sin(a*kappa) == _A2_z*e^(a*s) + _A1_z*e^(-a*s)
# eqEx_a = (_B1_x*cos(a*kappa) + _B2_x*sin(a*kappa))*epsilon == _A1_x*e^(a*s) + _A2_x*e^(-a*s)
# eqEx_ma = (_B1_x*cos(a*kappa) - _B2_x*sin(a*kappa))*epsilon == _A2_x*e^(a*s) + _A1_x*e^(-a*s)


# тангенциальная компонента напряженности магнитного поля
eqHy_a  = H_vacuum2_y (a)  == H_dielectric_y (a)
eqHy_ma = H_vacuum1_y (-a) == H_dielectric_y (-a)

print "eqHy_a =", eqHy_a
print "eqHy_ma =", eqHy_ma

# H_x = c / (I * omega) * (diff(E_z, y) - diff(E_y, z))
# H_z = c / (I * omega) * (diff(E_y, x) - diff(E_x, y))
# E_y = 0

# H_dielectric_x (x) = c / (I * omega) * (diff(E_dielectric_z (x), y))


res = solve([eq_rot_H_dielectric_x,
             eq_rot_H_dielectric_z,
             eq_rot_H_vacuum1_x,
             eq_rot_H_vacuum2_x,
             eq_rot_H_vacuum1_z,
             eq_rot_H_vacuum2_z,
             eqEz_a, eqEz_ma,
             eqEx_a, eqEx_ma,
             eqHy_a, eqHy_ma,
             eq_s, eq_kappa,
             _B1_z==0, _B2_x==0,
             s > 0, c > 0, a > 0, omega > 0
             ],
    #_B2_z,
    #_A_z,

    _B1_x,
    #_A_x,

    #kappa, s
    )
print "res =", res


res = solve([eqEz_a, eqEx_a, eqHy_a, _B1_z==0, _B2_x==0],
    #_B1_z,
    _B2_z,
    _A_z,

    _B1_x,
    _B2_x,
    _A_x,
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


res = solve([eqEz_a, eqEx_a, eqHy_a, _B1_z==0, _B2_x==0, eq_s, eq_kappa, s > 0, c > 0, a > 0, omega > 0],
    #_B1_z,
    _B2_z,
    _A_z,

    _B1_x,
    _B2_x,
    _A_x,
    kappa, s
    )
print "res =", res

'''
res = [
[
    I*_B2_x*k*cos(a*kappa)*e^(a*s) - _B1_z*kappa*cos(a*kappa)*e^(a*s) - I*_A2_x*k - _A2_z*s == 0,
    -c^2*k^2 + c^2*s^2 + omega^2 == 0,
    -c^2*k^2 - c^2*kappa^2 + epsilon*mu*omega^2 == 0,
    -_B1_z*e^(a*s)*sin(a*kappa) + _A2_z == 0,
    -_B2_x*epsilon*cos(a*kappa)*e^(a*s) + _A2_x == 0
]
]
'''

M = matrix(SR, 6, 6,
[
    [ epsilon*cos(a*kappa),  epsilon*sin(a*kappa),                   0,                   0,     -e^(-a*s),           0 ],
    [                    0,                     0,        cos(a*kappa),        sin(a*kappa),             0,   -e^(-a*s) ],
    [     I*k*cos(a*kappa),      I*k*sin(a*kappa),  kappa*sin(a*kappa), -kappa*cos(a*kappa), -I*k*e^(-a*s), -s*e^(-a*s) ],

    [ epsilon*cos(a*kappa), -epsilon*sin(a*kappa),                   0,                   0,     -e^(-a*s),           0 ],
    [                    0,                     0,        cos(a*kappa),       -sin(a*kappa),             0,   +e^(-a*s) ],
    [     I*k*cos(a*kappa),     -I*k*sin(a*kappa), -kappa*sin(a*kappa), -kappa*cos(a*kappa), -I*k*e^(-a*s), -s*e^(-a*s) ],

])

print M

mdet = M.det()
print mdet
# 8*I*epsilon^2*k*s*cos(a*kappa)^2*e^(-2*a*s)*sin(a*kappa)^2 - 8*I*epsilon^2*k*kappa*cos(a*kappa)*e^(-2*a*s)*sin(a*kappa)^3 - 8*I*epsilon*k*s*cos(a*kappa)^2*e^(-2*a*s)*sin(a*kappa)^2 + 8*I*epsilon*k*kappa*cos(a*kappa)*e^(-2*a*s)*sin(a*kappa)^3

print solve([mdet == 0], s)
# s == kappa*sin(a*kappa)/cos(a*kappa)
#it seems wrong

M = matrix(SR, 4, 4,
[
    [ epsilon*cos(a*kappa),                     0,     -e^(-a*s),           0 ],
    [                    0,          sin(a*kappa),             0,   -e^(-a*s) ],
    [     I*k*cos(a*kappa),   -kappa*cos(a*kappa), -I*k*e^(-a*s), -s*e^(-a*s) ],
    [                kappa,                  -I*k,             0,            0]
])
#                       Bx,                    Bx,            Ax,            Az

mdet = M.det()
print mdet
# -epsilon*k^2*cos(a*kappa)*e^(-2*a*s) + k^2*cos(a*kappa)*e^(-2*a*s) + kappa^2*cos(a*kappa)*e^(-2*a*s) + kappa*s*e^(-2*a*s)*sin(a*kappa)

print solve([mdet == 0], s)

# s == ((epsilon - 1)*k^2 - kappa^2)*cos(a*kappa)/(kappa*sin(a*kappa))


