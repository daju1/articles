import sys
try:
    reload(sys)
    sys.setdefaultencoding('utf8')
except:
    pass


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

de_dielectric = diff (E, x, 2) + kappa^2 * E == 0
E_dielectric = desolve(de_dielectric, E, ivar = x)
de_vacuum  = diff (E, x, 2) - s^2 * E == 0
E_vacuum = desolve(de_vacuum, E, ivar = x)

# получая таким образом выражения для электрического поля в диэлектрике и в вакууме
print (E_dielectric)
# _K2*cos(kappa*x) + _K1*sin(kappa*x)
print (E_vacuum)
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
# H_y = c / (I * omega * mu) * (diff(E_x, z) - diff(E_z, x))

H_dielectric_y (x) = c / (I * omega * mu) * (I*k*E_dielectric_x (x) - diff(E_dielectric_z (x), x))
H_vacuum1_y (x) = c / (I * omega) * (I*k*E_vacuum1_x (x) - diff(E_vacuum1_z (x), x))
H_vacuum2_y (x) = c / (I * omega) * (I*k*E_vacuum2_x (x) - diff(E_vacuum2_z (x), x))

H_dielectric_y (x) = c / ( omega * mu) * (k*E_dielectric_x (x) + I * diff(E_dielectric_z (x), x))
H_vacuum1_y (x) = c / ( omega) * (k*E_vacuum1_x (x) + I * diff(E_vacuum1_z (x), x))
H_vacuum2_y (x) = c / ( omega) * (k*E_vacuum2_x (x) + I * diff(E_vacuum2_z (x), x))

print ("H_dielectric_y (x) ", H_dielectric_y (x).simplify())
print ("H_vacuum1_y (x) =", H_vacuum1_y (x))
print ("H_vacuum2_y (x) =", H_vacuum2_y (x))

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

print ("rot_H_dielectric_x =", rot_H_dielectric_x)
print ("rot_H_vacuum1_x =", rot_H_vacuum1_x)
print ("rot_H_vacuum2_x =", rot_H_vacuum2_x)
print ("rot_H_dielectric_z =", rot_H_dielectric_z)
print ("rot_H_vacuum1_z =", rot_H_vacuum1_z)
print ("rot_H_vacuum2_z =", rot_H_vacuum2_z)

eq_rot_H_dielectric_x = rot_H_dielectric_x == -I * epsilon * omega / c * E_dielectric_x (x)
eq_rot_H_dielectric_z = rot_H_dielectric_z == -I * epsilon * omega / c * E_dielectric_z (x)
eq_rot_H_vacuum1_x     = rot_H_vacuum1_x     == -I * omega / c * E_vacuum1_x (x)
eq_rot_H_vacuum2_x     = rot_H_vacuum2_x     == -I * omega / c * E_vacuum2_x (x)
eq_rot_H_vacuum1_z     = rot_H_vacuum1_z     == -I * omega / c * E_vacuum1_z (x)
eq_rot_H_vacuum2_z     = rot_H_vacuum2_z     == -I * omega / c * E_vacuum2_z (x)

print ("eq_rot_H_dielectric_x =", eq_rot_H_dielectric_x.simplify())
print ("eq_rot_H_dielectric_z =", eq_rot_H_dielectric_z.simplify())


res = solve( eq_rot_H_dielectric_x.subs(_B1_z==0, _B2_x==0), _B1_x)
print ("res =", res)

res = solve( eq_kappa, kappa)
print ("res =", res)


res = solve( eq_rot_H_dielectric_z.subs(_B1_z==0, _B2_x==0), _B1_x)
print ("res =", res)

# граничные условия на гринице вакуума и диэлектрика
# тангенциальная компонента электрического поля
eqEz_a  = E_dielectric_z (a)  == E_vacuum2_z (a)
eqEz_ma = E_dielectric_z (-a) == E_vacuum1_z (-a)

print ("eqEz_a = ",  eqEz_a)
print ("eqEz_ma = ",  eqEz_ma)

# нормальная компонента электрического смещения
eqEx_a  = E_dielectric_x (a) * epsilon  == E_vacuum2_x (a)
eqEx_ma = E_dielectric_x (-a) * epsilon == E_vacuum1_x (-a)

print ("eqEx_a =", eqEx_a)
print ("eqEx_ma =", eqEx_ma)

# eqEz_a =  _B1_z*cos(a*kappa) + _B2_z*sin(a*kappa) == _A1_z*e^(a*s) + _A2_z*e^(-a*s)
# eqEz_ma =  _B1_z*cos(a*kappa) - _B2_z*sin(a*kappa) == _A2_z*e^(a*s) + _A1_z*e^(-a*s)
# eqEx_a = (_B1_x*cos(a*kappa) + _B2_x*sin(a*kappa))*epsilon == _A1_x*e^(a*s) + _A2_x*e^(-a*s)
# eqEx_ma = (_B1_x*cos(a*kappa) - _B2_x*sin(a*kappa))*epsilon == _A2_x*e^(a*s) + _A1_x*e^(-a*s)


# тангенциальная компонента напряженности магнитного поля
eqHy_a  = H_vacuum2_y (a)  == H_dielectric_y (a)
eqHy_ma = H_vacuum1_y (-a) == H_dielectric_y (-a)

print ("eqHy_a =", eqHy_a)
print ("eqHy_ma =", eqHy_ma)

# H_x = c / (I * omega) * (diff(E_z, y) - diff(E_y, z))
# H_z = c / (I * omega) * (diff(E_y, x) - diff(E_x, y))
# E_y = 0

# H_dielectric_x (x) = c / (I * omega) * (diff(E_dielectric_z (x), y))

def GenerateMatrix2(equsys, vars):
    A=matrix([[equ.lhs().coefficient(v) for v in vars] for equ in equsys])
    b=matrix([[equ.rhs()] for equ in equsys])
    return (A,b)

def GenerateMatrix(equsys, vars):
    A=matrix(SR, [[(equ.lhs() - equ.rhs()).coefficient(v) for v in vars] for equ in equsys])
    return A

def GenerateMatrixMult(equsys, vars):
    A=matrix(SR, [[(equ.lhs() - equ.rhs()).coefficient(v) * v for v in vars] for equ in equsys])
    return A

def reduceDependedRows(M):
    print ("")
    print ("M.nrows()", M.nrows())
    print ("M.ncols()", M.ncols())
    print ("M.rank()", M.rank())

    M_pivot_rows = M.pivot_rows()
    print ("M.pivot_rows() =", M_pivot_rows)

    M_rows = M.rows()
    reduced_list = []
    for r in M_pivot_rows:
        print ("M_rows[", r, "] =", M_rows[r])
        reduced_list.append(M_rows[r])

    reduced_M = matrix(SR, len(M_pivot_rows), M.ncols(), reduced_list)
    reduced_M_det = reduced_M.det()
    print ("")
    print("reduced_M.det() =", reduced_M_det)
    return reduced_M

eqsys_boundary_conditions = [eqEx_a,  eqEz_a,  eqHy_a*(-I)*omega/c,
                             eqEx_ma, eqEz_ma, eqHy_ma*(-I)*omega/c,
                            ]


vars = [_B1_x, _B2_x, _B1_z, _B2_z, _A_x, _A_z]

M = GenerateMatrix(eqsys_boundary_conditions, vars)
print ("")
print ("")
print ("M =")
print (M)
print ("")
print (M.det())
print ("")
print ("M.rows()", M.rows())
print ("")
print ("M.nrows()", M.nrows())
print ("M.ncols()", M.ncols())
print ("")
print ("M.rank()", M.rank())
print ("M.det()", M.det())
print ("M.pivot_rows() =", M.pivot_rows())



M = matrix(SR, 6, 6,
[
    [ epsilon*cos(a*kappa),  epsilon*sin(a*kappa),                   0,                   0,     -e^(-a*s),           0 ],
    [                    0,                     0,        cos(a*kappa),        sin(a*kappa),             0,   -e^(-a*s) ],
    [     I*k*cos(a*kappa),      I*k*sin(a*kappa),  kappa*sin(a*kappa), -kappa*cos(a*kappa), -I*k*e^(-a*s), -s*e^(-a*s) ],

    [ epsilon*cos(a*kappa), -epsilon*sin(a*kappa),                   0,                   0,     -e^(-a*s),           0 ],
    [                    0,                     0,        cos(a*kappa),       -sin(a*kappa),             0,   +e^(-a*s) ],
    [     I*k*cos(a*kappa),     -I*k*sin(a*kappa), -kappa*sin(a*kappa), -kappa*cos(a*kappa), -I*k*e^(-a*s), -s*e^(-a*s) ],

])
#                      B1x,                   B2x,                 B1z,                 B2z,            Ax,            Az


print (M)

mdet = M.det()
print ("mdet =", mdet)

print (solve([mdet == 0], s))


eqsys_rotH_dielectric = [eq_rot_H_dielectric_x*(-I)*omega/c,
                         eq_rot_H_dielectric_z*(-I)*omega/c,
                        ]

M = GenerateMatrix(eqsys_rotH_dielectric, vars)
print ("")
print ("")
print ("M =")
print (M)
print ("")
#print (M.det())
print ("")
print ("M.nrows()", M.nrows())
print ("M.ncols()", M.ncols())
print ("M.rank()", M.rank())
print ("M.pivot_rows() =", M.pivot_rows())


#eqsys = [eqEx_a,  eqEz_a,  eqHy_a*(-I)*omega/c,
#         eqEx_ma, eqEz_ma, eqHy_ma*(-I)*omega/c,
#         eq_rot_H_dielectric_x*(-I)*omega/c,
#         eq_rot_H_dielectric_z*(-I)*omega/c,
#        ]

M = GenerateMatrix(eqsys_boundary_conditions + eqsys_rotH_dielectric, vars)
print ("")
print ("")
print ("M =")
print (M)
print ("")
#print (M.det())
print ("")
print ("M.rows()", M.rows())
print ("")
print ("M.nrows()", M.nrows())
print ("M.ncols()", M.ncols())
print ("M.rank()", M.rank())
print ("M.pivot_rows() =", M.pivot_rows())

print ("")
reduced_M = reduceDependedRows(M)
reduced_M_det = reduced_M.det()
print ("")
print (solve([reduced_M_det == 0], s))



# чётные

even_vars_dielectric = [_B1_x, _B2_z]


even_M_rotH_dielectric = GenerateMatrix(eqsys_rotH_dielectric, even_vars_dielectric)
print ("")
print ("")
print ("even_M_rotH_dielectric =")
print (even_M_rotH_dielectric)
print ("")

x = vector(SR, even_vars_dielectric)
y = vector(SR, [0, 0])
eqns = [(even_M_rotH_dielectric*x)[index] == y[index] for index in [0,1]]
print(eqns)

for index in [0,1]:
    eq = eqns[index]
    # print(eq)
    print(solve(eq, even_vars_dielectric[0]))
    # print(solve(eq, even_vars_dielectric[1]))
    print(eq_kappa)

even_vars = [_B1_x, _B2_z, _A_x, _A_z]


even_M = GenerateMatrix(eqsys_boundary_conditions + eqsys_rotH_dielectric, even_vars)
print ("")
print ("")
print ("even_M =")
print (even_M)
print ("")
even_reduced_M = reduceDependedRows(even_M)
even_reduced_M_det = even_reduced_M.det()

print ("")
print("even_reduced_M.det() =", even_reduced_M_det)


M = matrix(SR, 4, 4,
[
    [                                 var("M11"),                     0,      var("M13"),              0 ],
    [                                          0,            var("M22"),               0,     var("M24") ],
    [                             I*k*var("M31"),            var("M32"), -I*k*var("M33"),  -s*var("M34") ],
    [ -k^2*var("a41")+epsilon*omega^2*var("b41"),       -I*k*var("M42"),               0,               0]
])
mdet = M.det()
print ("")
print (mdet)

print ("")
print (solve([even_reduced_M_det == 0], s))
# s == -(c^2*epsilon*k^2*kappa*cos(a*kappa) - epsilon*kappa*omega^2*cos(a*kappa))/(c^2*k^2*sin(a*kappa) - epsilon*omega^2*sin(a*kappa))


# sage: %display ascii_art

#  / 2          2                                         2             \ 
# -\c *epsilon*k *kappa*cos(a*kappa) - epsilon*kappa*omega *cos(a*kappa)/ 
# ------------------------------------------------------------------------
#              2  2                             2                         
#             c *k *sin(a*kappa) - epsilon*omega *sin(a*kappa) 





M = matrix(SR, 4, 4,
[
    [ epsilon*cos(a*kappa),                     0,     -e^(-a*s),           0 ],
    [                    0,          sin(a*kappa),             0,   -e^(-a*s) ],
    [     I*k*cos(a*kappa),   -kappa*cos(a*kappa), -I*k*e^(-a*s), -s*e^(-a*s) ],
    [                kappa,                  -I*k,             0,            0]
])
#                       Bx,                    Bz,            Ax,            Az

mdet = M.det()
print (mdet)
# -epsilon*k^2*cos(a*kappa)*e^(-2*a*s) + k^2*cos(a*kappa)*e^(-2*a*s) + kappa^2*cos(a*kappa)*e^(-2*a*s) + kappa*s*e^(-2*a*s)*sin(a*kappa)

print (solve([mdet == 0], s))

# s == ((epsilon - 1)*k^2 - kappa^2)*cos(a*kappa)/(kappa*sin(a*kappa))

# sage: %display ascii_art

# sage: ((epsilon - 1)*k^2 - kappa^2)*cos(a*kappa)/(kappa*sin(a*kappa))
# / 2                      2\             
# \k *(epsilon - 1) - kappa /*cos(a*kappa)
# ----------------------------------------
#            kappa*sin(a*kappa)    


# нечётные
odd_vars = [_B2_x, _B1_z, _A_x, _A_z]



odd_M = GenerateMatrix(eqsys_boundary_conditions + eqsys_rotH_dielectric, odd_vars)
print ("")
print ("")
print ("odd_M =")
print (odd_M)
print ("")
#print (odd_M.det())
print ("")

odd_reduced_M = reduceDependedRows(odd_M)
odd_reduced_M_det = odd_reduced_M.det()

print ("")
print("odd_reduced_M.det() =", odd_reduced_M_det)

print (solve([odd_reduced_M_det == 0], s))
# s == kappa*sin(a*kappa)/(mu*cos(a*kappa))



'''
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
    _B1_x,

    #kappa, s
    )
print ("res =", res)


'''





