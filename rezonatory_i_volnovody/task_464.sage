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

eq_kappa_s = eq_kappa.subs(solve(eq_s, k^2)).full_simplify()
print("eq_kappa_s =", eq_kappa_s)

eq_s_kappa = eq_s.subs(solve(eq_kappa, k^2)).full_simplify()
print("eq_s_kappa =", eq_s_kappa)

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
_A1_z = var ('_A1_z')
_A2_z = var ('_A2_z')

_B1_x = var ('_B1_x')
_B2_x = var ('_B2_x')
_A1_x = var ('_A1_x')
_A2_x = var ('_A2_x')

# E_dielectric_z = E_dielectric.subs(var('_K1'), _B1_z).subs(var('_K2'), _B2_z)

# выражения для электрического поля в диэлектрике и в вакууме
E_dielectric_z (x) = _B1_z*cos(kappa*x) + _B2_z*sin(kappa*x)
E_dielectric_x (x) = _B1_x*cos(kappa*x) + _B2_x*sin(kappa*x)

# x > a
E_vacuum2_z (x) = _A2_z*e^(-s*x)
E_vacuum2_x (x) = _A2_x*e^(-s*x)

# x < -a
E_vacuum1_z (x) = _A1_z*e^(s*x)
E_vacuum1_x (x) = _A1_x*e^(s*x)

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

def GenerateMatrixSubs(equsys, vars, vars_subs):
    A=matrix(SR, [[(equ.lhs() - equ.rhs()).subs(vars_subs).coefficient(v) for v in vars] for equ in equsys])
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
    for row in M_rows:
        print ("row =", row)

    reduced_list = []
    for r in M_pivot_rows:
        print ("M_rows[", r, "] =", M_rows[r])
        reduced_list.append(M_rows[r])

    reduced_M = matrix(SR, len(M_pivot_rows), M.ncols(), reduced_list)
    reduced_M_det = reduced_M.det()
    print ("")
    print("reduced_M.det() =", reduced_M_det)
    return reduced_M

eqsys_boundary_conditions = [eqEx_a,  eqEz_a,  # eqHy_a*(-I)*omega/c,
                             eqEx_ma, eqEz_ma, # eqHy_ma*(-I)*omega/c,
                            ]


vars = [_B1_x, _B2_x, _B1_z, _B2_z, _A1_x, _A2_x, _A1_z, _A2_z]

M = GenerateMatrix(eqsys_boundary_conditions, vars)
print ("")
print ("")
print ("M =")
print (M)
print ("")
# print (M.det())
print ("")
print ("M.rows()", M.rows())
print ("")
print ("M.nrows()", M.nrows())
print ("M.ncols()", M.ncols())
print ("")
print ("M.rank()", M.rank())
# print ("M.det()", M.det())
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

eqsys_rotH_vacuum = [ # eq_rot_H_vacuum1_x*(-I)*omega/c,
                      eq_rot_H_vacuum2_x*(-I)*omega/c,
                      # eq_rot_H_vacuum1_z*(-I)*omega/c,
                      eq_rot_H_vacuum2_z*(-I)*omega/c,
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

'''
M = GenerateMatrix(eqsys_boundary_conditions + eqsys_rotH_dielectric + eqsys_rotH_vacuum, vars)
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
# reduced_M_det = reduced_M.det()
print ("")
# print (solve([reduced_M_det == 0], s))
'''

def rot_H_solve(M_rotH, vars, eq_s_or_kappa, s_or_kappa):
    x = vector(SR, vars)
    y = vector(SR, [0, 0])
    eqns = [(M_rotH*x)[index] == y[index] for index in [0,1]]
    print(eqns)
    res = []

    for index in [0,1]:
        eq = eqns[index]
        # print(eq)
        sol = solve(eq, vars[0])
        print("sol=", sol[0])
        # print(solve(eq, vars[1]))
        print("eq_s_or_kappa =", eq_s_or_kappa)
        sol_s_or_kappa_2 = solve(eq_s_or_kappa, s_or_kappa^2)
        print("sol_s_or_kappa_2[0].rhs()=", sol_s_or_kappa_2[0].rhs())
        print("")
        s = sol[0].rhs().collect(vars[1])
        ans1 = s.subs(sol_s_or_kappa_2[0].rhs()*c^2==sol_s_or_kappa_2[0].lhs()*c^2)
        ans2 = s.subs(-sol_s_or_kappa_2[0].rhs()*c^2==-sol_s_or_kappa_2[0].lhs()*c^2)

        ans = ans1
        if ans1.number_of_operands() > ans2.number_of_operands():
            ans = ans2
        print(ans)
        if omega in ans.arguments():
            print(ans.arguments())
            sol_omega_2 = solve(eq_s_or_kappa, omega^2)
            print("sol_omega_2 =", sol_omega_2[0].rhs())
            som = ans.subs(omega^2 == sol_omega_2[0].rhs())
            print("som =", som)
            som = som.full_simplify()
            print("som =", som)
            res.append(vars[0] == som)
        else:
            res.append(vars[0] == ans)

        print("")
    return res


# чётные
even_vars = [_B1_x, _B2_z, _A2_x, _A2_z]
eqsys_even = [_A1_z == - _A2_z, _A1_x == _A2_x]

even_vars_dielectric = [_B1_x, _B2_z]
even_eqs_dielectric = [_B2_x == 0, _B1_z == 0]

even_M_rotH_dielectric = GenerateMatrix(eqsys_rotH_dielectric, even_vars_dielectric)
print ("")
print ("even_M_rotH_dielectric =")
print (even_M_rotH_dielectric)
print ("")

print ("")
even_res_rotH_dielectric = rot_H_solve(even_M_rotH_dielectric, even_vars_dielectric, eq_kappa, kappa)
print ("even_res_rotH_dielectric =", even_res_rotH_dielectric)
# even_res_rotH_dielectric = [_B1_x == I*_B2_z*k/kappa, _B1_x == I*_B2_z*k/kappa]

even_vars_vacuum = [_A2_x, _A2_z]

even_M_rotH_vacuum = GenerateMatrixSubs(eqsys_rotH_vacuum, even_vars_vacuum,  eqsys_even)
print ("")
print ("even_M_rotH_vacuum =")
print (even_M_rotH_vacuum)
print ("")

even_res_rotH_vacuum = rot_H_solve(even_M_rotH_vacuum, even_vars_vacuum, eq_s, s)
print ("even_res_rotH_vacuum =", even_res_rotH_vacuum)
# even_res_rotH_vacuum = [_A2_x == I*_A2_z*k/s, _A2_x == I*_A2_z*k/s]



even_M = GenerateMatrixSubs(eqsys_boundary_conditions + even_res_rotH_dielectric + even_res_rotH_vacuum, even_vars,  eqsys_even)
print ("")
print ("")
print ("even_M =")
print (even_M)
print ("")
even_reduced_M = reduceDependedRows(even_M)
even_reduced_M_det = even_reduced_M.det()


print ("")
even_disp_eq = solve([even_reduced_M_det == 0], s)
print ("even_disp_eq =", even_disp_eq)

# s == kappa*sin(a*kappa)/(epsilon*cos(a*kappa))

even_H_dielectric_y (x) = H_dielectric_y (x).subs(even_res_rotH_dielectric[0]).subs(even_eqs_dielectric).subs(solve(eq_kappa,k^2)).full_simplify()
even_H_vacuum1_y    (x) = H_vacuum1_y    (x).subs(eqsys_even).subs(even_res_rotH_vacuum[0]).subs(solve(eq_s,k^2)).full_simplify()
even_H_vacuum2_y    (x) = H_vacuum2_y    (x).subs(eqsys_even).subs(even_res_rotH_vacuum[0]).subs(solve(eq_s,k^2)).full_simplify()

print ("")
print ("even_H_dielectric_y (x) ", even_H_dielectric_y (x))
print ("even_H_vacuum1_y (x) =", even_H_vacuum1_y (x))
print ("even_H_vacuum2_y (x) =", even_H_vacuum2_y (x))
print ("")


even_E_dielectric_z (x) = E_dielectric_z (x).subs(even_res_rotH_dielectric[0]).subs(even_eqs_dielectric).subs(solve(eq_kappa,k^2)).full_simplify()
even_E_dielectric_x (x) = E_dielectric_x (x).subs(even_res_rotH_dielectric[0]).subs(even_eqs_dielectric).subs(solve(eq_kappa,k^2)).full_simplify()

# x > a
even_E_vacuum2_z (x) = E_vacuum2_z (x).subs(eqsys_even).subs(even_res_rotH_vacuum[0]).subs(solve(eq_s,k^2)).full_simplify()
even_E_vacuum2_x (x) = E_vacuum2_x (x).subs(eqsys_even).subs(even_res_rotH_vacuum[0]).subs(solve(eq_s,k^2)).full_simplify()

# x < -a
even_E_vacuum1_z (x) = E_vacuum1_z (x).subs(eqsys_even).subs(even_res_rotH_vacuum[0]).subs(solve(eq_s,k^2)).full_simplify()
even_E_vacuum1_x (x) = E_vacuum1_x (x).subs(eqsys_even).subs(even_res_rotH_vacuum[0]).subs(solve(eq_s,k^2)).full_simplify()


print ("even_E_dielectric_z (x) =", even_E_dielectric_z (x))
print ("even_E_dielectric_x (x) =", even_E_dielectric_x (x))
print ("")

print ("even_E_vacuum2_z (x) =", even_E_vacuum2_z (x))
print ("even_E_vacuum2_x (x) =", even_E_vacuum2_x (x))
print ("")

print ("even_E_vacuum1_z (x) =", even_E_vacuum1_z (x))
print ("even_E_vacuum1_x (x) =", even_E_vacuum1_x (x))
print ("")

# тангенциальная компонента напряженности магнитного поля
even_eqHy_a  = even_H_vacuum2_y (a)  == even_H_dielectric_y (a)
even_eqHy_ma = even_H_vacuum1_y (-a) == even_H_dielectric_y (-a)


even_A = solve (even_eqHy_a, _A2_z)
print ("even_A =", even_A)
even_A = ((even_A[0].rhs()/e^(a*s)).subs(even_disp_eq).full_simplify())*e^(a*s)
print ("even_A =", even_A)

even_eq_kappa = eq_kappa_s.subs(even_disp_eq).full_simplify()
print("even_eq_kappa =", even_eq_kappa)

even_sol_kappa = solve(even_eq_kappa, kappa)
print("even_sol_kappa =", even_sol_kappa)

# sgs
A = 0.1
kappa_A_max = 13*pi.n()/2
s_A_max = 20
digit_values = [omega == 2*pi.n()*10^11, mu == 1, epsilon == 5, c == 299792458 * 100]

even_sol_kappa_d = even_sol_kappa[0].subs(digit_values)
print("even_sol_kappa_d =", even_sol_kappa_d)

even_disp_eq_d = even_disp_eq[0].subs(digit_values)
print ("even_disp_eq_d =", even_disp_eq_d)

even_disp_eq_da = even_disp_eq_d * a
print ("even_disp_eq_da =", even_disp_eq_da)

print (even_disp_eq_d.rhs().subs(a * kappa == x).subs(kappa == x))

r = omega / c * sqrt(epsilon * mu - 1)
print ("r =", r)

r = r.subs(digit_values)
print ("r =", r)

even_disp_eq_fa = lambda x : even_disp_eq_d.rhs().subs(kappa * a == x).subs(kappa == x)
p = plot([even_disp_eq_fa, sqrt((r * A)^2 - x^2)], (x, 0, kappa_A_max), ymin = 0, ymax = s_A_max)
p.save("even_disp_eq_fa.png")


even_disp_eq_f = lambda x : even_disp_eq_d.rhs().subs(kappa == x).subs(a == A)
p = plot([even_disp_eq_f, sqrt((r)^2 - x^2)], (x, 0, kappa_A_max / A), ymin = 0, ymax = s_A_max / A)
p.save("even_disp_eq_f.png")


even_f = lambda x : even_disp_eq_f(x) - sqrt((r)^2 - x^2)
print ("even_f(x) =",even_f(x))
p = plot(even_f, (0, kappa_A_max / A), ymin = - s_A_max / A, ymax = s_A_max / A)
p.save("even_f.png")

even_kappa_sol = find_root(even_f(x) == 0, 0, pi/2/A)
print("even_kappa_sol", even_kappa_sol)

even_s_sol = even_disp_eq_d.rhs()
print("even_s_sol =", even_s_sol)

even_s_sol = even_s_sol.subs(kappa == even_kappa_sol)
print("even_s_sol =", even_s_sol)

even_s_sol = even_s_sol.subs(a == A)
print("even_s_sol =", even_s_sol)

even_k_sol = solve(eq_kappa, k)
print("even_k_sol =", even_k_sol)

even_k_sol = abs(even_k_sol[0].rhs()).subs(digit_values).subs(kappa == even_kappa_sol)
print("even_k_sol =", even_k_sol)

p = plot(even_E_dielectric_z (x).subs(_B2_z == 1, kappa == even_kappa_sol), (x, 0, A))
p.save("even_E_dielectric_z.png")
p = plot(even_E_dielectric_x (x).subs(_B2_z == 1, kappa == even_kappa_sol, k = even_k_sol).imag(), (x, 0, A))
p.save("even_E_dielectric_x.png")
p = plot(even_H_dielectric_y (x).subs(_B2_z == 1, kappa == even_kappa_sol).subs(digit_values).imag(), (x, 0, A))
p.save("even_H_dielectric_y.png")

# Abraham force
# (epsilon * mu - 1) / (4 * pi * c) * ([E x H]).diff(t)
# (epsilon * mu - 1) / (4 * pi * c) * ([E.diff(t) x H] + [E x H.diff(t)])
# ([a x b])_x = a_y * b_z - a_z * b_y
# f_x = (epsilon * mu - 1) / (4 * pi * c) * (-E_z.diff(t) * H_y - E_z * H_y.diff(t))
# E_z.diff(t) = -I * omega * E_z
# H_y.diff(t) = -I * omega * H_y
# f_x = (epsilon * mu - 1) / (4 * pi * c) * (E_z * H_y + E_z * H_y) * I * omega
# f_x = (epsilon * mu - 1) / (2 * pi * c) * (E_z * H_y) * I * omega

f_x(x) = (epsilon * mu - 1) / (2 * pi * c) * (even_E_dielectric_z (x) * even_H_dielectric_y (x)) * I * omega
print("f_x =", f_x(x))
# f_x = -1/2*(epsilon*mu - 1)*_B2_z^2*epsilon*omega^2*cos(kappa*x)*sin(kappa*x)/(pi*c^2*kappa)

f_x(x) = f_x(x).subs(digit_values)
print("f_x =", f_x(x))

f_x(x) = f_x(x).subs(_B2_z == 1, kappa == even_kappa_sol)
print("f_x =", f_x(x))
p = plot(f_x(x), (x, -A, A))
p.save("f_x.png")

F_x = integrate(f_x(x), (x, 0, A))
print("F_x =", F_x)
# F_x = -9.774789731929044/pi


# superconductor force
f_superconductor = - (even_H_dielectric_y (0)^2) / (8* pi)
print("f_superconductor =", f_superconductor)
# f_superconductor = 1/8*_B2_z^2*epsilon^2*omega^2/(pi*c^2*kappa^2)

f_superconductor = f_superconductor.subs(digit_values).subs(_B2_z == 1, kappa == even_kappa_sol)
print("f_superconductor =", f_superconductor)
# f_superconductor = 6.14483801317679/pi


# vacuum field impuls
# 1 / (4 * pi * c) * ([E x H])
# p_x = 1 / (4 * pi * c) * (- E_z * H_y)

p_x(x) = - 1 / (4 * pi * c) * even_E_vacuum2_z (x) * even_H_vacuum2_y (x)
print("p_x =", p_x(x))

p_x(x) = p_x(x).subs( _A2_z == even_A)
print("p_x =", p_x(x))

p_x(x) = p_x(x).subs( s == even_s_sol).subs(digit_values).subs(_B2_z == 1, kappa == even_kappa_sol, a == A)
print("p_x =", p_x(x))
P_x = integrate(p_x(x), (x, A, infinity))
print("P_x", P_x)
# P_x -(5.665029595810918e-14*I)/pi


# нечётные
odd_vars = [_B2_x, _B1_z, _A2_x, _A2_z]
eqsys_odd = [_A1_z == _A2_z, _A1_x == -_A2_x]

odd_vars_dielectric = [_B2_x, _B1_z]
odd_eqs_dielectric = [_B1_x == 0, _B2_z == 0]


odd_M_rotH_dielectric = GenerateMatrix(eqsys_rotH_dielectric, odd_vars_dielectric)
print ("")
print ("odd_M_rotH_dielectric =")
print (odd_M_rotH_dielectric)
print ("")

odd_res_rotH_dielectric = rot_H_solve(odd_M_rotH_dielectric, odd_vars_dielectric, eq_kappa, kappa)
print ("odd_res_rotH_dielectric =", odd_res_rotH_dielectric)
# odd_res_rotH_dielectric = [_B2_x == -I*_B1_z*k/kappa, _B2_x == -I*_B1_z*k/kappa]


odd_vars_vacuum = [_A2_x, _A2_z]

odd_M_rotH_vacuum = GenerateMatrixSubs(eqsys_rotH_vacuum, odd_vars_vacuum,  eqsys_odd)
print ("")
print ("odd_M_rotH_vacuum =")
print (odd_M_rotH_vacuum)
print ("")

odd_res_rotH_vacuum = rot_H_solve(odd_M_rotH_vacuum, odd_vars_vacuum, eq_s, s)
print ("odd_res_rotH_vacuum =", odd_res_rotH_vacuum)
# odd_res_rotH_vacuum = [_A2_x == I*_A2_z*k/s, _A2_x == I*_A2_z*k/s]



odd_M = GenerateMatrixSubs(eqsys_boundary_conditions + odd_res_rotH_dielectric + odd_res_rotH_vacuum, odd_vars, eqsys_odd)
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

odd_disp_eq = solve([odd_reduced_M_det == 0], s)
print ("odd_disp_eq =", odd_disp_eq)

# s == -kappa*cos(a*kappa)/(epsilon*sin(a*kappa))

odd_H_dielectric_y (x) = H_dielectric_y (x).subs(odd_res_rotH_dielectric[0]).subs(odd_eqs_dielectric).subs(solve(eq_kappa,k^2)).full_simplify()
odd_H_vacuum1_y    (x) = H_vacuum1_y    (x).subs(eqsys_odd).subs(odd_res_rotH_vacuum[0]).subs(solve(eq_s,k^2)).full_simplify()
odd_H_vacuum2_y    (x) = H_vacuum2_y    (x).subs(eqsys_odd).subs(odd_res_rotH_vacuum[0]).subs(solve(eq_s,k^2)).full_simplify()

print ("")
print ("odd_H_dielectric_y (x) ", odd_H_dielectric_y (x))
print ("odd_H_vacuum1_y (x) =", odd_H_vacuum1_y (x))
print ("odd_H_vacuum2_y (x) =", odd_H_vacuum2_y (x))
print ("")


odd_E_dielectric_z (x) = E_dielectric_z (x).subs(odd_res_rotH_dielectric[0]).subs(odd_eqs_dielectric).subs(solve(eq_kappa,k^2)).full_simplify()
odd_E_dielectric_x (x) = E_dielectric_x (x).subs(odd_res_rotH_dielectric[0]).subs(odd_eqs_dielectric).subs(solve(eq_kappa,k^2)).full_simplify()

# x > a
odd_E_vacuum2_z (x) = E_vacuum2_z (x).subs(eqsys_odd).subs(odd_res_rotH_vacuum[0]).subs(solve(eq_s,k^2)).full_simplify()
odd_E_vacuum2_x (x) = E_vacuum2_x (x).subs(eqsys_odd).subs(odd_res_rotH_vacuum[0]).subs(solve(eq_s,k^2)).full_simplify()

# x < -a
odd_E_vacuum1_z (x) = E_vacuum1_z (x).subs(eqsys_odd).subs(odd_res_rotH_vacuum[0]).subs(solve(eq_s,k^2)).full_simplify()
odd_E_vacuum1_x (x) = E_vacuum1_x (x).subs(eqsys_odd).subs(odd_res_rotH_vacuum[0]).subs(solve(eq_s,k^2)).full_simplify()


print ("odd_E_dielectric_z (x) =", odd_E_dielectric_z (x))
print ("odd_E_dielectric_x (x) =", odd_E_dielectric_x (x))
print ("")

print ("odd_E_vacuum2_z (x) =", odd_E_vacuum2_z (x))
print ("odd_E_vacuum2_x (x) =", odd_E_vacuum2_x (x))
print ("")

print ("odd_E_vacuum1_z (x) =", odd_E_vacuum1_z (x))
print ("odd_E_vacuum1_x (x) =", odd_E_vacuum1_x (x))
print ("")

# тангенциальная компонента напряженности магнитного поля
odd_eqHy_a  = odd_H_vacuum2_y (a)  == odd_H_dielectric_y (a)
odd_eqHy_ma = odd_H_vacuum1_y (-a) == odd_H_dielectric_y (-a)


odd_A = solve (odd_eqHy_a, _A2_z)
print ("odd_A =", odd_A)
odd_A = ((odd_A[0].rhs()/e^(a*s)).subs(odd_disp_eq).full_simplify())*e^(a*s)
print ("odd_A =", odd_A)




