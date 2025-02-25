import sys
reload(sys)
sys.setdefaultencoding('utf8')

import numpy as np

# расмотрим модель заряженной частицы имеющей сферическую симметрию распределения заряда с плотностью rho_q(r)
# задача : найти коэффициент самоиндукции такой частицы L
# 1/(J^2)*integral(integral(vector_j*vector_j/R, dV), dV)
# предположение : величина 1/c^2 * L * J^2 / 2 представляет собой кинетическую энергию частицы
# обоснование 1/c^2 * L * J^2 / 2 - энергия магнитного поля, которая равна кинетической энергии (по Николаеву)
# поскольку vector_j = rho_q(r) * vector_v, интеграл
# vector_v*vector_v * integral(rho_q(ra) * integral(rho_q(rq)/R, dVq), dVa)
# по своей форме очень напоминает выражение для кинетической энергии в случае малых скоростей
# m * v^2 / 2 = 1/c^2 * v * v * integral(integral(rho_q(ra)*rho_q(rq)/R, dVq), dVa) / 2
# в котором роль инертной массы играет
# m = 1/c^2 * integral(rho_q(ra) * integral(rho_q(rq)/R, dVq), dVa)
# таким образом отыскав выражение этого интеграла мы можем получить зависимость инертной массы частицы
# от объёмного распределения ее электрического заряда

# для случая больших скоростей в знаменателе нужно использовать радиус Лиенара Вихерта
# R_lw = R - vector_v / c * vector_R

# будем производить вычисления в сферической системе координат, расположенной таким образом,
# что вектор скорости частицы сонаправлен вектору z сферической системы координат

# следуя Тамму, обозначим индексом q координату заряда а индексом a координату точки наблюдения
# vector_R направлен от q к a
# vector_v / c * vector_R = v / c * vector_ort_z * vector_R = v / c * Rz = v / c * (za - zq)

# радиальная координната заряда
# rq = var("rq")
# радиальная координата точки наблюдения
# ra = var("ra")

# assume(rq>0)
# assume(ra>0)

# азимутальная координата заряда
# phi_q = var("theta_q")

# азимутальная координата точки наблюдения
# phi_a = var("theta_a")

# полярная координата заряда
# theta_q = var("theta_q")

# полярная координата точки наблюдения
# theta_a = var("theta_a")

# z - координаты заряда и точки наблюдения
zq = lambda rq, theta_q : rq*cos(theta_q)
za = lambda ra, theta_a : ra*cos(theta_a)

# физический смысл первого интегрирования по dVq есть отыскание векторного потенциала
# создаваемого движущимся распределённым зарядом частицы в ее объёме

# ввиду сферической симметрии векторный потенциал не зависит от угла phi_a
# следовательно при вычислении векторного потенциала мы можем положить phi_a = 0

# введём вспомогательные переменные - цилиндрический радиус
# как координата r точки при переходе в цилиндрическую систему коорденат с тем же направлением оси z

rcq = lambda rq, theta_q : rq*sin(theta_q)
rca = lambda ra, theta_a : ra*sin(theta_a)

# выражение для расстояния между точкой заряда и точкой наблюдения примет вид
# R0 = lambda ra, theta_a, rq, theta_q, phi_q : sqrt(rca^2+rcq^2+(za-zq)^2-2*rca*rcq*cos(phi_q))
R0 = lambda ra, theta_a, rq, theta_q, phi_q : sqrt(rca(ra, theta_a)^2+rcq(rq, theta_q)^2+(za(ra, theta_a)-zq(rq, theta_q))^2-2*rca(ra, theta_a)*rcq(rq, theta_q)*cos(phi_q))

# поскольку по условию задачи распределение плотности заряда частицы сферически симметрично,
# мы можем вынести rho_q(rq) за знак интегрирования по phi_q
# интегрируя 1/R по phi_q от 0 до 2*pi запишем

# rcqa2 = (rcq-rca)^2+(zq-za)^2
# module = - 4*rcq*rca / rcqa2
# Iphi=4*elliptic_kc(module) / sqrt(rcqa2)

rcqa2  = lambda ra, theta_a, rq, theta_q : (rcq(rq,theta_q)-rca(ra,theta_a))^2+(zq(rq, theta_q)-za(ra, theta_a))^2;
module = lambda ra, theta_a, rq, theta_q : - 4*rcq(rq,theta_q)*rca(ra,theta_a) / rcqa2(ra, theta_a, rq, theta_q);
Iphi   = lambda ra, theta_a, rq, theta_q : 4*elliptic_kc(module(ra, theta_a, rq, theta_q)) / sqrt(rcqa2(ra, theta_a, rq, theta_q));


# http://maths.cnam.fr/Membres/wilk/MathMax/help/Maxima/maxima_17.html
# Function: elliptic_kc (m)
# The complete elliptic integral of the first kind, defined as
# integrate(1/sqrt(1 - m*sin(x)^2), x, 0, %pi/2)

# https://www.math.aau.at/user/cheuberg/sage/doc/6.5.beta2-4-g2e57cf9/en/reference/functions/sage/functions/special.html
# class sage.functions.special.EllipticKC
# Bases: sage.functions.special.MaximaFunction
# This returns the value of the “complete elliptic integral of the first kind,”




# print "Iphi =", Iphi(ra, theta_a, rq, theta_q)
# Iphi = 4*elliptic_kc(-4*ra*rq*sin(theta_a)*sin(theta_q)/((ra*cos(theta_a) - rq*cos(theta_q))^2 + (ra*sin(theta_a) - rq*sin(theta_q))^2))/sqrt((ra*cos(theta_a) - rq*cos(theta_q))^2 + (ra*sin(theta_a) - rq*sin(theta_q))^2)

# sage: %display ascii_art
# sage: Iphi
#    /                    -4*ra*rq*sin(theta_a)*sin(theta_q)                     \
# 4*K|---------------------------------------------------------------------------|
#    |                                   2                                      2|
#    \(ra*sin(theta_a) - rq*sin(theta_q))  + (ra*cos(theta_a) - rq*cos(theta_q)) /
# --------------------------------------------------------------------------------
#    _____________________________________________________________________________
#   /                                    2                                      2 
# \/  (ra*sin(theta_a) - rq*sin(theta_q))  + (ra*cos(theta_a) - rq*cos(theta_q))  

# integral(Iphi * sin(theta_q), theta_q)
# integral(Iphi * sin(theta_q), (theta_q, 0, pi))
# RuntimeError: Encountered operator mismatch in maxima-to-sr translation

def my_numerical_integral1(f, a, b):
    print "f = ", f
    print "f(x) = ", f(x)
    print "a = ", a
    print "b = ", b
    integral = numerical_integral(f, a, b)
    print "integral = ", integral
    result = integral[0]
    print "result = ", result
    return result

def my_numerical_integral(f, a, b):
    #print "f = ", f
    #print "f(x) = ", f(x)
    #print "a = ", a
    #print "b = ", b
    from scipy import integrate
    #import ctypes
    try:
        integral = integrate.quad(f, a, b)
        #print "integral = ", integral
        result = integral[0]
        #print "result = ", result
        return result
    except Exception as ex:
        print "ex = ", str(ex)
        print "f = ", f
        print "f(x) = ", f(x)
        print "a = ", a
        print "b = ", b
        integral = numerical_integral(f, a, b)
        print "integral = ", integral
        result = integral[0]
        print "result = ", result
        return result

print exp(-1/x).nintegral(x, 1, 2)

print my_numerical_integral(lambda x : exp(-1/x), 1, 2)

print ( lambda phi_q : 1 / R0 (ra, theta_a, rq, theta_q, phi_q) )


def calc1_m():
    I1 = lambda ra, theta_a, rq, theta_q : ( 1 / R0 (ra, theta_a, rq, theta_q, phi_q) )   .nintegral(phi_q, 0, 2*pi)

    # I2 = lambda ra, theta_a, rq          : ( I1(ra, theta_a, rq, theta_q) * sin(theta_q) ).nintegral(theta_q, 0, pi)
    I2 = lambda ra, theta_a, rq          : ( Iphi(ra, theta_a, rq, theta_q) * sin(theta_q) ).nintegral(theta_q, 0, pi)

    # распределение заряда ядра приближённо выражается распределением Ферми
    # http://nuclphys.sinp.msu.ru/ndb/ndb102.htm

    rho_q = lambda rho0, Rq, aq, r : rho0 / (1 + exp( (r - Rq) / aq) )

    I3 = lambda rho0, Rq, aq, ra, theta_a : ( rho_q(rho0, Rq, aq, rq) * rq^2 * I2(ra, theta_a, rq) ).nintegral(rq, 0, infinity)

    I4 = lambda rho0, Rq, aq, ra, theta_a : 2 * pi * I3(rho0, Rq, aq, ra, theta_a)

    I5 = lambda rho0, Rq, aq, ra : ( I4(rho0, Rq, aq, ra, theta_a) * sin(theta_a)).nintegral(theta_a, 0, pi)

    I6 = lambda rho0, Rq, aq : ( rho_q(rho0, Rq, aq, ra) * ra^2 * I5(rho0, Rq, aq, ra)).nintegral(ra, 0, infinity)

    I6(1, 1, 1)

    #m = (mju_0 / (4 * pi)) * I6(rho0, Rq, aq)

def calc2_m():
    I1 = lambda ra, theta_a, rq, theta_q : my_numerical_integral( lambda phi_q : 1 / R0 (ra, theta_a, rq, theta_q, phi_q), 0, 2*pi)

    I2 = lambda ra, theta_a, rq          : my_numerical_integral( lambda theta_q : I1(ra, theta_a, rq, theta_q) * sin(theta_q), 0, pi)
    # I2 = lambda ra, theta_a, rq          : my_numerical_integral( lambda theta_q : Iphi(ra, theta_a, rq, theta_q) * sin(theta_q), 0, pi)

    # распределение заряда ядра приближённо выражается распределением Ферми
    # http://nuclphys.sinp.msu.ru/ndb/ndb102.htm

    rho_q = lambda rho0, Rq, aq, r : rho0 / (1 + exp( (r - Rq) / aq) )

    I3 = lambda rho0, Rq, aq, ra, theta_a : my_numerical_integral( lambda rq : rho_q(rho0, Rq, aq, rq) * rq^2 * I2(ra, theta_a, rq),  0, infinity)

    I4 = lambda rho0, Rq, aq, ra, theta_a : 2 * pi * I3(rho0, Rq, aq, ra, theta_a)

    I5 = lambda rho0, Rq, aq, ra : my_numerical_integral( lambda theta_a : I4(rho0, Rq, aq, ra, theta_a) * sin(theta_a), 0, pi)

    I6 = lambda rho0, Rq, aq : my_numerical_integral( lambda ra : rho_q(rho0, Rq, aq, ra) * ra^2 * I5(rho0, Rq, aq, ra), 0, infinity)

    # I6(rho0, Rq, aq)

    I7 = I6(1, 1, 1)
    print "I6(1, 1, 1) = ", I7

    # using Iphi
    # integral =  (820.1392952900268, 0.0004965284196251068)
    # result =  820.13929529
    # I6(1, 1, 1) =  820.13929529
    # m = (mju_0 / (4 * pi)) * I6(rho0, Rq, aq)

    # using I1
    # float division by zero
    # ...
    # float division by zero
    # integral =  (17.077164899361996, 0.0721483434946714)
    # result =  17.0771648994
    # ex =  float division by zero
    # f =  <function <lambda> at 0x7fcf9034a7d0>
    # f(x) =  1/sqrt(-10.595854716776243*cos(x) + 10.595854716776243)
    # a =  0
    # b =  2*pi
    # float division by zero
    # ....
    # float division by zero
    # integral =  (17.07716489912273, 0.07214834349503378)
    # result =  17.0771648991
    # math range error
    # math range error
    # math range error
    # Killed

def calc3_scalar_potential():
    # распределение заряда ядра приближённо выражается распределением Ферми
    # http://nuclphys.sinp.msu.ru/ndb/ndb102.htm

    rho_q = lambda rho0, Rq, aq, r : rho0 / (1 + np.exp( (r - Rq) / aq) )

    Ir = lambda rho0, Rq, aq, theta_a, ra, phi_q, theta_q : my_numerical_integral( lambda rq : rho_q(rho0, Rq, aq, rq) * rq^2 * sin(theta_q) / R0 (ra, theta_a, rq, theta_q, phi_q),  0, infinity)

    I2 = lambda rho0, Rq, aq, theta_a, ra, phi_q : my_numerical_integral( lambda theta_q : Ir(rho0, Rq, aq, theta_a, ra, phi_q, theta_q), 0, pi)

    I3 = lambda rho0, Rq, aq, theta_a, ra : my_numerical_integral( lambda phi_q : I2(rho0, Rq, aq, theta_a, ra, phi_q), 0, 2*pi)

    # scalar_potential = I3(1, 1, 1, 0.5, pi/4)
    # print "I3(1, 1, 1, 0.5, pi/4) = ", scalar_potential
    # I3(1, 1, 1, 0.5, pi/4) =  21.8634904207

    I4 = lambda rho0, Rq, aq, theta_a, ra : 2 * pi * I3(rho0, Rq, aq, theta_a, ra)

    I5 = lambda rho0, Rq, aq, theta_a : my_numerical_integral( lambda ra : rho_q(rho0, Rq, aq, ra) * I4(rho0, Rq, aq, theta_a, ra) * sin(theta_a) * ra^2, 0, infinity)

    I6 = lambda rho0, Rq, aq : my_numerical_integral( lambda theta_a : I5 (rho0, Rq, aq, theta_a), 0, pi)

    # I6(rho0, Rq, aq)

    I7 = I6(1, 1, 1)
    print "I6(1, 1, 1) = ", I7
    # sage/src/bin/sage-ipython:220: RuntimeWarning: overflow encountered in exp
    # I6(1, 1, 1) =  820.139143519


def calc_proton_mass():
    # http://www.actaphys.uj.edu.pl/fulltext?series=Reg&vol=30&page=119
    # the rms charge radius of the proton being
    # rp_rms = 0.8
    # r0 = 2/3 * rp_rms
    # the charge distribution of the proton
    # g(r) = exp(-r^2/r__0^2)/(r__0^3*sqrt(pi)^3)

    rho_q = lambda r0, r : exp(-r^2/r0^2)/(r0^3*sqrt(pi)^3)

    Ir = lambda r0, theta_a, ra, phi_q, theta_q : my_numerical_integral( lambda rq : rho_q(r0, rq) * rq^2 * sin(theta_q) / R0 (ra, theta_a, rq, theta_q, phi_q),  0, infinity)

    I2 = lambda r0, theta_a, ra, phi_q : my_numerical_integral( lambda theta_q : Ir(r0, theta_a, ra, phi_q, theta_q), 0, pi)

    I3 = lambda r0, theta_a, ra : my_numerical_integral( lambda phi_q : I2(r0, theta_a, ra, phi_q), 0, 2*pi)

    I4 = lambda r0, theta_a, ra : 2 * pi * I3(r0, theta_a, ra)

    I5 = lambda r0, theta_a : my_numerical_integral( lambda ra : rho_q(r0, ra) * I4(r0, theta_a, ra) * sin(theta_a) * ra^2, 0, infinity)

    I6 = lambda r0 : my_numerical_integral( lambda theta_a : I5 (r0, theta_a), 0, pi)

    I7 = I6(2/3 * 0.8)
    print "I6(2/3 * 0.8) = ", I7
    # Killed

def calc_neutron_mass():
    # http://www.actaphys.uj.edu.pl/fulltext?series=Reg&vol=30&page=119

    # the rms value of g__n
    # r2 := -0.113

    # The third parameter r1 is a scaling parameter, which is necessary to define
    # a dimensionless quantity (r/r1) in the Gaussian exponent. The results of
    # QCD-calculations of the charge density distribution inside the neutron [2]
    # are best reproduced by choosing:
    # r1 = 0.71*sqrt(2/5) # fm

    # the charge density distribution within the neutron
    # gn(r) = (-2/3)*(r2 / (r1^2 * (r1*sqrt(pi))^3)) * (r/r1)^2 * (1 - (2/5)*r^2/r1^2) * exp(-r^2/r1^2)

    rho_q = lambda r1, r2, r : (-2/3)*(r2 / (r1^2 * (r1*sqrt(pi))^3)) * (r/r1)^2 * (1 - (2/5)*r^2/r1^2) * exp(-r^2/r1^2)

    Ir = lambda r1, r2, theta_a, ra, phi_q, theta_q : my_numerical_integral( lambda rq : rho_q(r1, r2, rq) * rq^2 * sin(theta_q) / R0 (ra, theta_a, rq, theta_q, phi_q),  0, infinity)

    I2 = lambda r1, r2, theta_a, ra, phi_q : my_numerical_integral( lambda theta_q : Ir(r1, r2, theta_a, ra, phi_q, theta_q), 0, pi)

    I3 = lambda r1, r2, theta_a, ra : my_numerical_integral( lambda phi_q : I2(r1, r2, theta_a, ra, phi_q), 0, 2*pi)

    I4 = lambda r1, r2, theta_a, ra : 2 * pi * I3(r1, r2, theta_a, ra)

    I5 = lambda r1, r2, theta_a : my_numerical_integral( lambda ra : rho_q(r1, r2, ra) * I4(r1, r2, theta_a, ra) * sin(theta_a) * ra^2, 0, infinity)

    I6 = lambda r1, r2 : my_numerical_integral( lambda theta_a : I5 (r1, r2, theta_a), 0, pi)

    I7 = I6(0.71*sqrt(2/5), -0.113)

    print "I6(0.71*sqrt(2/5), -0.113) = ", I7
    # Killed

nquad_default_opts = \
           { 'epsabs' : 1.49e-08,
             'epsrel' : 1.49e-08,
             'limit'  : 50 }

nquad_opts = { 'epsabs' : 1.49e-03,
               'epsrel' : 1.49e-03,
               'limit'  : 20 }

def calc_proton_mass2():
    # http://www.actaphys.uj.edu.pl/fulltext?series=Reg&vol=30&page=119
    # the rms charge radius of the proton being
    # rp_rms = 0.8
    # r0 = 2/3 * rp_rms
    # the charge distribution of the proton
    # g(r) = exp(-r^2/r__0^2)/(r__0^3*sqrt(pi)^3)

    file = open('calc_proton_mass2.txt', 'a')
    file.write('http://www.actaphys.uj.edu.pl/fulltext?series=Reg&vol=30&page=119\n')
    file.write('the rms charge radius of the proton being\n')
    file.write('rp_rms = 0.8\n')
    file.write('r0 = 2/3 * rp_rms\n')
    file.write('the charge distribution of the proton\n')
    file.write('g(r) = exp(-r^2/r__0^2)/(r__0^3*sqrt(pi)^3)\n')
    file.close()

    rho_q = lambda r0, r : exp(-r^2/r0^2)/(r0^3*sqrt(pi)^3)

    Ir = lambda r0, theta_a, ra, phi_q, theta_q, rq : rho_q(r0, rq) * rq^2 * sin(theta_q) / R0 (ra, theta_a, rq, theta_q, phi_q)    # (rq,  0, infinity)

    I2 = lambda r0, theta_a, ra, phi_q, theta_q, rq : Ir(r0, theta_a, ra, phi_q, theta_q, rq)                                       # (theta_q, 0, pi)

    I3 = lambda r0, theta_a, ra, phi_q, theta_q, rq : I2(r0, theta_a, ra, phi_q, theta_q, rq)                                       # (phi_q, 0, 2*pi)

    I4 = lambda r0, theta_a, ra, phi_q, theta_q, rq : 2 * pi * I3(r0, theta_a, ra, phi_q, theta_q, rq)

    I5 = lambda r0, theta_a, ra, phi_q, theta_q, rq : rho_q(r0, ra) * I4(r0, theta_a, ra, phi_q, theta_q, rq) * sin(theta_a) * ra^2 # (ra, 0, infinity)

    I6 = lambda r0, theta_a, ra, phi_q, theta_q, rq : I5 (r0, theta_a, ra, phi_q, theta_q, rq)                                      # (theta_a, 0, pi)

    I7 = lambda theta_a, ra, phi_q, theta_q, rq : I6(2/3 * 0.8, theta_a, ra, phi_q, theta_q, rq)

    from scipy import integrate
    #                              theta_a, ra,            phi_q,      theta_q,  rq
    answer = integrate.nquad(I7, [ [0, pi], [0, infinity], [0, 2*pi],  [0, pi],  [0, infinity]], opts=nquad_opts)

    print "I6(2/3 * 0.8) = ", answer
    # I6(2/3 * 0.8) =  (1.4960348943817992, 0.0026474827067254846)

    # Elementary charge, coulombs
    e = 1.60217662e-19
    mju0 = 4*pi*10^(-7) # H/m
    k = mju0 / (4 * pi) * e^2 / 10^-15
    m = k * answer[0]
    print "m = ", m
    # 3.84027657565375e-30

    me = 9.1093837015e-31
    print "m / me = ", m / me
    # 4.21573698231790

    file = open('calc_proton_mass2.txt', 'a')
    file.write('\n')
    file.write('I6(2/3 * 0.8) = ')
    file.write(str(answer))
    file.write('\n')
    file.write('m = ')
    file.write(str(m))
    file.write('\n\n')
    file.close()

def calc_neutron_mass2():
    # http://www.actaphys.uj.edu.pl/fulltext?series=Reg&vol=30&page=119

    # the rms value of g__n
    # r2 := -0.113

    # The third parameter r1 is a scaling parameter, which is necessary to define
    # a dimensionless quantity (r/r1) in the Gaussian exponent. The results of
    # QCD-calculations of the charge density distribution inside the neutron [2]
    # are best reproduced by choosing:
    # r1 = 0.71*sqrt(2/5) # fm

    # the charge density distribution within the neutron
    # gn(r) = (-2/3)*(r2 / (r1^2 * (r1*sqrt(pi))^3)) * (r/r1)^2 * (1 - (2/5)*r^2/r1^2) * exp(-r^2/r1^2)

    file = open('calc_neutron_mass2.txt', 'a')
    file.write('http://www.actaphys.uj.edu.pl/fulltext?series=Reg&vol=30&page=119\n')
    file.write('the rms value of g__n\n')
    file.write('r2 := -0.113\n')
    file.write('\n')
    file.write('The third parameter r1 is a scaling parameter, which is necessary to define\n')
    file.write('a dimensionless quantity (r/r1) in the Gaussian exponent. The results of\n')
    file.write('QCD-calculations of the charge density distribution inside the neutron [2]\n')
    file.write('are best reproduced by choosing:\n')
    file.write('r1 = 0.71*sqrt(2/5) # fm\n')
    file.write('\n')
    file.write('the charge density distribution within the neutron\n')
    file.write('gn(r) = (-2/3)*(r2 / (r1^2 * (r1*sqrt(pi))^3)) * (r/r1)^2 * (1 - (2/5)*r^2/r1^2) * exp(-r^2/r1^2)\n')
    file.close()

    rho_q = lambda r1, r2, r : (-2/3)*(r2 / (r1^2 * (r1*sqrt(pi))^3)) * (r/r1)^2 * (1 - (2/5)*r^2/r1^2) * exp(-r^2/r1^2)

    Ir = lambda r1, r2, theta_a, ra, phi_q, theta_q, rq : rho_q(r1, r2, rq) * rq^2 * sin(theta_q) / R0 (ra, theta_a, rq, theta_q, phi_q)              # (rq, 0, infinity)

    I2 = lambda r1, r2, theta_a, ra, phi_q, theta_q, rq : Ir(r1, r2, theta_a, ra, phi_q, theta_q, rq)                                                     # (theta_q, 0, pi)

    I3 = lambda r1, r2, theta_a, ra, phi_q, theta_q, rq : I2(r1, r2, theta_a, ra, phi_q, theta_q, rq)                                                              # (phi_q, 0, 2*pi)

    I4 = lambda r1, r2, theta_a, ra, phi_q, theta_q, rq : 2 * pi * I3(r1, r2, theta_a, ra, phi_q, theta_q, rq)

    I5 = lambda r1, r2, theta_a, ra, phi_q, theta_q, rq : rho_q(r1, r2, ra) * I4(r1, r2, theta_a, ra, phi_q, theta_q, rq) * sin(theta_a) * ra^2       # (ra, 0, infinity)

    I6 = lambda r1, r2, theta_a, ra, phi_q, theta_q, rq : I5 (r1, r2, theta_a, ra, phi_q, theta_q, rq )                                               # (theta_a, 0, pi)

    I7 = lambda theta_a, ra, phi_q, theta_q, rq  : I6(0.71*sqrt(2/5), -0.113, theta_a, ra, phi_q, theta_q, rq)

    from scipy import integrate
    #                              theta_a, ra,            phi_q,      theta_q,  rq
    answer = integrate.nquad(I7, [ [0, pi], [0, infinity], [0, 2*pi],  [0, pi],  [0, infinity]], opts=nquad_opts)

    print "I6(0.71*sqrt(2/5), -0.113) = ", answer

    file = open('calc_neutron_mass2.txt', 'a')
    file.write('\n')
    file.write('I6(0.71*sqrt(2/5), -0.113) = ')
    file.write(str(answer))
    file.write('\n\n')
    file.close()
    # Killed

def calc_sphere_mass():

    file = open('calc_sphere_mass.txt', 'a')
    file.close()

    rho_q = lambda r0, q : 3*q / (4*pi*r0^3)

    Ir = lambda r0, q, theta_a, ra, phi_q, theta_q, rq : rho_q(r0, q) * rq^2 * sin(theta_q) / R0 (ra, theta_a, rq, theta_q, phi_q)    # (rq,  0, infinity)

    I2 = lambda r0, q, theta_a, ra, phi_q, theta_q, rq : Ir(r0, q, theta_a, ra, phi_q, theta_q, rq)                                       # (theta_q, 0, pi)

    I3 = lambda r0, q, theta_a, ra, phi_q, theta_q, rq : I2(r0, q, theta_a, ra, phi_q, theta_q, rq)                                       # (phi_q, 0, 2*pi)

    I4 = lambda r0, q, theta_a, ra, phi_q, theta_q, rq : 2 * pi * I3(r0, q, theta_a, ra, phi_q, theta_q, rq)

    I5 = lambda r0, q, theta_a, ra, phi_q, theta_q, rq : rho_q(r0, q) * I4(r0, q, theta_a, ra, phi_q, theta_q, rq) * sin(theta_a) * ra^2 # (ra, 0, infinity)

    I6 = lambda r0, q, theta_a, ra, phi_q, theta_q, rq : I5 (r0, q, theta_a, ra, phi_q, theta_q, rq)                                      # (theta_a, 0, pi)

    Iq1 = lambda r0, theta_a, ra, phi_q, theta_q, rq : I6(r0, 1, theta_a, ra, phi_q, theta_q, rq)

    from scipy import integrate
    #                                                                                                                      theta_a, ra,      phi_q,      theta_q,  rq
    Ir0 = lambda r0 : integrate.nquad(lambda theta_a, ra, phi_q, theta_q, rq : Iq1(r0, theta_a, ra, phi_q, theta_q, rq), [ [0, pi], [0, r0], [0, 2*pi],  [0, pi],  [0, r0]], opts=nquad_opts)
    # answer = Ir0(0.1)
    # print "Ir0(0.1) = ", answer
    # Killed

    answer = Ir0(1.0)
    print "Ir0(1) = ", answer
    # Ir0(1) =  (1.1999998963704677, 0.005411660477958158)

    file = open('calc_sphere_mass.txt', 'a')
    file.write('\n')
    file.write('Ir0(1) = ')
    file.write(str(answer))
    file.close()

def legendre_summ(l, theta_q, phi_q, theta_a, phi_a):
    from sage.functions.special import spherical_harmonic
    f = spherical_harmonic(l, m, theta_q, phi_q) * conjugate(spherical_harmonic(l, m, theta_a, phi_a))
    return 4 * pi / (2*l + 1) * sum(f, m, -l, l)

def legendre_summ_of_mass_of_proton(l):
    # http://www.actaphys.uj.edu.pl/fulltext?series=Reg&vol=30&page=119
    # the rms charge radius of the proton being
    # rp_rms = 0.8
    # r0 = 2/3 * rp_rms
    # the charge distribution of the proton
    # g(r) = exp(-r^2/r__0^2)/(r__0^3*sqrt(pi)^3)

    theta_q, phi_q = var('theta_q, phi_q')
    assume(theta_q, 'real')
    assume(phi_q, 'real')

    theta_a, phi_a = var('theta_a, phi_a')
    assume(theta_a, 'real')
    assume(phi_a, 'real')

    r_q, r_a = var("r_q, r_a")
    R = var("R")

    rho_q = lambda r0, r : exp(-r^2/r0^2)/(r0^3*sqrt(pi)^3)

    print latex(rho_q(var("r0"), r_q))

    assume(r_a>0)
    assume(r_a<R)
    assume(R>0)

    R = infinity
    rp_rms = 0.8
    r0 = 2/3 * rp_rms

    # if r_q < r_a
    potential_1 = ((1/r_a)*((r_q/r_a)^l)*rho_q(r0, r_q)*sin(theta_q)*(r_q^2)*legendre_summ(l, theta_q, phi_q, theta_a, phi_a)).integrate(theta_q, 0, pi).integrate(phi_q, 0, 2*pi).integrate(r_q, 0, r_a)
    # if r_a < r_q
    potential_2 = ((1/r_q)*((r_a/r_q)^l)*rho_q(r0, r_q)*sin(theta_q)*(r_q^2)*legendre_summ(l, theta_q, phi_q, theta_a, phi_a)).integrate(theta_q, 0, pi).integrate(phi_q, 0, 2*pi).integrate(r_q, r_a, R)

    potential = potential_1 + potential_2

    return (potential * rho_q(r0, r_a) * sin(theta_a) * r_a^2).integrate(r_a, 0, R).integrate(theta_a, 0, pi).integrate(phi_a, 0, 2*pi)

def legendre_summ_of_mass_of_neutron(l):
    # http://www.actaphys.uj.edu.pl/fulltext?series=Reg&vol=30&page=119

    # the rms value of g__n
    r2 = -0.113

    # The third parameter r1 is a scaling parameter, which is necessary to define
    # a dimensionless quantity (r/r1) in the Gaussian exponent. The results of
    # QCD-calculations of the charge density distribution inside the neutron [2]
    # are best reproduced by choosing:
    r1 = 0.71*sqrt(2/5) # fm

    # the charge density distribution within the neutron
    # gn(r) = (-2/3)*(r2 / (r1^2 * (r1*sqrt(pi))^3)) * (r/r1)^2 * (1 - (2/5)*r^2/r1^2) * exp(-r^2/r1^2)


    theta_q, phi_q = var('theta_q, phi_q')
    assume(theta_q, 'real')
    assume(phi_q, 'real')

    theta_a, phi_a = var('theta_a, phi_a')
    assume(theta_a, 'real')
    assume(phi_a, 'real')

    r_q, r_a = var("r_q, r_a")
    R = var("R")

    rho_q = lambda r1, r2, r : (-2/3)*(r2 / (r1^2 * (r1*sqrt(pi))^3)) * (r/r1)^2 * (1 - (2/5)*r^2/r1^2) * exp(-r^2/r1^2)

    print latex(rho_q(var("r1"), var("r2"), r_q))

    assume(r_a>0)
    assume(r_a<R)
    assume(R>0)

    R = infinity

    # if r_q < r_a
    potential_1 = ((1/r_a)*((r_q/r_a)^l)*rho_q(r1, r2, r_q)*sin(theta_q)*(r_q^2)*legendre_summ(l, theta_q, phi_q, theta_a, phi_a)).integrate(theta_q, 0, pi).integrate(phi_q, 0, 2*pi).integrate(r_q, 0, r_a)
    # if r_a < r_q
    potential_2 = ((1/r_q)*((r_a/r_q)^l)*rho_q(r1, r2, r_q)*sin(theta_q)*(r_q^2)*legendre_summ(l, theta_q, phi_q, theta_a, phi_a)).integrate(theta_q, 0, pi).integrate(phi_q, 0, 2*pi).integrate(r_q, r_a, R)

    potential = potential_1 + potential_2

    return (potential * rho_q(r1, r2, r_a) * sin(theta_a) * r_a^2).integrate(r_a, 0, R).integrate(theta_a, 0, pi).integrate(phi_a, 0, 2*pi)

def legendre_summ_of_inductivity_of_sphere(l):
    theta_q, phi_q = var('theta_q, phi_q')
    assume(theta_q, 'real')
    assume(phi_q, 'real')

    theta_a, phi_a = var('theta_a, phi_a')
    assume(theta_a, 'real')
    assume(phi_a, 'real')

    r_q, r_a = var("r_q, r_a")
    R = var("R")
    rho = var("rho")

    assume(r_a>0)
    assume(r_a<R)
    assume(R>0)

    # if r_q < r_a
    A1 = ((1/r_a)*((r_q/r_a)^l)*rho*v*sin(theta_q)*(r_q^2)*legendre_summ(l, theta_q, phi_q, theta_a, phi_a)).integrate(theta_q, 0, pi).integrate(phi_q, 0, 2*pi).integrate(r_q, 0, r_a)
    # if r_a < r_q
    A2 = ((1/r_q)*((r_a/r_q)^l)*rho*v*sin(theta_q)*(r_q^2)*legendre_summ(l, theta_q, phi_q, theta_a, phi_a)).integrate(theta_q, 0, pi).integrate(phi_q, 0, 2*pi).integrate(r_q, r_a, R)

    A = A1 + A2

    V = 4/3*pi*R^3
    q = rho*V
    # J = (rho*v).integrate(dS) -> (mean) -> (rho*v).integrate(dV) / (2*R) = q*v / (2 * R)
    J = q*v / (2 * R)

    return var("mju0")/(4*pi) * (A * rho * v * sin(theta_a) * r_a^2).integrate(r_a, 0, R).integrate(theta_a, 0, pi).integrate(phi_a, 0, 2*pi) / (J^2)

def legendre_summ_of_mass_of_sphere(l):
    theta_q, phi_q = var('theta_q, phi_q')
    assume(theta_q, 'real')
    assume(phi_q, 'real')

    theta_a, phi_a = var('theta_a, phi_a')
    assume(theta_a, 'real')
    assume(phi_a, 'real')

    r_q, r_a = var("r_q, r_a")
    R = var("R")
    rho = var("rho")

    assume(r_a>0)
    assume(r_a<R)
    assume(R>0)

    # if r_q < r_a
    potential_1 = ((1/r_a)*((r_q/r_a)^l)*rho*sin(theta_q)*(r_q^2)*legendre_summ(l, theta_q, phi_q, theta_a, phi_a)).integrate(theta_q, 0, pi).integrate(phi_q, 0, 2*pi).integrate(r_q, 0, r_a)
    # if r_a < r_q
    potential_2 = ((1/r_q)*((r_a/r_q)^l)*rho*sin(theta_q)*(r_q^2)*legendre_summ(l, theta_q, phi_q, theta_a, phi_a)).integrate(theta_q, 0, pi).integrate(phi_q, 0, 2*pi).integrate(r_q, r_a, R)

    potential = potential_1 + potential_2

    return (potential * rho * sin(theta_a) * r_a^2).integrate(r_a, 0, R).integrate(theta_a, 0, pi).integrate(phi_a, 0, 2*pi)

def legendre_summ_of_inductivity_of_spherical_shell(l):
    theta_q, phi_q = var('theta_q, phi_q')
    assume(theta_q, 'real')
    assume(phi_q, 'real')

    theta_a, phi_a = var('theta_a, phi_a')
    assume(theta_a, 'real')
    assume(phi_a, 'real')

    R = var("R")
    # rho = var("q")/(4*pi*R^2)*dirac_delta(r_a-R)
    sigma = var("q")/(4*pi*R^2)

    assume(R>0)

    A = (1/R)*(sigma*v*sin(theta_q)*(R^2)*legendre_summ(l, theta_q, phi_q, theta_a, phi_a)).integrate(theta_q, 0, pi).integrate(phi_q, 0, 2*pi)

    # S = 4*pi*R^2
    # q = rho*V

    # J = (sigma*v).integrate(dl) -> (mean) -> (sigma*v).integrate(dS) / (2*R)
    J = q * v / (2 * R)
    J = (sigma*v*sin(theta_q)*(R^2)).integrate(theta_q, 0, pi).integrate(phi_q, 0, 2*pi) / (2 * R)
    print "J =", J

    return var("mju0")/(4*pi) * (A * sigma * v * sin(theta_a) * R^2).integrate(theta_a, 0, pi).integrate(phi_a, 0, 2*pi) / (J^2)

def legendre_summ_of_mass_of_spherical_shell(l):
    theta_q, phi_q = var('theta_q, phi_q')
    assume(theta_q, 'real')
    assume(phi_q, 'real')

    theta_a, phi_a = var('theta_a, phi_a')
    assume(theta_a, 'real')
    assume(phi_a, 'real')

    R = var("R")
    # rho = var("q")/(4*pi*R^2)*dirac_delta(r_a-R)
    sigma = var("q")/(4*pi*R^2)

    assume(R>0)

    potential = ((1/R)*sigma*sin(theta_q)*(R^2)*legendre_summ(l, theta_q, phi_q, theta_a, phi_a)).integrate(theta_q, 0, pi).integrate(phi_q, 0, 2*pi)

    return (potential * sigma * sin(theta_a) * R^2).integrate(theta_a, 0, pi).integrate(phi_a, 0, 2*pi)

def make_spherical_plot3d(f):
    # Color plots on surface of sphere
    var('x,y')
    cm = colormaps.spring
    cf = lambda x,y: (sin(x) + cos(y)) % 1
    return spherical_plot3d(f, (phi_a,0,2*pi), (theta_a,0,pi), color=(cf,cm))


def calc_inductivity_of_sphere():
    from sage.functions.special import spherical_harmonic

    theta_q, phi_q = var('theta_q, phi_q')
    assume(theta_q, 'real')
    assume(phi_q, 'real')

    theta_a, phi_a = var('theta_a, phi_a')
    assume(theta_a, 'real')
    assume(phi_a, 'real')

    l,m = var('l,m')
    f = spherical_harmonic(l, m, theta_q, phi_q) * conjugate(spherical_harmonic(l, m, theta_a, phi_a))
    print sum(f, m, -l, l)
    print legendre_summ(l, theta_q, phi_q, theta_a, phi_a)

    for l in range(0, 4):
        print "l =", str(l), " legendre_summ = ", legendre_summ(l, theta_q, phi_q, theta_a, phi_a).simplify()
        print "l =", str(l), " legendre_summ.integral = ", (sin(theta_q)*legendre_summ(l, theta_q, phi_q, theta_a, phi_a)).integrate(theta_q, 0, pi)
        print "l =", str(l), " legendre_summ.integral = ", (sin(theta_a)*(sin(theta_q)*legendre_summ(l, theta_q, phi_q, theta_a, phi_a)).integrate(theta_q, 0, pi)).integrate(theta_a, 0, pi)
        print "l =", str(l), " legendre_summ.integral = ", (sin(theta_a)*(sin(theta_q)*legendre_summ(l, theta_q, phi_q, theta_a, phi_a)).integrate(theta_q, 0, pi)).integrate(theta_a, 0, pi).integrate(phi_q, 0, 2*pi)
        print "l =", str(l), " legendre_summ.integral = ", (sin(theta_a)*(sin(theta_q)*legendre_summ(l, theta_q, phi_q, theta_a, phi_a)).integrate(theta_q, 0, pi)).integrate(theta_a, 0, pi).integrate(phi_q, 0, 2*pi).integrate(phi_a, 0, 2*pi)
        print "legendre_summ_of_inductivity_of_sphere(", l, ") =", legendre_summ_of_inductivity_of_sphere(l)
        L = legendre_summ_of_inductivity_of_spherical_shell(l)
        print "legendre_summ_of_inductivity_of_spherical_shell(", l, ") =", L

    # legendre_summ_of_inductivity_of_sphere( 0 ) = 6/5*R*mju0/pi
    # legendre_summ_of_inductivity_of_spherical_shell( 0 ) = R*mju0/pi

def calc_mass_of_sphere():
    print "\ncalc_mass_of_sphere"
    # https://en.wikipedia.org/wiki/Electromagnetic_mass
    m = legendre_summ_of_mass_of_sphere(0)
    print "m =", m
    # m = 32/15*pi^2*R^5*rho^2

    V = 4/3*pi*R^3
    m = m.subs(rho = q/V)
    print "m =", m
    # m = 6/5*q^2/R

    m = var("mju0")/(4*pi) * m
    print "m =", m
    # m = 3/10*mju0*q^2/(pi*R)

    # Classical electron radius
    # https://en.wikipedia.org/wiki/
    # Классический радиус электрона равен радиусу полой сферы, на которой равномерно распределён заряд, если этот заряд равен заряду электрона, а потенциальная энергия электростатического поля {\displaystyle U_{0}\ } U_{0}\ полностью эквивалентна половине массы электрона (без учета квантовых эффектов):
    # {\displaystyle U_{0}={\frac {1}{2}}{\frac {1}{4\pi \varepsilon _{0}}}\cdot {\frac {e^{2}}{r_{0}}}={\frac {1}{2}}m_{0}c^{2}}.
    r_e = 1/(4*pi*var("epsilon_0"))*q^2/(var("m_e")*var("c")^2)
    print "r_e =", r_e

    # собственная электрическая энергия заряженного шара радиуса r_0 равна 3/5 * e^2/r_0, если заряд равномерно распределён по всему объёму шара (Тамм, стр.88)
    # для такой модели электрона классический радиус будет иметь уже другое значение
    # U = 3/5 * 1/(4*pi*epsilon_0) * e^2/r_0 = 1/2*m*c^2

    r_e = 3/5 * 1/(4*pi*var("epsilon_0"))*q^2/((1/2)*var("m_e")*var("c")^2)
    print "r_e =", r_e

    m = m.subs(R = r_e)
    print "m =", m
    # m = 6/5*c^2*epsilon_0*m_e*mju0
    # для классического радиуса перерасчитанного в модели заряда распределённого по всему объёму шара
    # m = c^2*epsilon_0*m_e*mju0

    # https://en.wikipedia.org/wiki/Vacuum_permittivity
    m = m.subs(epsilon_0 = 1/(mju0*c^2))
    print "m =", m
    # m = 6/5*m_e
    # для классического радиуса перерасчитанного в модели заряда распределённого по всему объёму шара
    # m = m_e

def calc_mass_of_spherical_shell():
    print "\ncalc_mass_of_spherical_shell"
    # https://en.wikipedia.org/wiki/Electromagnetic_mass
    m = legendre_summ_of_mass_of_spherical_shell(0)
    print "m =", m
    # m = q^2/R

    m = var("mju0")/(4*pi) * m
    print "m =", m
    # m = 1/4*mju0*q^2/(pi*R)

    # Classical electron radius
    # https://en.wikipedia.org/wiki/
    # Классический радиус электрона равен радиусу полой сферы, на которой равномерно распределён заряд, если этот заряд равен заряду электрона, а потенциальная энергия электростатического поля {\displaystyle U_{0}\ } U_{0}\ полностью эквивалентна половине массы электрона (без учета квантовых эффектов):
    # {\displaystyle U_{0}={\frac {1}{2}}{\frac {1}{4\pi \varepsilon _{0}}}\cdot {\frac {e^{2}}{r_{0}}}={\frac {1}{2}}m_{0}c^{2}}.
    # U = 1/2 * 1/(4*pi*epsilon_0) * e^2/r_0 = 1/2*m*c^2
    # собственная электрическая энергия заряженного шара радиуса r_0 равна 1/2 * e^2/r_0, если заряд распределён по поверхности шара (Тамм, стр.88)

    r_e = 1/(4*pi*var("epsilon_0"))*q^2/(var("m_e")*var("c")^2)
    print "r_e =", r_e

    m = m.subs(R = r_e)
    print "m =", m
    # m = c^2*epsilon_0*m_e*mju0

    # https://en.wikipedia.org/wiki/Vacuum_permittivity
    m = m.subs(epsilon_0 = 1/(mju0*c^2))
    print "m =", m
    # m = m_e

def calc_mass_of_proton():
    print "\ncalc_mass_of_proton"
    m0 = legendre_summ_of_mass_of_proton(0)
    print "legendre_summ_of_mass_of_proton(0) =", m0, m0.n()

    m1 = legendre_summ_of_mass_of_proton(1)
    print "legendre_summ_of_mass_of_proton(1) =", m1

    m2 = legendre_summ_of_mass_of_proton(2)
    print "legendre_summ_of_mass_of_proton(2) =", m2

    print latex(m0)

    # legendre_summ_of_mass_of_proton(0) = 1.875*sqrt(2)/sqrt(pi) 1.49603355150537
    # legendre_summ_of_mass_of_proton(1) = 0
    # legendre_summ_of_mass_of_proton(2) = 0
    # VEGAS RESULT:   1.49695439 +- 0.00345851        p = 0.006
    # SUAVE RESULT:   1.49424975 +- 0.00149158        p = 1.000
    # scipy integrate.quad calc_proton_mass results (for comparing):
    # I6(2/3 * 0.8) =  (1.4960348943817992, 0.0026474827067254846)

    # Elementary charge, coulombs
    e = 1.60217662e-19
    mju0 = 4*pi*10^(-7) # H/m
    k = mju0 / (4 * pi) * e^2 / 10^-15
    m = k * m0
    print "m = ", m, m.n()
    print latex(m)
    # 3.84027312853036e-30
    # 3.84027657565375e-30 - scipy integrate.nquad result


    me = 9.1093837015e-31
    print "m / me = ", m / me, (m / me).n()
    # 4.21573319817234
    # 4.21573698231790 - scipy integrate.nquad result

def calc_mass_of_neutron():
    print "\ncalc_mass_of_neutron"
    m0 = legendre_summ_of_mass_of_neutron(0)
    print "legendre_summ_of_mass_of_neutron(0) =", m0
    print "legendre_summ_of_mass_of_neutron(0) =", m0.n()

    m1 = legendre_summ_of_mass_of_neutron(1)
    print "legendre_summ_of_mass_of_neutron(1) =", m1

    m2 = legendre_summ_of_mass_of_neutron(2)
    print "legendre_summ_of_mass_of_neutron(2) =", m2

    print latex(m0)

    # legendre_summ_of_mass_of_neutron(0) = (4.48918680252563e-28)*sqrt(5)*sqrt(2)*(1824320471*sqrt(10)*sqrt(5) + 14369256122481640385175552*sqrt(2))/sqrt(pi)
    # legendre_summ_of_mass_of_neutron(0) = 0.0162757880193542
    # legendre_summ_of_mass_of_neutron(1) = 0
    # legendre_summ_of_mass_of_neutron(2) = 0
    # VEGAS RESULT:   0.00183803 +- 0.00000179        p = 0.000
    # SUAVE RESULT:   0.00183606 +- 0.00000183        p = 1.000
    # Elementary charge, coulombs
    e = 1.60217662e-19
    mju0 = 4*pi*10^(-7) # H/m
    k = mju0 / (4 * pi) * e^2 / 10^-15
    m = k * m0
    print "m = ", m
    print "m = ", m.n()
    # m =  (1.1523607494861974e-57)*sqrt(5)*sqrt(2)*(1824320471*sqrt(10)*sqrt(5) + 14369256122481640385175552*sqrt(2))/sqrt(pi) 4.17794582972344e-32
    # m / me =  (1.2650260294738088e-27)*sqrt(5)*sqrt(2)*(1824320471*sqrt(10)*sqrt(5) + 14369256122481640385175552*sqrt(2))/sqrt(pi) 0.0458641985739987

    print latex(m)

    me = 9.1093837015e-31
    print "m / me = ", m / me, (m / me).n()



def legendre_summ_of_vector_potencial_of_rotated_sphere(l):
    theta_q, phi_q = var('theta_q, phi_q')
    assume(theta_q, 'real')
    assume(phi_q, 'real')

    theta_a, phi_a = var('theta_a, phi_a')
    assume(theta_a, 'real')
    assume(phi_a, 'real')
    phi_a = 0

    r_q, r_a = var("r_q, r_a")
    R = var("R")
    v = var("v")
    rho = var("rho")

    assume(r_a>0)
    assume(r_a<R)
    assume(R>0)

    # if r_q < r_a
    A1 = ((1/r_a)*((r_q/r_a)^l)*rho*v*cos(phi_q)*sin(theta_q)*(r_q^2)*legendre_summ(l, theta_q, phi_q, theta_a, phi_a)).integrate(theta_q, 0, pi).integrate(phi_q, 0, 2*pi).integrate(r_q, 0, r_a)
    # if r_a < r_q
    A2 = ((1/r_q)*((r_a/r_q)^l)*rho*v*cos(phi_q)*sin(theta_q)*(r_q^2)*legendre_summ(l, theta_q, phi_q, theta_a, phi_a)).integrate(theta_q, 0, pi).integrate(phi_q, 0, 2*pi).integrate(r_q, r_a, R)

    return A1 + A2

def legendre_summ_of_vector_potencial_of_rotated_solid_sphere(l):
    theta_q, phi_q = var('theta_q, phi_q')
    assume(theta_q, 'real')
    assume(phi_q, 'real')

    theta_a, phi_a = var('theta_a, phi_a')
    assume(theta_a, 'real')
    assume(phi_a, 'real')

    phi_a = 0

    r_q = var("r_q")
    omega = var("omega")
    R = var("R")
    v = omega*r_q*sin(theta_q)
    rho = var("rho")

    assume(r_a>0)
    assume(r_a<R)
    assume(R>0)

    # if r_q < r_a
    A1 = ((1/r_a)*((r_q/r_a)^l)*rho*v*cos(phi_q)*sin(theta_q)*(r_q^2)*legendre_summ(l, theta_q, phi_q, theta_a, phi_a)).integrate(theta_q, 0, pi).integrate(phi_q, 0, 2*pi).integrate(r_q, 0, r_a)
    # if r_a < r_q
    A2 = ((1/r_q)*((r_a/r_q)^l)*rho*v*cos(phi_q)*sin(theta_q)*(r_q^2)*legendre_summ(l, theta_q, phi_q, theta_a, phi_a)).integrate(theta_q, 0, pi).integrate(phi_q, 0, 2*pi).integrate(r_q, r_a, R)

    return A1 + A2


def magnetic_moment(v, rho, R):
    return 2*pi * (r_a^2 * (sin(theta_a))^3 * v * rho).integrate(r_a, 0, R).integrate(theta_a, 0, pi)
def angular_momentum(A, rho, R):
    return 2*pi * (r_a^2 * (sin(theta_a))^3 * A * rho).integrate(r_a, 0, R).integrate(theta_a, 0, pi)
def orbital_angular_momentum(v, rho_m, R):
    return 2*pi * (r_a^2 * (sin(theta_a))^3 * v * rho_m).integrate(r_a, 0, R).integrate(theta_a, 0, pi)

def calc_gyromagnetic_ratio_of_sphere():
    # the ratio of its magnetic moment to its angular momentum
    #for l in range(0, 4):
    #    print "l =", str(l), " legendre_summ = ", legendre_summ(l, theta_q, phi_q, theta_a, phi_a).simplify()
    #    print "l =", str(l), " legendre_summ.integral theta_q = ", (cos(phi_q)*sin(theta_q)*legendre_summ(l, theta_q, phi_q, theta_a, phi_a)).integrate(theta_q, 0, pi)
    #    print "l =", str(l), " legendre_summ.integral phi_q = ", (cos(phi_q)*sin(theta_q)*legendre_summ(l, theta_q, phi_q, theta_a, phi_a)).integrate(theta_q, 0, pi).integrate(phi_q, 0, 2*pi)

    l_max = 5

    A = legendre_summ_of_vector_potencial_of_rotated_sphere(0)
    for l in range(1, l_max+1):
        print "l = ", l
        dA = legendre_summ_of_vector_potencial_of_rotated_sphere(l)
        print "dA = ", dA
        A += dA
    print "A = ", A
    # make_spherical_plot3d(A.real()).show(aspect_ratio=(1,1,1))
    A = A.subs(phi_a=0)
    print "A = ", A

    rho = var("rho")

    mass = legendre_summ_of_mass_of_sphere(0)
    print "mass =", mass

    rho_m = mass/(4/3*pi*R^3)
    print "rho_m =", rho_m

    M_orbital = orbital_angular_momentum(v, rho_m, R)
    print "M_orbital =", M_orbital

    M = angular_momentum(A, rho, R)
    print "M =", M

    m = magnetic_moment(v, rho, R)
    print "m =", m

    g = m/M
    print "g =", g

    g_orbital = m/M_orbital
    print "g_orbital =", g_orbital

    # plot3d(A.subs(v = 1, rho = 1, R = 1), (r_a, 0, 1), (theta_a,0,pi)).show(aspect_ratio=(1,1,1))
    p = plot(A.subs(v = 1, rho = 1, R = 1, theta_a=pi/2), (r_a, 0, 1))
    p.save("A.png")


    A_solid = legendre_summ_of_vector_potencial_of_rotated_solid_sphere(0)
    for l in range(1, l_max+1):
        print "l = ", l
        dA_solid = legendre_summ_of_vector_potencial_of_rotated_solid_sphere(l)
        print "dA_solid = ", dA_solid
        A_solid += dA_solid
    print "A_solid = ", A_solid
    # make_spherical_plot3d(A_solid.real()).show(aspect_ratio=(1,1,1))
    A_solid = A_solid.subs(phi_a=0)
    print "A_solid = ", A_solid

    v_solid = omega*r_a*sin(theta_a)
    Mo_solid = orbital_angular_momentum(v_solid, rho_m, R)
    print "Mo_solid =", Mo_solid

    M_orbital_solid = orbital_angular_momentum(v_solid, rho_m, R)
    print "M_orbital_solid =", M_orbital_solid

    M_solid = angular_momentum(A_solid, rho, R)
    print "M_solid =", M_solid

    m_solid = magnetic_moment(v_solid, rho, R)
    print "m_solid =", m_solid

    g_solid = m_solid/M_solid
    print "g_solid =", g_solid

    g_orbital_solid = m_solid/M_orbital_solid
    print "g_orbital_solid =", g_orbital_solid

    print "g =", g
    print "g_solid =", (g_solid*rho*R^2).n()/(rho*R^2)
    print "g =", (g*rho*R^2).n()/(rho*R^2)
    print "g_orbital_solid =", (g_orbital_solid*rho*R^2).n()/(rho*R^2)
    print "g_orbital =", (g_orbital*rho*R^2).n()/(rho*R^2)

    rho_uniform = var("q")/(4/3*pi*R^3)
    g = g.subs(rho = rho_uniform)
    g_solid = g_solid.subs(rho = rho_uniform)
    g_orbital = g_orbital.subs(rho = rho_uniform)
    g_orbital_solid = g_orbital_solid.subs(rho = rho_uniform)

    print "g/g_orbital =", g/g_orbital
    print "g_solid/g_orbital_solid =", g_solid/g_orbital_solid

    print "g_solid =", g_solid
    print "g =", g
    print "g_orbital_solid =", g_orbital_solid
    print "g_orbital =", g_orbital

    print "g_solid =", (g_solid*q/R).n()/(q/R)
    print "g =", (g*q/R).n()/(q/R)
    print "g_orbital_solid =", (g_orbital_solid*q/R).n()/(q/R)
    print "g_orbital =", (g_orbital*q/R).n()/(q/R)

    print "g/g_orbital =", (g/g_orbital).n()
    print "g_solid/g_orbital_solid =", (g_solid/g_orbital_solid).n()

    # plot3d(A_solid.subs(omega = 1, rho = 1, R = 1), (r_a, 0, 1), (theta_a,0,pi)).show(aspect_ratio=(1,1,1))
    p = plot(A_solid.subs(omega = 1, rho = 1, R = 1, theta_a=pi/2), (r_a, 0, 1))
    p.save("A_solid.png")

    # p = make_spherical_plot3d(spherical_harmonic(1, 0, theta_a, phi_a).real())
    #p = make_spherical_plot3d(legendre_summ_of_vector_potencial_of_rotated_sphere(0).real())
    #p += make_spherical_plot3d(legendre_summ_of_vector_potencial_of_rotated_sphere(1).real())
    #p += make_spherical_plot3d(legendre_summ_of_vector_potencial_of_rotated_sphere(2).real())
    #p += make_spherical_plot3d(legendre_summ_of_vector_potencial_of_rotated_sphere(3).real())
    #p += make_spherical_plot3d(legendre_summ_of_vector_potencial_of_rotated_sphere(4).real())
    #p += make_spherical_plot3d(legendre_summ_of_vector_potencial_of_rotated_sphere(5).real())
    #p.show(aspect_ratio=(1,1,1))

Sq = function('Sq') # x^2
Cb = function('Cb') # x^3
Qu = function('Qu') # x^4
Fi = function('Fi') # x^5
Si = function('Si') # x^6
Se = function('Si') # x^7
Ei = function('Ei') # x^8
Ni = function('Ni') # x^9

def replace_powers(res):
    return res.subs({
        R_0^2 : Sq(R_0),
        dz^2  : Sq(dz),
        c^2   : Sq(c),
        v^2   : Sq(v),
        a^2   : Sq(a),

        1/(a^2)  : (1/Sq(a)),

        R_0^3 : Cb(R_0),
        dz^3  : Cb(dz),
        c^3   : Cb(c),
        v^3   : Cb(v),
        a^3   : Cb(a),

        1/(a^3)  : (1/Cb(a)),

        R_0^4 : Qu(R_0),
        dz^4  : Qu(dz),
        c^4   : Qu(c),
        v^4   : Qu(v),
        a^4   : Qu(a),

        1/(a^4)  : (1/Qu(a)),

        R_0^5 : Fi(R_0),
        dz^5  : Fi(dz),
        c^5   : Fi(c),
        v^5   : Fi(v),
        a^5   : Fi(a),

        1/(a^5)  : (1/Fi(a)),

        R_0^6 : Si(R_0),
        dz^6  : Si(dz),
        c^6   : Si(c),
        v^6   : Si(v),
        a^6   : Si(a),

        1/(a^6)  : (1/Si(a)),

        R_0^7 : Se(R_0),
        dz^7  : Se(dz),
        c^7   : Se(c),
        v^7   : Se(v),
        a^7   : Se(a),

        1/(a^7)  : (1/Se(a)),

        R_0^8 : Ei(R_0),
        dz^8  : Ei(dz),
        c^8   : Ei(c),
        v^8   : Ei(v),
        a^8   : Ei(a),

        1/(a^8)  : (1/Ei(a)),

        R_0^9 : Ni(R_0),
        dz^9  : Ni(dz),
        c^9   : Ni(c),
        v^9   : Ni(v),
        a^9   : Ni(a),
      })

def replase_integer_frac(r):
    s = str(r)
    import re
    digits_re = r'\d+\/'
    s_double_frac = re.sub(r'([1-9]+)\/([1-9]+)', r'\1.0/\2.0', s)
    return s_double_frac

def solve_lagging():
    # Чтобы учесть запаздывание следует решить систему уравнений
    # s = v*(t-t') + a/2*(t-t')^2
    # и
    # R = c*(t-t')
    # Учитывая, что по теореме косинусов
    # R^2 = R_0^2 + s^2 - 2*R_0*s*cos(alpha)
    #     = R_0^2 + s^2 - 2*R_0*s*(z_q-z_a) / R_0
    #     = R_0^2 + s^2 + 2*s*(z_a-z_q)
    # где R_0 расстояние от точки источника заряда к точке наблюдения без учёта запаздывания.
    # R^2 = c^2*( t-t')^2 = R_0^2 + s^2 + 2*s*(z_a-z_q)

    # (t-t') = dt
    # (z_a-z_q) = dz


    dt = var("dt")
    dz = var("dz")
    a = var("a")
    v = var("v")
    c = var("c")
    z_q, z_a = var("z_q, z_a")
    R_0 = var("R_0")
    s = var("s")

    assume(dz <= R_0)

    # eq = c^2*dt^2 - (R_0^2 + s^2 + 2*s*(z_a-z_q))
    eq = c^2*dt^2 - (R_0^2 + s^2 + 2*s*dz)
    print "eq = ", eq
    eq2 = eq.subs(s = (v*dt + a/2*dt^2))
    print "eq2 = ", eq2

    res2 = solve(eq2, dt)
    # print "res2 = ", res2

    eq3 = eq2.subs(v = 0)
    # print "eq3 = ", eq3
    print "eq3 = ", latex(eq3)

    res3 = solve(eq3, dt)
    # print "res3 = ", res3

    for i in range (0, len(res3)):
        # print "\n"
        r = simplify(c*res3[i])
        # print r
        # print latex(r)
        r3 = replace_powers(r)
        s3 = replase_integer_frac(r3)
        print "res3 =",  s3

    eq4 = eq2.subs(a = 0)
    print "eq4 = ", eq4
    print "eq4 = ", latex(eq4)


    print latex((z_a-z_q) - dz)

    res4 = solve(eq4, dt)
    print "res4 = ", res4
    print "res4 = ", latex(res4)

    for i in range (0, len(res4)):
        r4 = replace_powers(c*res4[i])
        s4 = replase_integer_frac(r4)
        print s4

    for i in range (0, len(res2)):
        # print "\n"
        # print c*res2[i]
        print "\n"
        r2 = replace_powers(c*res2[i])
        s2 = replase_integer_frac(r2)
        # print s2


#calc1_m()
#calc2_m()
#test()
#calc3_scalar_potential()

# calc_proton_mass2()
# calc_neutron_mass2()
# calc_sphere_mass()

#calc_inductivity_of_sphere()
#calc_gyromagnetic_ratio_of_sphere()
#calc_mass_of_sphere()
#calc_mass_of_spherical_shell()

#calc_mass_of_proton()
#calc_mass_of_neutron()

solve_lagging()
