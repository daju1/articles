import sys
reload(sys)
sys.setdefaultencoding('utf8')

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

    #I2 = lambda ra, theta_a, rq          : my_numerical_integral( lambda theta_q : I1(ra, theta_a, rq, theta_q) * sin(theta_q), 0, pi)
    I2 = lambda ra, theta_a, rq          : my_numerical_integral( lambda theta_q : Iphi(ra, theta_a, rq, theta_q) * sin(theta_q), 0, pi)

    # распределение заряда ядра приближённо выражается распределением Ферми
    # http://nuclphys.sinp.msu.ru/ndb/ndb102.htm

    rho_q = lambda rho0, Rq, aq, r : rho0 / (1 + exp( (r - Rq) / aq) )

    I3 = lambda rho0, Rq, aq, ra, theta_a : my_numerical_integral( lambda rq : rho_q(rho0, Rq, aq, rq) * rq^2 * I2(ra, theta_a, rq),  0, infinity)

    I4 = lambda rho0, Rq, aq, ra, theta_a : 2 * pi * I3(rho0, Rq, aq, ra, theta_a)

    I5 = lambda rho0, Rq, aq, ra : my_numerical_integral( lambda theta_a : I4(rho0, Rq, aq, ra, theta_a) * sin(theta_a), 0, pi)

    I6 = lambda rho0, Rq, aq : my_numerical_integral( lambda ra : rho_q(rho0, Rq, aq, ra) * ra^2 * I5(rho0, Rq, aq, ra), 0, infinity)

    # I6(rho0, Rq, aq)

    I7 = I6(1, 1, 1)

    print "I7 = ", I7
    #m = (mju_0 / (4 * pi)) * I6(rho0, Rq, aq)

def test():
    f = lambda k,xx,yy,zz : k * xx^2 + yy^3 + zz^4;

    I1 = lambda k,xx,yy : my_numerical_integral(lambda zz : f(k,xx,yy,zz), 0, 3 );

    I2 = lambda k,xx   : my_numerical_integral(lambda yy : I1(k,xx,yy), 0, 2 );

    I3 = lambda k     : my_numerical_integral(lambda xx : I2(k,xx), 0, 1 );

    I4 = I3(1.0)

    print "I4 = ", I4

    from sympy import integrate, Symbol
    k = Symbol('k')
    xx = Symbol('xx')
    yy = Symbol('yy')
    zz = Symbol('zz')
    print integrate(f(k,xx,yy,zz), (xx), (yy), (zz))
    print integrate(f(k,xx,yy,zz), (xx, 0, 1), (yy, 0, 2), (zz, 0, 3))
    print integrate(f(k,xx,yy,zz), (xx, 0, 1), (yy, 0, 2), (zz, 0, 3))
    print integrate(f(1,xx,yy,zz), (xx, 0, 1), (yy, 0, 2), (zz, 0, 3)).n()

#calc1_m()
calc2_m()
#test()