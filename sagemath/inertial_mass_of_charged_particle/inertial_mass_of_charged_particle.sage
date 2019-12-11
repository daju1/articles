import sys
reload(sys)
sys.setdefaultencoding('utf8')

# расмотрим модель заряженной частицы имеющей сферическую симметрию распределения заряда с плотностью rho_q(r)
# задача : найти коэффициент самоиндукции такой частицы L
# 1/(J^2)*integral(integral(vector_j*vector_j/R, dV), dV)
# предположение : величина L * J^2 представляет собой кинетическую энергию частица
# поскольку vector_j = rho_q(r) * vector_v, интеграл
# vector_v*vector_v * integral(rho_q(ra) * integral(rho_q(rq)/R, dVq), dVa)
# по своей форме очень напоминает выражение для кинетической энергии в случае малых скоростей
# m * v^2 / 2 = v * v * integral(integral(rho_q(ra)*rho_q(rq)/R, dVq), dVa)
# в котором роль инертной массы играет
# m = 2 * integral(rho_q(ra) * integral(rho_q(rq)/R, dVq), dVa)
# таким образом отыскав выражение этого интеграла мы можем получить зависимость инертной массы частицы
# от объёмного распределения ее электрического заряда

# для случая больших скоростей в знаменателе нужно использовать радиус Лиенара Вихерта
#

# будем производить вычисления в сферической системе координат, расположенной таким образом,
# что вектор скорости частицы сонаправлен вектору z сферической системы координат

# следуя Тамму, обозначим индексом q координату заряда а индексом a координату точки наблюдения

# радиальная координната заряда
rq = var("rq")
# радиальная координата точки наблюдения
ra = var("ra")


# азимутальная координата заряда
phi_q = var("theta_q")

# азимутальная координата точки наблюдения
phi_a = var("theta_a")

# полярная координата заряда
theta_q = var("theta_q")

# полярная координата точки наблюдения
theta_a = var("theta_a")

# z - координаты заряда и точки наблюдения
zq = rq*cos(theta_q)
za = ra*cos(theta_a)

# физический смысл первого интегрирования по dVq есть отыскание векторного потенциала
# создаваемого движущимся распределённым зарядом частицы в ее объёме

# ввиду сферической симметрии векторный потенциал не зависит от угла phi_a
# следовательно при вычислении векторного потенциала мы можем положить phi_a = 0

# введём вспомогательные переменные - цилиндрический радиус
# как координата r точки при переходе в цилиндрическую систему коорденат с тем же направлением оси z

rca = ra*sin(theta_a)
rcq = rq*sin(theta_q)

# выражение для расстояния между точкой заряда и точкой наблюдения примет вид
R = sqrt(rca^2+rcq^2+(za-zq)^2-2*rca*rcq*cos(phi_q))

# поскольку по условию задачи распределение плотности заряда частицы сферически симметрично,
# мы можем вынести rho_q(rq) за знак интегрирования по phi_q
# интегрируя 1/R по phi_q от 0 до 2*pi запишем

rcqa2 = (rcq-rca)^2+(zq-za)^2
module = - 4*rcq*rca / rcqa2
Iphi=4*elliptic_kc(module) / sqrt(rcqa2)

# http://maths.cnam.fr/Membres/wilk/MathMax/help/Maxima/maxima_17.html
# Function: elliptic_kc (m)
# The complete elliptic integral of the first kind, defined as
# integrate(1/sqrt(1 - m*sin(x)^2), x, 0, %pi/2)

# https://www.math.aau.at/user/cheuberg/sage/doc/6.5.beta2-4-g2e57cf9/en/reference/functions/sage/functions/special.html
# class sage.functions.special.EllipticKC
# Bases: sage.functions.special.MaximaFunction
# This returns the value of the “complete elliptic integral of the first kind,”




print "Iphi =", Iphi
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



I2=integral(Iphi * sin(theta_q), (theta_q, 0, pi), algorithm="giac")
print I2

rho_q(rq) = 1/(rq^2)

I3 = integral(rho_q(rq) * rq^2 * I2, (rq, 0, infinity), algorithm="giac")
print I3

I4 = 2 * pi * I3

I5 = integral(I4 * sin(theta_a), (theta_a, 0, pi), algorithm="giac")
print I5

I6 = integral(rho_q(ra) * ra^2 * I5, (ra, 0, infinity), algorithm="giac")
print I6