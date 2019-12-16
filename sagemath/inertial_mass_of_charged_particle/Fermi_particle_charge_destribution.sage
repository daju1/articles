import sys
reload(sys)
sys.setdefaultencoding('utf8')

# распределение заряда ядра приближённо выражается распределением Ферми
# http://nuclphys.sinp.msu.ru/ndb/ndb102.htm

rho_q = lambda rho0, Rq, aq, r : rho0 / (1 + exp( (r - Rq) / aq) )

Rq = var("Rq")
aq = var("aq")
Rmax = var("Rmax")
assume (Rmax > 0);
assume (Rq > 0);
assume (aq > 0);
# float not a valid assumption, must be one of [integer, noninteger, even, odd, rational, irrational, real, imaginary, complex, analytic, increasing, decreasing, oddfun, evenfun, posfun, constant, commutative, lassociative, rassociative, symmetric, antisymmetric, integervalued, one_to_one]
assume(Rq, 'noninteger');
assume(aq, 'noninteger');

integral = integrate(rho_q(1, Rq, aq, x), (x, 0, Rmax))
print "integral = ", integral

# integral =  aq*(Rmax/aq - log((e^(Rmax/aq) + e^(Rq/aq))*e^(-Rq/aq)) + log((e^(Rq/aq) + 1)*e^(-Rq/aq)))
integral_q = lambda Rq, aq, Rmax : aq*(Rmax/aq - log((e^(Rmax/aq) + e^(Rq/aq))*e^(-Rq/aq)) + log((e^(Rq/aq) + 1)*e^(-Rq/aq)))

derivRmax = integral_q(Rq, aq, Rmax).derivative(Rmax)
print "derivRmax = ", derivRmax
# derivRmax =  aq*(1/aq - e^(Rmax/aq)/(aq*(e^(Rmax/aq) + e^(Rq/aq))))

derivRmax_q = lambda Rq, aq, Rmax : aq*(1/aq - e^(Rmax/aq)/(aq*(e^(Rmax/aq) + e^(Rq/aq))))

fsolve (derivRmax_q(Rq, aq, x) - 0.001)
#Rq = 10
#aq = 0.5

#p = plot(derivRmax_q(Rq, aq, x), (x, 0, 100))

#p = plot(integral_q(Rq, aq, x), (x, 0, 100))

#p = plot(rho_q(1, Rq, aq, x), (x, 0, 100))
#p.show()