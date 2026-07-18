from sage.calculus.calculus import var
from sage.symbolic.assumptions import assume

k_z = var ('k_z')
assume(k_z, 'complex')

kz = var ('kz')
assume(kz, 'real')

sz = var ('sz')
assume(sz, 'real')
assume(sz>0)


omega = var('omega')
assume(omega, 'real')
assume(omega>0)