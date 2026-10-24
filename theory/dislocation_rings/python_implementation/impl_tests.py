# -*- coding: utf-8 -*-
# impl_tests.py
# Воспроизведение всех численных результатов сессии из независимая Python-реализацияов.
# Запуск: python3 impl_tests.py
import numpy as np
from impl_v4_field import (make_electron, make_mobius_particle, force, elem_charge,
                           FERMI_PLAIN, FERMI_FILM_RETARDED, FERMI_FILM_INSTANT)
from impl_v4_stability import ring_equilibrium_radius, analytic_a_eq
from impl_mobius_analysis import scaling_exponent, dens_mobius, dens_dipole, kink_density

def test_dipole_dipole():
    print("=== Диполь-диполь: F*R^2 (три режима) ===")
    c=1.0
    for R in [5,10,20,40,80,160]:
        A=make_electron((0,0,0)); B=make_electron((R,0,0))
        vals=[np.linalg.norm(force(A,B,0.0,c,m))*R*R
              for m in (FERMI_PLAIN,FERMI_FILM_RETARDED,FERMI_FILM_INSTANT)]
        print("  R={:3d}: PLAIN={:.5f} RET={:.5f} INST={:.5f}".format(R, vals[0], vals[1], vals[2]))

def test_mobius_vs_electron():
    print("=== Нетто-заряд: электрон vs мёбиус ===")
    for name,obj in [("электрон",make_electron((0,0,0))),
                     ("мёбиус",make_mobius_particle((0,0,0)))]:
        q=sum(elem_charge(R,k,0.0) for R in obj for k in range(R['n']))
        print("  {}: sum(q)={:+.5f}".format(name, q))

def test_scaling():
    print("=== Асимптотика поля (показатель p) ===")
    Rs=[100,300,1000,3000,10000,30000]
    print("  мёбиус p={:.3f} (1/R)".format(scaling_exponent(1.0,dens_mobius,Rs)))
    print("  диполь p={:.3f} (1/R^2)".format(scaling_exponent(1.0,dens_dipole,Rs)))
    for delta in [0.5,0.2,0.05]:
        p=scaling_exponent(1.0,lambda th:kink_density(th,delta),[24,48,96,192])
        print("  кинк delta={}: p={:.2f}".format(delta, p))

def test_stability():
    print("=== Устойчивость (натяжение) ===")
    for T in [0.02,0.05,0.1,0.2,0.5]:
        a_eq,_,st=ring_equilibrium_radius(T=T)
        print("  T={:.2f}: a_eq={:.3f} a_eq*sqrt(T)={:.3f} stable={}".format(T, a_eq, a_eq*np.sqrt(T), st))

if __name__=="__main__":
    test_dipole_dipole(); test_mobius_vs_electron(); test_scaling(); test_stability()
