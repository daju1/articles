# -*- coding: utf-8 -*-
# impl_v4_stability.py
# независимая Python-реализация функций устойчивости из dislocation_rings_v4.c
# (ring_total_energy, ring_equilibrium_radius): баланс кулоновской
# самоэнергии и линейного натяжения T=mu*b^2.
import numpy as np

def ring_total_energy(a, q=1.0, r0=1e-3, k=1.0, T=0.05):
    """E(a) = (k q^2/(pi a)) ln(8a/r0) + T*2*pi*a. Возвращает (E,E_coul,E_tens)."""
    E_coul = (k*q*q/(np.pi*a))*np.log(8.0*a/r0)
    E_tens = T*2.0*np.pi*a
    return E_coul+E_tens, E_coul, E_tens

def ring_equilibrium_radius(q=1.0, r0=1e-3, k=1.0, T=0.05, amin=1e-3, amax=50.0, n=20000):
    ag = np.linspace(amin, amax, n)
    E = np.array([ring_total_energy(a, q, r0, k, T)[0] for a in ag])
    i = int(np.argmin(E))
    stable = 1 if (0 < i < n-1 and (E[i+1]-2*E[i]+E[i-1]) > 0) else 0
    return ag[i], E[i], stable

def analytic_a_eq(q=1.0, k=1.0, T=0.05):
    """a_eq = sqrt(k q^2/(2 pi T)) (без лог-поправки)."""
    return np.sqrt(k*q*q/(2.0*np.pi*T))

if __name__ == "__main__":
    a_eq, E, st = ring_equilibrium_radius()
    print("a_eq={:.4f} E={:.4f} stable={} (аналит~{:.4f})".format(a_eq, E, st, analytic_a_eq()))
    for T in [0.02, 0.05, 0.1, 0.2, 0.5]:
        a_eq, _, _ = ring_equilibrium_radius(T=T)
        print("  T={:.2f}: a_eq={:.3f}  a_eq*sqrt(T)={:.3f}".format(T, a_eq, a_eq*np.sqrt(T)))
