# -*- coding: utf-8 -*-
# impl_mobius_analysis.py
# независимая Python-реализация анализа мёбиусовой плотности: асимптотика поля, мультиполи,
# регуляризация особой точки кинком. Подтверждает закон 1/R (p~1.00).
import numpy as np
from numpy import pi

def ring_potential(R_obs, a, density_fn, n=20000, axis='x'):
    """Статический потенциал кольца радиуса a с плотностью density_fn(theta)
    в точке на оси axis на расстоянии R_obs."""
    th = np.linspace(0, 2*pi, n, endpoint=False)
    xs = a*np.cos(th); ys = a*np.sin(th); zs = 0*th
    obs = np.array([R_obs, 0, 0.]) if axis == 'x' else np.array([0, 0, R_obs])
    d = np.sqrt((obs[0]-xs)**2 + (obs[1]-ys)**2 + (obs[2]-zs)**2)
    return np.sum(density_fn(th)/d)*(2*pi/n)

dens_mobius = lambda th: np.cos(0.5*th)       # период 4pi (спинор)
dens_dipole = lambda th: np.cos(th)           # обычный диполь (2pi)

def kink_density(th, delta):
    """Регуляризованный скачок знака (кинк) ширины delta в особой точке."""
    return np.tanh((th-pi)/delta)

def scaling_exponent(a, density_fn, Rs):
    V = np.array([abs(ring_potential(R, a, density_fn)) for R in Rs])
    A = np.vstack([np.ones(len(Rs)), np.log(Rs)]).T
    return -np.linalg.lstsq(A, np.log(V), rcond=None)[0][1]

if __name__ == "__main__":
    Rs = [100, 300, 1000, 3000, 10000, 30000]
    print("мёбиус p = {:.3f} (ожид. 1.0)".format(scaling_exponent(1.0, dens_mobius, Rs)))
    print("диполь p = {:.3f} (ожид. 2.0)".format(scaling_exponent(1.0, dens_dipole, Rs)))
    for delta in [0.5, 0.2, 0.05]:
        p = scaling_exponent(1.0, lambda th: kink_density(th, delta), [24,48,96,192])
        print("  кинк delta={}: p={:.2f}".format(delta, p))
