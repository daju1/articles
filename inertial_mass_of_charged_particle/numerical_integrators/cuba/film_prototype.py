#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
film_prototype.py: плёнка, скользящая по мировым трубкам двух зарядов (c = 1).
Только стандартная библиотека, совместим с Python 3.5.

Модель: плёнка (позвоночник) - кубическая кривая Эрмита между точками P1, P2 на
двух трубках. Касательные в концах жёстко ортогональны 4-скорости трубки.
Энергия = длина кривой в метрике (+,+,-). Положение на трубках - СТАЦИОНАРНАЯ
точка энергии (не минимум, см. F'' в тестах) при заданной фазе.
Калибровка T: T = (t1 + t2)/2, варьируется d = t2 - t1.
Калибровка A: t1 задано, варьируется t2.

PYTHONIOENCODING=utf-8 python3 ./film_prototype.py
"""
import math


def dot(a, b):                         # метрика (x, y, t): (+, +, -)
    return a[0]*b[0] + a[1]*b[1] - a[2]*b[2]


def axpy(a, b, k):                     # a + k*b
    return (a[0] + k*b[0], a[1] + k*b[1], a[2] + k*b[2])


class Tube(object):
    def __init__(self, pos, vel):
        self.pos, self.vel = pos, vel

    def point(self, t):
        x, y = self.pos(t)
        return (x, y, t)

    def gamma(self, t):
        vx, vy = self.vel(t)
        return 1.0/math.sqrt(1.0 - vx*vx - vy*vy)

    def u(self, t):                    # 4-скорость, u.u = -1
        vx, vy = self.vel(t)
        g = self.gamma(t)
        return (g*vx, g*vy, g)


def hyperbolic(g, shift=0.0):          # равноускоренное движение, x(0) = 1/g + shift
    r2 = 1.0/(g*g)
    return Tube(lambda t: (math.sqrt(r2 + t*t) + shift, 0.0),
                lambda t: (t/math.sqrt(r2 + t*t), 0.0))


def circular(rho, omega, phi0):
    return Tube(lambda t: (rho*math.cos(omega*t + phi0), rho*math.sin(omega*t + phi0)),
                lambda t: (-rho*omega*math.sin(omega*t + phi0),
                           rho*omega*math.cos(omega*t + phi0)))


def film_length(P1, u1, P2, u2, n=200):
    """Длина кубической плёнки, ортогональной u1 в P1 и u2 в P2 (n чётное)."""
    D = (P2[0]-P1[0], P2[1]-P1[1], P2[2]-P1[2])
    DD = dot(D, D)
    if DD <= 0.0:
        return float('nan')            # хорда не пространственноподобна
    S = math.sqrt(DD)
    T1 = axpy(D, u1, dot(D, u1))       # проекция хорды на u1-перпендикуляр
    T2 = axpy(D, u2, dot(D, u2))
    q1, q2 = dot(T1, T1), dot(T2, T2)
    if q1 <= 0.0 or q2 <= 0.0:
        return float('nan')
    t1 = tuple(S*c/math.sqrt(q1) for c in T1)      # |касательная| = длина хорды
    t2 = tuple(S*c/math.sqrt(q2) for c in T2)
    h = 1.0/n
    total = 0.0
    for k in range(n + 1):
        l = k*h
        a = 6.0*l - 6.0*l*l
        b = 3.0*l*l - 4.0*l + 1.0
        c = 3.0*l*l - 2.0*l
        dC = tuple(a*D[i] + b*t1[i] + c*t2[i] for i in range(3))
        q = dot(dC, dC)
        if q <= 0.0:
            return float('nan')
        w = 1.0 if (k == 0 or k == n) else (4.0 if k % 2 == 1 else 2.0)
        total += w*math.sqrt(q)
    return total*h/3.0


def energy(A, B, t1, t2):
    return film_length(A.point(t1), A.u(t1), B.point(t2), B.u(t2))


def stationary(f, x0, h=3e-4, tol=1e-9, itmax=80, maxstep=0.05):
    """Корень f'(x) = 0 (Ньютон, конечные разности, ограничение шага). x, f''(x)."""
    x = x0
    for _ in range(itmax):
        g = (f(x + h) - f(x - h))/(2.0*h)
        g2 = (f(x + 2*h) - 2.0*f(x) + f(x - 2*h))/(4.0*h*h)
        dx = -g/g2
        dx = max(-maxstep, min(maxstep, dx))
        x += dx
        if abs(dx) < tol:
            g2 = (f(x + 2*h) - 2.0*f(x) + f(x - 2*h))/(4.0*h*h)
            return x, g2
    raise RuntimeError("Ньютон не сошёлся")


def solve_T(A, B, T, d0=0.0):
    f = lambda d: energy(A, B, T - 0.5*d, T + 0.5*d)
    return stationary(f, d0)


def march_T(A, B, T, steps=40):
    """Идём по T от 0 до T, на каждом шаге стартуем с предыдущего решения."""
    d = 0.0
    for k in range(1, steps + 1):
        d, _ = solve_T(A, B, T*k/steps, d)
    return d


def sliding_T(A, B, T, dT=1e-2):
    """ds2/ds1, dt1/dT, dt2/dT вдоль семейства плёнок (калибровка T)."""
    d = march_T(A, B, T)
    dm, _ = solve_T(A, B, T - dT, d)
    dp, _ = solve_T(A, B, T + dT, d)
    dt1 = 1.0 - 0.5*(dp - dm)/(2*dT)
    dt2 = 1.0 + 0.5*(dp - dm)/(2*dT)
    t1, t2 = T - 0.5*d, T + 0.5*d
    return (dt2/B.gamma(t2))/(dt1/A.gamma(t1)), dt1, dt2, d


def march_A(A, B, t1, steps=40):
    t2 = 0.0
    for k in range(1, steps + 1):
        tk = t1*k/steps
        t2, _ = stationary(lambda s: energy(A, B, tk, s), t2)
    return t2


def sliding_A(A, B, t1, dt=1e-2):
    t2 = march_A(A, B, t1)
    a, _ = stationary(lambda s: energy(A, B, t1 - dt, s), t2)
    b, _ = stationary(lambda s: energy(A, B, t1 + dt, s), t2)
    dt2 = (b - a)/(2*dt)
    return dt2*A.gamma(t1)/B.gamma(t2), dt2, t2


if __name__ == "__main__":
    print("T1 Риндлер g1=1.0, g2=0.8: ожидание ds2/ds1 = g1/g2 = 1.25, F = 0.25")
    A, B = hyperbolic(1.0), hyperbolic(0.8)
    for T in (0.3, 1.5):
        r, a, b, d = sliding_T(A, B, T)
        F = energy(A, B, T - d/2, T + d/2)
        _, g2 = solve_T(A, B, T, d)
        print("  T=%.1f ds2/ds1=%.6f dt1/dT=%.5f dt2/dT=%.5f F=%.8f F''=%.3f"
              % (T, r, a, b, F, g2))
    print("T2 Белл: g=1 у обоих, сдвиг 0.5 (не жёсткое тело)")
    A, B = hyperbolic(1.0), hyperbolic(1.0, 0.5)
    for T in (0.3, 1.0):
        r, a, b, d = sliding_T(A, B, T)
        print("  калибровка T: T=%.1f ds2/ds1=%.5f dt2/dT=%.5f d=%.5f" % (T, r, b, d))
    for t1 in (0.3, 1.0):
        r, dt2, t2 = sliding_A(A, B, t1)
        print("  калибровка A: t1=%.1f ds2/ds1=%.5f dt2/dt1=%.5f t2-t1=%.5f"
              % (t1, r, dt2, t2 - t1))
    print("T3 двойная звезда: rho1=1.0, rho2=0.5, omega=0.4 (v1=0.4, v2=0.2)")
    A, B = circular(1.0, 0.4, 0.0), circular(0.5, 0.4, math.pi)
    g1, g2_ = A.gamma(0.0), B.gamma(0.0)
    print("  ожидание (винтовая симметрия): dt_i/dT = 1, ds2/ds1 = g1/g2 = %.6f" % (g1/g2_))
    for T in (0.0, 1.7, 5.0):
        r, a, b, d = sliding_T(A, B, T)
        _, F2 = solve_T(A, B, T, d)
        print("  T=%.1f ds2/ds1=%.6f dt1/dT=%.5f dt2/dT=%.5f d*=%.6f F=%.8f F''=%.4f"
              % (T, r, a, b, d, energy(A, B, T - d/2, T + d/2), F2))