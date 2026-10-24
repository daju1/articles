# -*- coding: utf-8 -*-
# impl_v4_field.py
# Чистый независимая Python-реализация dislocation_rings_v4.c (ядро поля + каналы central/vortex
# + режимы Ферми + мёбиусов знак). Используется для проверки логики C без
# компиляции (on-device песочница блокирует ctypes/компилятор).
# c=1 в безразмерных единицах, если не указано иначе.
import numpy as np

FERMI_PLAIN, FERMI_FILM_RETARDED, FERMI_FILM_INSTANT, FERMI_NONE = 0, 1, 2, 3

def ring_kin(R, k, t):
    """Положение/скорость/ускорение k-го элемента кольца R в момент t."""
    th = R['phase0'] + 2*np.pi*k/R['n'] + R['omega']*t
    c1, s1 = np.cos(th), np.sin(th)
    radial = R['e1']*c1 + R['e2']*s1
    pos = R['center'] + R['radius']*radial
    vel = R['radius']*R['omega']*(-R['e1']*s1 + R['e2']*c1)
    acc = -R['omega']**2 * R['radius']*radial
    return pos, vel, acc

def mobius_sign(R, theta):
    """Мёбиусов знак величины заряда: cos(theta/2), период 4pi."""
    return np.cos(0.5*theta) if R.get('mobius', 0) else 1.0

def elem_charge(R, k, t):
    """Скалярная величина заряда k-го элемента (central-канал)."""
    th = R['phase0'] + 2*np.pi*k/R['n'] + R['omega']*t
    return R['bmag']*R['sign']*mobius_sign(R, th)/R['n']

def elem_burgers(R, k, t):
    """Полный вектор Бюргерса k-го элемента (ориентация поворачивается с кольцом)."""
    th = R['phase0'] + 2*np.pi*k/R['n'] + R['omega']*t
    ang = th + R.get('bangle', 0.0)
    mag = R['bmag']*R['sign']*mobius_sign(R, th)
    return mag*np.cos(ang)*R['e1'] + mag*np.sin(ang)*R['e2']

def retarded_time(fp, tf, R, k, c):
    pos, vel, acc = ring_kin(R, k, tf)
    Rd = np.linalg.norm(fp-pos); tp = tf - Rd/c
    for _ in range(20):
        pos, vel, acc = ring_kin(R, k, tp)
        diff = fp-pos; Rn = np.linalg.norm(diff)
        f = tp + Rn/c - tf
        vr = np.dot(diff, vel)/(Rn+1e-30); df = 1.0 - vr/c
        if abs(df) < 1e-12: break
        d = -f/df; tp += d
        if abs(d) < 1e-14: break
    return tp

def field_from_elem(fp, R, k, tf, c, a_extra, mode):
    """Поле элемента: central (кулон) + vortex (аналог B). Соответствует
    compute_field_from_element из v4.c."""
    tp = retarded_time(fp, tf, R, k, c)
    pos, vel, acc = ring_kin(R, k, tp)
    Rvec = fp-pos; Rd = np.linalg.norm(Rvec)
    if Rd < 1e-9: return np.zeros(3), np.zeros(3)
    n = Rvec/Rd
    Rstar = Rd*(1 - np.dot(vel, n)/c); Rstar = max(Rstar, 1e-9)
    a_eff = acc.copy()
    if mode != FERMI_NONE and a_extra is not None:
        a_eff = a_eff + a_extra
    fermi = 1.0 + np.dot(a_eff, Rvec)/c**2 if mode != FERMI_NONE else 1.0
    v2c2 = np.dot(vel, vel)/c**2; Rn_a = np.dot(Rvec, acc)/c**2
    common = 1.0/Rstar**3; bracket = fermi*(1 + Rn_a - v2c2)
    central = common*(bracket*Rvec - bracket*Rd*vel/c - fermi*acc*Rstar*Rd/c**2)
    vortex = np.cross(n, central)
    return central, vortex

def film_accel(fp, teval, rings, c):
    """Мгновенная/запаздывающая плёнка: взвешенное по 1/d среднее ускорений
    чужих элементов. Соответствует film_effective_accel."""
    wsum = 0.0; asum = np.zeros(3)
    for R in rings:
        for k in range(R['n']):
            pos, vel, acc = ring_kin(R, k, teval)
            d = np.linalg.norm(fp-pos); d = max(d, 1e-9); w = 1.0/d
            wsum += w; asum += w*acc
    return asum/wsum if wsum > 1e-30 else np.zeros(3)

def force(field_rings, source_rings, t, c, mode=FERMI_NONE):
    """Полная сила между наборами колец. Соответствует compute_ring_dipole_force."""
    F = np.zeros(3)
    for FR in field_rings:
        for fk in range(FR['n']):
            qf = elem_charge(FR, fk, t)
            fpos, fvel, _ = ring_kin(FR, fk, t)
            a_extra = None
            if mode == FERMI_FILM_INSTANT:
                a_extra = film_accel(fpos, t, source_rings, c)
            elif mode == FERMI_FILM_RETARDED:
                a_extra = film_accel(fpos, t, source_rings, c)  # упрощённо; точно — на tp
            for SR in source_rings:
                for sk in range(SR['n']):
                    qs = elem_charge(SR, sk, t)
                    central, vortex = field_from_elem(fpos, SR, sk, t, c, a_extra, mode)
                    vxB = np.cross(fvel, vortex)/c
                    F += qf*qs*central + qf*vxB
    return F

def make_ring(center, a=1.0, omega=0.0, n=24, bangle=0.0, bmag=1.0, sign=1,
              phase0=0.0, mobius=0):
    return dict(center=np.array(center, float),
                e1=np.array([1., 0, 0]), e2=np.array([0, 1., 0]),
                radius=a, omega=omega, phase0=phase0,
                bangle=bangle, bmag=bmag, n=n, sign=sign, mobius=mobius)

def make_electron(center, a=1.0, omega=0.2, n=24, dphi=0.0):
    """Электрон = диполь из двух колец противоположного знака (mobius off)."""
    return [make_ring((center[0], center[1], center[2]+0.3), a, omega, n, sign=+1),
            make_ring((center[0], center[1], center[2]-0.3), a, omega, n, sign=-1,
                      phase0=dphi)]

def make_mobius_particle(center, a=1.0, omega=0.2, n=24, phase0=0.0):
    """Мёбиус-частица (Уровень 1) = одно кольцо с мёбиусовым знаком (период 4pi)."""
    return [make_ring(center, a, omega, n, sign=+1, phase0=phase0, mobius=1)]
