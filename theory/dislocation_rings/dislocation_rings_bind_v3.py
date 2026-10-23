# dislocation_rings_bind_v3_fixed.py
# ИСПРАВЛЕННЫЙ биндинг v3: структура ТОЧНО соответствует dislocation_rings.h
# (скалярный total_charge, 7-параметровая сигнатура, без glide).
# Оригинал dislocation_rings_bind_v3.py оставлен как задел под v4
# (там burgers_* поля) — см. dislocation_rings_bind_v4.py для рабочей v4.

import ctypes, math
lib = ctypes.CDLL("dislocation_rings_v3.so")

class RingConfig(ctypes.Structure):
    _fields_ = [
        ("center", ctypes.c_double*3),
        ("e1", ctypes.c_double*3),
        ("e2", ctypes.c_double*3),
        ("radius", ctypes.c_double),
        ("omega", ctypes.c_double),
        ("phase0", ctypes.c_double),
        ("total_charge", ctypes.c_double),   # <-- совпадает с .h (было burgers_angle+magnitude)
        ("n_elem", ctypes.c_int),
        ("charge_sign", ctypes.c_int)]

lib.compute_ring_dipole_force.argtypes = [
    ctypes.POINTER(RingConfig), ctypes.c_int,
    ctypes.POINTER(RingConfig), ctypes.c_int,
    ctypes.c_double, ctypes.c_double,
    ctypes.c_int,                              # <-- один int (mode), без лишнего glide
    ctypes.c_double*3]
    
lib.compute_ring_dipole_force.restype = None

def make_ring_b(
        center,
        radius,
        total_charge,
        n_elem=48,
        sign=1,
        omega=0.0,
        phase0=0.0):
    r = RingConfig()
    r.center[:] = center;
    r.e1[:] = (1.0,0.0,0.0);
    r.e2[:] = (0.0,1.0,0.0)
    r.radius = radius;
    r.omega = omega;
    r.phase0 = phase0
    r.total_charge = total_charge
    r.n_elem = n_elem;
    r.charge_sign = sign
    return r

def force_v3(field_rings, source_rings, t, c, mode):
    F=(ctypes.c_double*3)()
    FR=(RingConfig*len(field_rings))(*field_rings)
    SR=(RingConfig*len(source_rings))(*source_rings)
    lib.compute_ring_dipole_force(
        FR,len(field_rings),
        SR,len(source_rings),
        t,c,
        mode,
        F)
    return [F[0],F[1],F[2]]

FERMI_PLAIN, FERMI_FILM_RETARDED, FERMI_FILM_INSTANT, FERMI_NONE = 0,1,2,3
