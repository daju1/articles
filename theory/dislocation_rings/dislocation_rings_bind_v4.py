# dislocation_rings_bind_v4.py
# Биндинг v4: векторный заряд Бюргерса + мёбиусов знак.
# Структура ТОЧНО соответствует dislocation_rings_v4.h (никакого ABI-сдвига).
import ctypes, math

lib = ctypes.CDLL("./dislocation_rings_v4.so")

class RingConfig(ctypes.Structure):
    _fields_ = [
        ("center", ctypes.c_double*3),
        ("e1", ctypes.c_double*3),
        ("e2", ctypes.c_double*3),
        ("radius", ctypes.c_double),
        ("omega", ctypes.c_double),
        ("phase0", ctypes.c_double),
        ("burgers_angle", ctypes.c_double),
        ("burgers_magnitude", ctypes.c_double),
        ("n_elem", ctypes.c_int),
        ("charge_sign", ctypes.c_int),
        ("mobius_enable", ctypes.c_int)]

lib.compute_ring_dipole_force.argtypes = [
    ctypes.POINTER(RingConfig), ctypes.c_int,
    ctypes.POINTER(RingConfig), ctypes.c_int,
    ctypes.c_double, ctypes.c_double,
    ctypes.c_int,
    ctypes.c_double*3]
lib.compute_ring_dipole_force.restype = None

lib.ring_element_burgers.argtypes = [
    ctypes.POINTER(RingConfig), ctypes.c_int, ctypes.c_double, ctypes.c_double*3]
lib.ring_element_burgers.restype = None

def make_ring_b(center, radius, burgers_angle, burgers_magnitude,
                n_elem=48, sign=1, omega=0.0, phase0=0.0, mobius=0):
    r = RingConfig()
    r.center[:] = center
    r.e1[:] = (1.0,0.0,0.0)
    r.e2[:] = (0.0,1.0,0.0)
    r.radius = radius
    r.omega = omega
    r.phase0 = phase0
    r.burgers_angle = burgers_angle
    r.burgers_magnitude = burgers_magnitude
    r.n_elem = n_elem
    r.charge_sign = sign
    r.mobius_enable = mobius
    return r

def make_dipole(center, a=1.0, omega=0.3, n=24, bmag=1.0, dphi=0.0, mobius=0):
    """Диполь = два коаксиальных кольца противоположного знака вдоль z."""
    import numpy as np
    cN = (center[0], center[1], center[2]+0.3)
    cS = (center[0], center[1], center[2]-0.3)
    rN = make_ring_b(cN, a, 0.0, bmag, n_elem=n, sign=+1, omega=omega, phase0=0.0, mobius=mobius)
    rS = make_ring_b(cS, a, 0.0, bmag, n_elem=n, sign=-1, omega=omega, phase0=dphi, mobius=mobius)
    return [rN, rS]

def force_v4(field_rings, source_rings, t, c, mode):
    F=(ctypes.c_double*3)()
    FR=(RingConfig*len(field_rings))(*field_rings)
    SR=(RingConfig*len(source_rings))(*source_rings)
    lib.compute_ring_dipole_force(FR,len(field_rings),SR,len(source_rings),t,c,mode,F)
    return [F[0],F[1],F[2]]

def burgers_of(ring, k, t):
    b=(ctypes.c_double*3)()
    lib.ring_element_burgers(ctypes.byref(ring), k, t, b)
    return [b[0],b[1],b[2]]

# Режимы Ферми (как в v3)
FERMI_PLAIN, FERMI_FILM_RETARDED, FERMI_FILM_INSTANT, FERMI_NONE = 0,1,2,3
