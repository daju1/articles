# -*- coding: utf-8 -*-
# binding_v4.py
# НАСТОЯЩИЙ Python-биндинг: линкует СКОМПИЛИРОВАННЫЙ dislocation_rings_v4.so
# через ctypes. Структура и сигнатуры ТОЧНО совпадают с dislocation_rings_v4.h.
# Это проверка ИМЕННО вашего C-кода (в отличие от python_implementation/,
# где формулы переписаны на Python с нуля).
#
# Сборка .so (см. build.sh):
#   gcc -O3 -shared -fPIC dislocation_rings_v4.c -o dislocation_rings_v4.so -lm
import ctypes, os

_here = os.path.dirname(os.path.abspath(__file__))
# .so ищем рядом с C-файлами (на уровень выше python_binding/)
_so = os.path.join(_here, "..", "dislocation_rings_v4.so")
lib = ctypes.CDLL(os.path.abspath(_so))

# ---- struct RingConfig: ТОЧНОЕ соответствие dislocation_rings_v4.h ----
class RingConfig(ctypes.Structure):
    _fields_ = [
        ("center",            ctypes.c_double*3),
        ("e1",                ctypes.c_double*3),
        ("e2",                ctypes.c_double*3),
        ("radius",            ctypes.c_double),
        ("omega",             ctypes.c_double),
        ("phase0",            ctypes.c_double),
        ("burgers_angle",     ctypes.c_double),   # theta_b: ориентация вектора Бюргерса
        ("burgers_magnitude", ctypes.c_double),   # |b|
        ("n_elem",            ctypes.c_int),
        ("charge_sign",       ctypes.c_int),
        ("mobius_enable",     ctypes.c_int),       # 1 -> знак cos(theta/2), период 4pi
    ]

# ---- сигнатуры (7-параметровая force, БЕЗ лишнего glide) ----
lib.compute_ring_dipole_force.argtypes = [
    ctypes.POINTER(RingConfig), ctypes.c_int,
    ctypes.POINTER(RingConfig), ctypes.c_int,
    ctypes.c_double, ctypes.c_double, ctypes.c_int,
    ctypes.c_double*3]
lib.compute_ring_dipole_force.restype = None

lib.ring_element_burgers.argtypes = [
    ctypes.POINTER(RingConfig), ctypes.c_int, ctypes.c_double, ctypes.c_double*3]
lib.ring_element_burgers.restype = None

lib.ring_total_energy.argtypes = [
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.POINTER(ctypes.c_double), ctypes.POINTER(ctypes.c_double)]
lib.ring_total_energy.restype = ctypes.c_double

lib.ring_equilibrium_radius.argtypes = [
    ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
    ctypes.c_double, ctypes.c_double, ctypes.c_int, ctypes.POINTER(ctypes.c_int)]
lib.ring_equilibrium_radius.restype = ctypes.c_double

# Режимы Ферми (как в заголовке)
FERMI_PLAIN, FERMI_FILM_RETARDED, FERMI_FILM_INSTANT, FERMI_NONE = 0, 1, 2, 3

def make_ring_b(center, radius, burgers_angle, burgers_magnitude,
                n_elem=48, sign=1, omega=0.0, phase0=0.0, mobius=0):
    """Создать кольцо v4 (векторный Бюргерс). Сигнатура совпадает с тем, что
    зовут ноутбуки: make_ring_b(center, radius, burgers_angle, burgers_magnitude, ...)."""
    r = RingConfig()
    r.center[:] = center
    r.e1[:] = (1.0, 0.0, 0.0)
    r.e2[:] = (0.0, 1.0, 0.0)
    r.radius = radius
    r.omega = omega
    r.phase0 = phase0
    r.burgers_angle = burgers_angle
    r.burgers_magnitude = burgers_magnitude
    r.n_elem = n_elem
    r.charge_sign = sign
    r.mobius_enable = mobius
    return r

def make_dipole(center, a=1.0, omega=0.3, n=24, bmag=1.0, dphi=0.0, mobius=0,
                sep=0.3):
    """Диполь = два коаксиальных кольца противоположного знака вдоль z."""
    cN = (center[0], center[1], center[2]+sep)
    cS = (center[0], center[1], center[2]-sep)
    return [make_ring_b(cN, a, 0.0, bmag, n_elem=n, sign=+1, omega=omega,
                        phase0=0.0, mobius=mobius),
            make_ring_b(cS, a, 0.0, bmag, n_elem=n, sign=-1, omega=omega,
                        phase0=dphi, mobius=mobius)]

def _force(field_rings, source_rings, t, c, mode):
    F = (ctypes.c_double*3)()
    FR = (RingConfig*len(field_rings))(*field_rings)
    SR = (RingConfig*len(source_rings))(*source_rings)
    lib.compute_ring_dipole_force(FR, len(field_rings), SR, len(source_rings),
                                  t, c, mode, F)
    return [F[0], F[1], F[2]]

# Имена-алиасы под ноутбуки:
def force_v4(field_rings, source_rings, t, c, mode):
    """7-параметровая сила (mode без glide). Так ДОЛЖЕН зваться v4.ipynb."""
    return _force(field_rings, source_rings, t, c, mode)

def force_multi(field_rings, source_rings, t, c, mode):
    """Алиас, который зовут ноутбуки для диполь-диполь теста."""
    return _force(field_rings, source_rings, t, c, mode)

def burgers_of(ring, k, t):
    b = (ctypes.c_double*3)()
    lib.ring_element_burgers(ctypes.byref(ring), k, t, b)
    return [b[0], b[1], b[2]]

def equilibrium_radius(q=1.0, r0=1e-3, k=1.0, T=0.05, amin=1e-3, amax=50.0, n=20000):
    st = ctypes.c_int(0)
    a_eq = lib.ring_equilibrium_radius(q, r0, k, T, amin, amax, n, ctypes.byref(st))
    return a_eq, bool(st.value)

def total_energy(a, q=1.0, r0=1e-3, k=1.0, T=0.05):
    ec = ctypes.c_double(0); et = ctypes.c_double(0)
    E = lib.ring_total_energy(a, q, r0, k, T, ctypes.byref(ec), ctypes.byref(et))
    return E, ec.value, et.value
