# -*- coding: utf-8 -*-
# binding_v3.py
# Настоящий ctypes-биндинг для dislocation_rings_v3.so (И для базового
# dislocation_rings.so — у них ОДИНАКОВАЯ структура из dislocation_rings.h:
# скалярный total_charge, 7-параметровая compute_ring_dipole_force, БЕЗ glide).
#
# Сборка:
#   gcc -O3 -shared -fPIC dislocation_rings_v3.c -o dislocation_rings_v3.so -lm
import ctypes, os

def _load(so_name):
    here=os.path.dirname(os.path.abspath(__file__))
    return ctypes.CDLL(os.path.abspath(os.path.join(here,"..",so_name)))

# по умолчанию грузим v3; для базового передайте so="dislocation_rings.so"
class RingConfig(ctypes.Structure):
    _fields_ = [
        ("center",       ctypes.c_double*3),
        ("e1",           ctypes.c_double*3),
        ("e2",           ctypes.c_double*3),
        ("radius",       ctypes.c_double),
        ("omega",        ctypes.c_double),
        ("phase0",       ctypes.c_double),
        ("total_charge", ctypes.c_double),   # СКАЛЯРНЫЙ заряд (не burgers_*!)
        ("n_elem",       ctypes.c_int),
        ("charge_sign",  ctypes.c_int),
    ]

FERMI_PLAIN, FERMI_FILM_RETARDED, FERMI_FILM_INSTANT, FERMI_NONE = 0, 1, 2, 3

class Lib:
    """Обёртка над одной .so (v3 или базовая)."""
    def __init__(self, so_name="dislocation_rings_v3.so"):
        self.lib=_load(so_name)
        self.lib.compute_ring_dipole_force.argtypes=[
            ctypes.POINTER(RingConfig), ctypes.c_int,
            ctypes.POINTER(RingConfig), ctypes.c_int,
            ctypes.c_double, ctypes.c_double, ctypes.c_int,
            ctypes.c_double*3]
        self.lib.compute_ring_dipole_force.restype=None

    def make_ring(self, center, radius, total_charge, n_elem=48, sign=1,
                  omega=0.0, phase0=0.0):
        r=RingConfig()
        r.center[:]=center; r.e1[:]=(1.,0.,0.); r.e2[:]=(0.,1.,0.)
        r.radius=radius; r.omega=omega; r.phase0=phase0
        r.total_charge=total_charge; r.n_elem=n_elem; r.charge_sign=sign
        return r

    # алиас под ноутбук v3, где make_ring_b(center, radius, total_charge, sign=...)
    def make_ring_b(self, center, radius, total_charge, n_elem=48, sign=1,
                    omega=0.0, phase0=0.0):
        return self.make_ring(center, radius, total_charge, n_elem, sign, omega, phase0)

    def make_dipole(self, center, a=1.0, omega=0.3, n=24, q=1.0, dphi=0.0, sep=0.3):
        cN=(center[0],center[1],center[2]+sep); cS=(center[0],center[1],center[2]-sep)
        return [self.make_ring(cN,a,q,n_elem=n,sign=+1,omega=omega,phase0=0.0),
                self.make_ring(cS,a,q,n_elem=n,sign=-1,omega=omega,phase0=dphi)]

    def force(self, field_rings, source_rings, t, c, mode):
        F=(ctypes.c_double*3)()
        FR=(RingConfig*len(field_rings))(*field_rings)
        SR=(RingConfig*len(source_rings))(*source_rings)
        self.lib.compute_ring_dipole_force(FR,len(field_rings),SR,len(source_rings),
                                           t,c,mode,F)
        return [F[0],F[1],F[2]]
    # алиасы под ноутбуки
    force_v3=force
    force_multi=force

# Готовые экземпляры
v3  = Lib("dislocation_rings_v3.so")
base= Lib("dislocation_rings.so")
