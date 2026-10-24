# -*- coding: utf-8 -*-
# dislocation_rings_v3_vector_bind.py (Python 3.5+)
# gcc -O3 -shared -fPIC dislocation_rings_v3_vector.c -o dislocation_rings_v3_vector.so -lm
import ctypes, os
_h=os.path.dirname(os.path.abspath(__file__))
lib=ctypes.CDLL(os.path.join(_h,"dislocation_rings_v3_vector.so"))

class RingV(ctypes.Structure):
    _fields_=[("center",ctypes.c_double*3),("a",ctypes.c_double),
              ("omega",ctypes.c_double),("phase0",ctypes.c_double),
              ("b_vec",ctypes.c_double*3),("delta_b",ctypes.c_double),
              ("mu_shear",ctypes.c_double),("nu",ctypes.c_double),
              ("n_elem",ctypes.c_int)]
D3=ctypes.c_double*3
lib.special_point.argtypes=[ctypes.POINTER(RingV),ctypes.c_double,D3,D3]; lib.special_point.restype=None
lib.pair_dipole.argtypes=[ctypes.POINTER(RingV),ctypes.POINTER(RingV),ctypes.c_double,D3]; lib.pair_dipole.restype=None
lib.pair_moments.argtypes=[ctypes.POINTER(RingV),ctypes.POINTER(RingV),ctypes.c_double,D3,D3]; lib.pair_moments.restype=None
lib.pair_magnetic_dipole.argtypes=[ctypes.POINTER(RingV),ctypes.POINTER(RingV)]; lib.pair_magnetic_dipole.restype=ctypes.c_double
lib.loop_mass.argtypes=[ctypes.POINTER(RingV),ctypes.c_double,ctypes.c_double]; lib.loop_mass.restype=ctypes.c_double
lib.pair_energy.argtypes=[ctypes.c_double]*4; lib.pair_energy.restype=ctypes.c_double
lib.total_burgers.argtypes=[ctypes.POINTER(RingV),ctypes.POINTER(RingV),D3]; lib.total_burgers.restype=None

def make(center,a=1.0,omega=1.0,phase0=0.0,b_vec=(1.0,0,0),delta_b=1.0,mu_shear=1.0,nu=0.25,n_elem=48):
    R=RingV(); R.center[:]=center; R.a=a; R.omega=omega; R.phase0=phase0
    R.b_vec[:]=b_vec; R.delta_b=delta_b; R.mu_shear=mu_shear; R.nu=nu; R.n_elem=n_elem
    return R
def special_point(R,t):
    pos=D3(); vel=D3(); lib.special_point(ctypes.byref(R),t,pos,vel); return list(pos),list(vel)
def pair_dipole(top,bot,t):
    o=D3(); lib.pair_dipole(ctypes.byref(top),ctypes.byref(bot),t,o); return list(o)
def pair_moments(top,bot,c=1.0):
    mu=D3(); J=D3(); lib.pair_moments(ctypes.byref(top),ctypes.byref(bot),c,mu,J); return list(mu),list(J)
def loop_mass(R,c=1.0,k=1.0): return lib.loop_mass(ctypes.byref(R),c,k)
def pair_energy(d,Jcons,m_eff,A): return lib.pair_energy(d,Jcons,m_eff,A)
def pair_magnetic_dipole(top,bot): return lib.pair_magnetic_dipole(ctypes.byref(top),ctypes.byref(bot))
def total_burgers(top,bot):
    o=D3(); lib.total_burgers(ctypes.byref(top),ctypes.byref(bot),o); return list(o)
