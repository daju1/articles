# -*- coding: utf-8 -*-
# mobius_dipole_v3B_bind.py (Python 3.5+)
# gcc -O3 -shared -fPIC mobius_dipole_v3B.c -o mobius_dipole_v3B.so -lm
import ctypes, os
_h=os.path.dirname(os.path.abspath(__file__))
lib=ctypes.CDLL(os.path.join(_h,"mobius_dipole_v3B.so"))
class MobiusRing(ctypes.Structure):
    _fields_=[("center",ctypes.c_double*3),("a",ctypes.c_double),
              ("omega",ctypes.c_double),("phase0",ctypes.c_double),
              ("delta_b",ctypes.c_double),("chirality",ctypes.c_int),
              ("mu_shear",ctypes.c_double),("nu",ctypes.c_double),("n_elem",ctypes.c_int)]
D3=ctypes.c_double*3
lib.mb_special_point.argtypes=[ctypes.POINTER(MobiusRing),ctypes.c_double,D3,D3]; lib.mb_special_point.restype=None
lib.mb_burgers_sign.argtypes=[ctypes.POINTER(MobiusRing),ctypes.c_double]; lib.mb_burgers_sign.restype=ctypes.c_double
lib.mb_pair_period.argtypes=[ctypes.POINTER(MobiusRing),ctypes.POINTER(MobiusRing)]; lib.mb_pair_period.restype=ctypes.c_double
lib.mb_pair_magnetic_dipole.argtypes=[ctypes.POINTER(MobiusRing),ctypes.POINTER(MobiusRing)]; lib.mb_pair_magnetic_dipole.restype=ctypes.c_double
lib.mb_pair_energy.argtypes=[ctypes.c_double]*4; lib.mb_pair_energy.restype=ctypes.c_double
lib.mb_twist_unwind_barrier.argtypes=[ctypes.POINTER(MobiusRing)]; lib.mb_twist_unwind_barrier.restype=ctypes.c_double

def make(center,a=1.0,omega=1.0,phase0=0.0,delta_b=1.0,chirality=1,mu_shear=1.0,nu=0.25,n_elem=48):
    R=MobiusRing(); R.center[:]=center; R.a=a; R.omega=omega; R.phase0=phase0
    R.delta_b=delta_b; R.chirality=chirality; R.mu_shear=mu_shear; R.nu=nu; R.n_elem=n_elem
    return R
def special_point(R,t):
    p=D3(); v=D3(); lib.mb_special_point(ctypes.byref(R),t,p,v); return list(p),list(v)
def burgers_sign(R,t): return lib.mb_burgers_sign(ctypes.byref(R),t)
def pair_period(N,S): return lib.mb_pair_period(ctypes.byref(N),ctypes.byref(S))
def pair_magnetic_dipole(N,S): return lib.mb_pair_magnetic_dipole(ctypes.byref(N),ctypes.byref(S))
def pair_energy(d,J,m,A): return lib.mb_pair_energy(d,J,m,A)
def twist_unwind_barrier(R): return lib.mb_twist_unwind_barrier(ctypes.byref(R))
