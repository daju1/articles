# -*- coding: utf-8 -*-
# disclination_fermion_v6_bind.py (Python 3.5+)
# Сборка: gcc -O3 -shared -fPIC disclination_fermion_v6.c -o disclination_fermion_v6.so -lm
import ctypes, os
_here=os.path.dirname(os.path.abspath(__file__))
lib=ctypes.CDLL(os.path.join(_here,"disclination_fermion_v6.so"))

class Disclination(ctypes.Structure):
    _fields_=[("center",ctypes.c_double*3),("e1",ctypes.c_double*3),("e2",ctypes.c_double*3),
              ("a",ctypes.c_double),("omega",ctypes.c_double),("phase0",ctypes.c_double),
              ("b0",ctypes.c_double),("mu",ctypes.c_double),("nu",ctypes.c_double),
              ("gamma",ctypes.c_double),("n_elem",ctypes.c_int)]

CH_COULOMB,CH_EXCHANGE,CH_GRAVITY=0,1,2

lib.burgers_twist.argtypes=[ctypes.POINTER(Disclination),ctypes.c_int,ctypes.c_double]
lib.burgers_twist.restype=ctypes.c_double
lib.twist_multipoles.argtypes=[ctypes.POINTER(Disclination),ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double)]
lib.twist_multipoles.restype=None
lib.disclination_spin_ratio.argtypes=[ctypes.POINTER(Disclination)]
lib.disclination_spin_ratio.restype=ctypes.c_double
lib.disclination_g_factor.argtypes=[ctypes.POINTER(Disclination)]
lib.disclination_g_factor.restype=ctypes.c_double
lib.disclination_pair_force.argtypes=[ctypes.POINTER(Disclination),ctypes.POINTER(Disclination),
    ctypes.c_double,ctypes.c_double,ctypes.c_int,ctypes.c_double*3]
lib.disclination_pair_force.restype=None
lib.disclination_pair_force_avg.argtypes=[ctypes.POINTER(Disclination),ctypes.POINTER(Disclination),
    ctypes.c_double,ctypes.c_int,ctypes.c_int,ctypes.c_double*3]
lib.disclination_pair_force_avg.restype=None

def make(center,a=1.0,omega=0.3,phase0=0.0,b0=1.0,n_elem=48,mu=1.0,nu=0.25,gamma=0.1):
    D=Disclination(); D.center[:]=center; D.e1[:]=(1.0,0.0,0.0); D.e2[:]=(0.0,1.0,0.0)
    D.a=a; D.omega=omega; D.phase0=phase0; D.b0=b0; D.n_elem=n_elem
    D.mu=mu; D.nu=nu; D.gamma=gamma; return D

def multipoles(D,t=0.0):
    m=ctypes.c_double(0); d=ctypes.c_double(0)
    lib.twist_multipoles(ctypes.byref(D),t,ctypes.byref(m),ctypes.byref(d))
    return m.value,d.value
def spin_ratio(D): return lib.disclination_spin_ratio(ctypes.byref(D))
def g_factor(D):   return lib.disclination_g_factor(ctypes.byref(D))
def pair_force(F1,F2,t,c,ch):
    o=(ctypes.c_double*3)()
    lib.disclination_pair_force(ctypes.byref(F1),ctypes.byref(F2),t,c,ch,o)
    return [o[0],o[1],o[2]]

def pair_force_avg(F1,F2,c,ch,nT=16):
    o=(ctypes.c_double*3)()
    lib.disclination_pair_force_avg(ctypes.byref(F1),ctypes.byref(F2),c,ch,nT,o)
    return [o[0],o[1],o[2]]
