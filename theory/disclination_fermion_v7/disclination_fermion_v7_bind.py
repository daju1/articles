# -*- coding: utf-8 -*-
# disclination_fermion_v7_bind.py (Python 3.5+)
# gcc -O3 -shared -fPIC disclination_fermion_v7.c -o disclination_fermion_v7.so -lm
import ctypes, os
_here=os.path.dirname(os.path.abspath(__file__))
lib=ctypes.CDLL(os.path.join(_here,"disclination_fermion_v7.so"))

PLANE_T_R,PLANE_R_Z,PLANE_T_Z,PLANE_TILTED=0,1,2,3
CH_COULOMB,CH_GRAVITY,CH_EXCHANGE,CH_MAGNETIC,CH_WEAK,CH_STRONG=0,1,2,3,4,5

class Disclination(ctypes.Structure):
    _fields_=[("center",ctypes.c_double*3),
              ("a",ctypes.c_double),("omega",ctypes.c_double),("phase0",ctypes.c_double),
              ("b0",ctypes.c_double),("winding",ctypes.c_int),("chirality",ctypes.c_int),
              ("plane",ctypes.c_int),("normal_beta",ctypes.c_double),
              ("mu",ctypes.c_double),("nu",ctypes.c_double),("gamma",ctypes.c_double),
              ("n_elem",ctypes.c_int)]

D3=ctypes.c_double*3
lib.burgers_vector.argtypes=[ctypes.POINTER(Disclination),ctypes.c_double,D3]; lib.burgers_vector.restype=None
lib.frank_vector.argtypes=[ctypes.POINTER(Disclination),ctypes.c_double,D3]; lib.frank_vector.restype=None
lib.magnetic_moment.argtypes=[ctypes.POINTER(Disclination),D3]; lib.magnetic_moment.restype=None
lib.angular_momentum.argtypes=[ctypes.POINTER(Disclination),ctypes.c_double,D3]; lib.angular_momentum.restype=None
lib.gyromagnetic_g.argtypes=[ctypes.POINTER(Disclination),ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double)]; lib.gyromagnetic_g.restype=ctypes.c_double
lib.spin_value.argtypes=[ctypes.POINTER(Disclination)]; lib.spin_value.restype=ctypes.c_double
lib.twist_period_2pi_factor.argtypes=[ctypes.POINTER(Disclination)]; lib.twist_period_2pi_factor.restype=ctypes.c_int
lib.pair_force_avg.argtypes=[ctypes.POINTER(Disclination),ctypes.POINTER(Disclination),
    ctypes.c_double,ctypes.c_int,ctypes.c_int,D3]; lib.pair_force_avg.restype=None

def make(center,a=1.0,omega=0.3,phase0=0.0,b0=1.0,winding=1,chirality=1,
         plane=PLANE_T_R,normal_beta=0.0,n_elem=48,mu=1.0,nu=0.25,gamma=0.1):
    D=Disclination(); D.center[:]=center
    D.a=a;D.omega=omega;D.phase0=phase0;D.b0=b0;D.winding=winding;D.chirality=chirality
    D.plane=plane;D.normal_beta=normal_beta;D.n_elem=n_elem;D.mu=mu;D.nu=nu;D.gamma=gamma
    return D

def burgers(D,theta):
    o=D3(); lib.burgers_vector(ctypes.byref(D),theta,o); return [o[0],o[1],o[2]]
def frank(D,theta):
    o=D3(); lib.frank_vector(ctypes.byref(D),theta,o); return [o[0],o[1],o[2]]
def mag_moment(D):
    o=D3(); lib.magnetic_moment(ctypes.byref(D),o); return [o[0],o[1],o[2]]
def ang_momentum(D,c=1.0):
    o=D3(); lib.angular_momentum(ctypes.byref(D),c,o); return [o[0],o[1],o[2]]
def gyro_g(D,c=1.0):
    mz=ctypes.c_double(0); sz=ctypes.c_double(0)
    g=lib.gyromagnetic_g(ctypes.byref(D),c,ctypes.byref(mz),ctypes.byref(sz))
    return g,mz.value,sz.value
def spin(D): return lib.spin_value(ctypes.byref(D))
def period_factor(D): return lib.twist_period_2pi_factor(ctypes.byref(D))
def pair_force_avg(F1,F2,c,ch,nT=16):
    o=D3(); lib.pair_force_avg(ctypes.byref(F1),ctypes.byref(F2),c,ch,nT,o); return [o[0],o[1],o[2]]
