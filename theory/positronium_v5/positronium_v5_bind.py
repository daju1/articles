# -*- coding: utf-8 -*-
# positronium_v5_bind.py — ctypes-биндинг к positronium_v5.so (Python 3.5+).
# Сборка: gcc -O3 -shared -fPIC positronium_v5.c -o positronium_v5.so -lm
import ctypes, os
_here=os.path.dirname(os.path.abspath(__file__))
lib=ctypes.CDLL(os.path.join(_here,"positronium_v5.so"))

class Monopole(ctypes.Structure):
    _fields_=[("center",ctypes.c_double*3),("e1",ctypes.c_double*3),("e2",ctypes.c_double*3),
              ("a",ctypes.c_double),("omega",ctypes.c_double),("phase0",ctypes.c_double),
              ("qm",ctypes.c_double),("sign",ctypes.c_int),
              ("mu",ctypes.c_double),("nu",ctypes.c_double),("gamma",ctypes.c_double)]
class Fermion(ctypes.Structure):
    _fields_=[("plus",Monopole),("minus",Monopole)]

# задачи и режимы заряда
TASK_STABILITY,TASK_COULOMB,TASK_GRAVITY,TASK_EXCHANGE,TASK_BOHR=0,1,2,3,4
ECHARGE_DARWIN,ECHARGE_FERMI=0,1

lib.fermion_self_energy.argtypes=[ctypes.POINTER(Fermion),ctypes.c_double,
    ctypes.POINTER(ctypes.c_double),ctypes.POINTER(ctypes.c_double)]
lib.fermion_self_energy.restype=ctypes.c_double
lib.fermion_equilibrium_radius.argtypes=[ctypes.POINTER(Fermion),ctypes.c_double,
    ctypes.c_double,ctypes.c_int,ctypes.POINTER(ctypes.c_int)]
lib.fermion_equilibrium_radius.restype=ctypes.c_double
lib.fermion_pair_force.argtypes=[ctypes.POINTER(Fermion),ctypes.POINTER(Fermion),
    ctypes.c_double,ctypes.c_double,ctypes.c_int,ctypes.c_int,ctypes.c_double*3]
lib.fermion_pair_force.restype=None
lib.bohr_orbit_residual.argtypes=[ctypes.POINTER(Fermion),ctypes.c_double,ctypes.c_int]
lib.bohr_orbit_residual.restype=ctypes.c_double

def make_monopole(center,a=1.0,omega=0.3,phase0=0.0,qm=1.0,sign=1,
                  mu=1.0,nu=0.25,gamma=0.1):
    m=Monopole(); m.center[:]=center; m.e1[:]=(1.0,0.0,0.0); m.e2[:]=(0.0,1.0,0.0)
    m.a=a; m.omega=omega; m.phase0=phase0; m.qm=qm; m.sign=sign
    m.mu=mu; m.nu=nu; m.gamma=gamma; return m

def make_fermion(center,a=1.0,omega=0.3,phase0=0.0,qm=1.0,**kw):
    F=Fermion()
    F.plus =make_monopole(center,a,omega,phase0,qm,sign=+1,**kw)
    F.minus=make_monopole(center,a,omega,phase0,qm,sign=-1,**kw)
    return F

def self_energy(F,a):
    ea=ctypes.c_double(0); ep=ctypes.c_double(0)
    E=lib.fermion_self_energy(ctypes.byref(F),a,ctypes.byref(ea),ctypes.byref(ep))
    return E,ea.value,ep.value

def equilibrium_radius(F,amin=0.05,amax=20.0,n=20000):
    st=ctypes.c_int(0)
    a=lib.fermion_equilibrium_radius(ctypes.byref(F),amin,amax,n,ctypes.byref(st))
    return a,bool(st.value)

def pair_force(F1,F2,t,c,task,emode=ECHARGE_DARWIN):
    out=(ctypes.c_double*3)()
    lib.fermion_pair_force(ctypes.byref(F1),ctypes.byref(F2),t,c,task,emode,out)
    return [out[0],out[1],out[2]]

def bohr_residual(F,c,n):
    return lib.bohr_orbit_residual(ctypes.byref(F),c,n)
