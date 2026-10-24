# -*- coding: utf-8 -*-
# unified_channels_bind.py (Python 3.5+)
# Сборка: gcc -O3 -shared -fPIC unified_channels.c -o unified_channels.so -lm
import ctypes, os
_h=os.path.dirname(os.path.abspath(__file__))
lib=ctypes.CDLL(os.path.join(_h,"unified_channels.so"))

PARTICLE_GLIDE_DIPOLE,PARTICLE_MOBIUS_DIPOLE,PARTICLE_SINGLE_RING,PARTICLE_QUARK_PARTIAL=0,1,2,3
CH_COULOMB,CH_MAGNETIC,CH_GRAVITY,CH_EXCHANGE,CH_WEAK,CH_STRONG=0,1,2,3,4,5

PARTICLE_NAMES={0:"glide-диполь(бозон)",1:"мёбиус-диполь(электрон)",2:"одиночное кольцо",3:"кварк(b/3)"}
CHANNEL_NAMES={0:"кулон",1:"магнитное",2:"гравитация",3:"обмен",4:"слабое",5:"сильное"}

class Particle(ctypes.Structure):
    _fields_=[("center",ctypes.c_double*3),("a",ctypes.c_double),("d",ctypes.c_double),
              ("omega",ctypes.c_double),("phase0",ctypes.c_double),("b0",ctypes.c_double),
              ("chirality",ctypes.c_int),("mu_shear",ctypes.c_double),("R_weak",ctypes.c_double),
              ("sigma_string",ctypes.c_double),("type",ctypes.c_int),("n_elem",ctypes.c_int)]

D3=ctypes.c_double*3
lib.make_particle.argtypes=[ctypes.c_int,D3,ctypes.c_double,ctypes.c_double,ctypes.c_int]
lib.make_particle.restype=Particle
lib.particle_charge.argtypes=[ctypes.POINTER(Particle)]; lib.particle_charge.restype=ctypes.c_double
lib.particle_magnetic_dipole.argtypes=[ctypes.POINTER(Particle)]; lib.particle_magnetic_dipole.restype=ctypes.c_double
lib.particle_mass.argtypes=[ctypes.POINTER(Particle),ctypes.c_double]; lib.particle_mass.restype=ctypes.c_double
lib.particle_period.argtypes=[ctypes.POINTER(Particle)]; lib.particle_period.restype=ctypes.c_double
lib.particle_spin.argtypes=[ctypes.POINTER(Particle)]; lib.particle_spin.restype=ctypes.c_double
lib.channel_force.argtypes=[ctypes.POINTER(Particle),ctypes.POINTER(Particle),
    ctypes.c_double,ctypes.c_double,ctypes.c_int,ctypes.c_int]; lib.channel_force.restype=ctypes.c_double
lib.channel_coupling.argtypes=[ctypes.POINTER(Particle),ctypes.POINTER(Particle),
    ctypes.c_double,ctypes.c_int]; lib.channel_coupling.restype=ctypes.c_double

def make(ptype,center=(0,0,0),a=1.0,omega=0.3,chirality=1):
    return lib.make_particle(ptype,D3(*center),a,omega,chirality)
def charge(p): return lib.particle_charge(ctypes.byref(p))
def magnetic_dipole(p): return lib.particle_magnetic_dipole(ctypes.byref(p))
def mass(p,c=1.0): return lib.particle_mass(ctypes.byref(p),c)
def period(p): return lib.particle_period(ctypes.byref(p))
def spin(p): return lib.particle_spin(ctypes.byref(p))
def force(p1,p2,R,c=1.0,ch=0,nT=24): return lib.channel_force(ctypes.byref(p1),ctypes.byref(p2),R,c,ch,nT)
def coupling(p1,p2,c=1.0,ch=0): return lib.channel_coupling(ctypes.byref(p1),ctypes.byref(p2),c,ch)
