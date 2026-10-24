# -*- coding: utf-8 -*-
# unified_channels_py.py — НЕЗАВИСИМАЯ Python-реализация тех же формул
# (зеркало C-кода для перепроверки). Не требует .so. Py3.5+, без f-strings.
import math

GLIDE,MOBIUS,SINGLE,QUARK=0,1,2,3
CH_COULOMB,CH_MAGNETIC,CH_GRAVITY,CH_EXCHANGE,CH_WEAK,CH_STRONG=0,1,2,3,4,5
PNAMES={0:"glide-диполь(бозон)",1:"мёбиус-диполь(электрон)",2:"одиночное кольцо",3:"кварк(b/3)"}
CNAMES={0:"кулон",1:"магнитное",2:"гравитация",3:"обмен",4:"слабое",5:"сильное"}

def make(ptype,center=(0,0,0),a=1.0,omega=0.3,chirality=1):
    b0=1.0/3.0 if ptype==QUARK else 1.0
    return {"center":list(center),"a":a,"d":a,"omega":omega,"phase0":0.0,"b0":b0,
            "chirality":chirality,"mu_shear":1.0,"R_weak":0.05,"sigma_string":0.05,
            "type":ptype,"n_elem":48}

def charge(p): return p["b0"]
def magnetic_dipole(p): return p["b0"]*p["d"]
def mass(p,c=1.0): return p["mu_shear"]*p["b0"]**2*p["a"]/c**2
def period(p): return 4*math.pi if p["type"]==MOBIUS else 2*math.pi
def spin(p): return 0.5 if p["type"]==MOBIUS else 0.0

def _twist(p,theta):
    if p["type"]==MOBIUS: return math.cos(p["chirality"]*theta/2.0)
    return 1.0
def _accel(p): return p["omega"]**2*p["a"]

def _force_inst(p1,p2,R,c,ch,t):
    q1,q2=p1["b0"],p2["b0"]
    ph1=p1["phase0"]+p1["omega"]*t
    ph2=p2["phase0"]+p2["omega"]*t+math.pi
    dphi=ph1-ph2
    if ch==CH_COULOMB:
        tw=_twist(p1,ph1)*_twist(p2,ph2)
        return abs(q1*q2*tw)/R**2
    if ch==CH_MAGNETIC:
        return abs((q1*p1["d"])*(q2*p2["d"]))/R**4
    if ch==CH_GRAVITY:
        W1,W2=_accel(p1),_accel(p2)
        s1=p1["mu_shear"]*q1*q1*W1*W1/c**2
        s2=p2["mu_shear"]*q2*q2*W2*W2/c**2
        return s1*s2/(4*math.pi*p1["mu_shear"]*R**2)
    if ch==CH_EXCHANGE:
        return abs(q1*q1*q2*q2*math.cos(2*dphi))/R**3
    if ch==CH_WEAK:
        chir=p1["chirality"]*p2["chirality"]; proj=(1.0+chir)/2.0
        return proj*abs(q1*q2)*math.exp(-R/p1["R_weak"])/R
    if ch==CH_STRONG:
        frac=1.0 if (p1["b0"]<0.9 and p2["b0"]<0.9) else 0.0
        return frac*p1["sigma_string"]*R
    return 0.0

def force(p1,p2,R,c=1.0,ch=0,nT=24):
    T=2*math.pi/p1["omega"]; s=0.0
    for i in range(nT):
        s+=_force_inst(p1,p2,R,c,ch,T*i/nT)
    return s/nT

if __name__=="__main__":
    print("Python-зеркало: свойства частиц")
    for pt in [GLIDE,MOBIUS,SINGLE,QUARK]:
        p=make(pt)
        print("  {:26} q={:.3f} μ={:.3f} m={:.4f} период={:.2f}π спин={:.1f}".format(
            PNAMES[pt],charge(p),magnetic_dipole(p),mass(p),period(p)/math.pi,spin(p)))
    print("\nвсе каналы, мёбиус-пара, R=10:")
    a=make(MOBIUS,(0,0,0)); b=make(MOBIUS,(10,0,0))
    for ch in range(6):
        print("  {:12}: |F|={:.4e}".format(CNAMES[ch],force(a,b,10.0,1.0,ch)))
