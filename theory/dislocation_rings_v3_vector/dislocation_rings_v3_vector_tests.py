# -*- coding: utf-8 -*-
# gcc -O3 -shared -fPIC dislocation_rings_v3_vector.c -o dislocation_rings_v3_vector.so -lm
# PYTHONIOENCODING=utf-8 python3 dislocation_rings_v3_vector_tests.py
import math
import dislocation_rings_v3_vector_bind as V

def mk_pair(omega=1.0,d=1.0,a=1.0,db=1.0):
    top=V.make((0,0,+d/2),a=a,omega=omega,phase0=0.0,     b_vec=(1,0,0),delta_b=db)
    bot=V.make((0,0,-d/2),a=a,omega=omega,phase0=math.pi, b_vec=(1,0,0),delta_b=db) # противофаза
    return top,bot

def t_dipole():
    print("=== Дипольный момент: осевой p_z=q_m*d const, поперечный прецессирует ===")
    top,bot=mk_pair()
    for deg in [0,90,180,270]:
        t=math.radians(deg)
        p=V.pair_dipole(top,bot,t)
        print("  wt={:3d}: p=({:+.2f},{:+.2f},{:+.2f}) p_z={:+.2f} p_xy={:.2f}".format(
            deg,p[0],p[1],p[2],p[2],math.hypot(p[0],p[1])))

def t_moments():
    print("=== mu(ток), J, дипольный момент ===")
    top,bot=mk_pair()
    mu,J=V.pair_moments(top,bot,1.0)
    pdip=V.pair_magnetic_dipole(top,bot)
    print("  mu_z(ток)={:+.3f} (ожид ~0: токи сокращаются при противофазе)".format(mu[2]))
    print("  J_z={:+.3f} (механический, складывается)".format(J[2]))
    print("  магнитный ДИПОЛЬ q_m*d={:+.3f} (физический момент пары)".format(pdip))

def t_stability():
    print("=== Стабильность: центробежный барьер J^2/d^2, d_eq=J^2/(mA) ===")
    for Jc in [0.5,1.0,2.0]:
        ds=[0.05*i for i in range(2,400)]
        Es=[(d,V.pair_energy(d,Jc,1.0,1.0)) for d in ds]  # m=1,A=1
        dmin,Emin=min(Es,key=lambda x:x[1])
        print("  J={:.1f}: d_eq={:.3f} (ожид J^2/(mA)={:.3f})  E_min={:.4f}".format(Jc,dmin,Jc*Jc,Emin))

def t_annihilation():
    print("=== Аннигиляция: суммарный Бюргерс пары (=0) ===")
    top,bot=mk_pair()
    tb=V.total_burgers(top,bot)
    print("  b(e-)+b(e+)=({:+.1f},{:+.1f},{:+.1f}) |сумма|={:.3f} -> {}".format(
        tb[0],tb[1],tb[2],math.sqrt(sum(x*x for x in tb)),
        "ПОЛНАЯ аннигиляция" if abs(sum(tb))<1e-9 else "остаток"))

def t_mass():
    print("=== Масса из самодействия m~mu_shear*q_m^2*a/c^2 ===")
    top,_=mk_pair()
    for c in [1.0,2.0]:
        print("  c={:.1f}: m={:.3f}".format(c,V.loop_mass(top,c,1.0)))

if __name__=="__main__":
    for t in (t_dipole,t_moments,t_stability,t_annihilation,t_mass):
        try: t()
        except OSError as e: print("[{}] .so нет — соберите gcc ({})".format(t.__name__,e))
        except Exception as e: print("[{}] ошибка: {}".format(t.__name__,e))
        print("")
