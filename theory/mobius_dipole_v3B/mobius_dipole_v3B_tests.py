# -*- coding: utf-8 -*-
# gcc -O3 -shared -fPIC mobius_dipole_v3B.c -o mobius_dipole_v3B.so -lm
# PYTHONIOENCODING=utf-8 python3 mobius_dipole_v3B_tests.py
import math
import mobius_dipole_v3B_bind as V

def mk_pair(omega=1.0,d=1.0,a=1.0,db=1.0):
    N=V.make((0,0,+d/2),a=a,omega=omega,phase0=0.0,     delta_b=db,chirality=+1) # мёбиус N
    S=V.make((0,0,-d/2),a=a,omega=omega,phase0=math.pi, delta_b=db,chirality=-1) # анти S
    return N,S

def t_period():
    print("=== Период пары (ожид 4π = ФЕРМИОН) ===")
    N,S=mk_pair()
    T=V.pair_period(N,S)
    print("  период пары = {:.2f}π → {}".format(T/math.pi,
        "ФЕРМИОН (спин-1/2)" if abs(T-4*math.pi)<1e-9 else "бозон"))

def t_twist():
    print("=== Твист знака Бюргерса (мёбиус, период 4π) ===")
    N,S=mk_pair()
    for deg in [0,180,360,540,720]:
        t=math.radians(deg)
        print("  θ={:3d}°: sign_N={:+.3f} sign_S={:+.3f}".format(deg,V.burgers_sign(N,t),V.burgers_sign(S,t)))

def t_dipole():
    print("=== Магнитный диполь пары (q_m·d) ===")
    N,S=mk_pair()
    print("  μ_dipole={:+.3f} (ожид q_m·d=1.0)".format(V.pair_magnetic_dipole(N,S)))

def t_stability():
    print("=== Защита прецессией: d_eq=J²/(mA), ω=0→аннигиляция ===")
    for J in [0.0,0.5,1.0]:
        ds=[0.02*i for i in range(2,500)]
        Es=[(d,V.pair_energy(d,J,1.0,1.0)) for d in ds]
        dmin,Emin=min(Es,key=lambda x:x[1])
        state="АННИГИЛЯЦИЯ" if J<1e-9 else "стабильно d_eq={:.2f}".format(J*J)
        print("  J={:.1f}: {}".format(J,state if J<1e-9 else "d_eq={:.3f} (ожид {:.3f})".format(dmin,J*J)))

def t_annih_barrier():
    print("=== Двухступенчатая аннигиляция: барьер разворота твиста ===")
    N,S=mk_pair()
    print("  барьер ∫(db/dθ)²≈{:.3f} (>0 → двухступенчатая, в отличие от glide)".format(
        V.twist_unwind_barrier(N)))

if __name__=="__main__":
    for t in (t_period,t_twist,t_dipole,t_stability,t_annih_barrier):
        try: t()
        except OSError as e: print("[{}] .so нет — соберите gcc ({})".format(t.__name__,e))
        except Exception as e: print("[{}] ошибка: {}".format(t.__name__,e))
        print("")
