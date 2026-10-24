# -*- coding: utf-8 -*-
# positronium_v5_tests.py — проверка пяти задач. Запуск: python3 positronium_v5_tests.py
# Сначала: gcc -O3 -shared -fPIC positronium_v5.c -o positronium_v5.so -lm
import math
import positronium_v5_bind as P

def t1_stability():
    print("=== ЗАДАЧА 1: устойчивость фермиона от самоаннигиляции ===")
    F=P.make_fermion((0,0,0),omega=0.3,qm=1.0,gamma=0.1)
    a_eq,stable=P.equilibrium_radius(F)
    E,ea,ep=P.self_energy(F,a_eq)
    print("  a_eq={:.4f}  E={:.4f} (притяж={:.4f}, прецессия={:.4f})  устойчиво={}".format(
        a_eq,E,ea,ep,stable))
    print("  -> прецессия создаёт барьер против стягивания пары +q_m,-q_m" if stable
          else "  -> равновесие не найдено")

def t2_coulomb():
    print("=== ЗАДАЧА 2: кулон двух фермионов (F*R^2, знак от Δφ) ===")
    c=1.0
    for R in [10,20,40,80]:
        A=P.make_fermion((0,0,0),omega=0.2,phase0=0.0)
        B=P.make_fermion((R,0,0),omega=0.2,phase0=0.0)
        F=P.pair_force(A,B,0.0,c,P.TASK_COULOMB)
        mag=math.sqrt(sum(x*x for x in F))
        print("  R={:3d}: |F|*R^2={:.5f}".format(R,mag*R*R))

def t3_gravity():
    print("=== ЗАДАЧА 3: гравитация (всегда притяжение, F*R^2) ===")
    c=1.0
    for R in [10,20,40,80]:
        A=P.make_fermion((0,0,0),omega=0.2)
        B=P.make_fermion((R,0,0),omega=0.2)
        F=P.pair_force(A,B,0.0,c,P.TASK_GRAVITY)
        mag=math.sqrt(sum(x*x for x in F))
        print("  R={:3d}: |F|*R^2={:.3e}".format(R,mag*R*R))

def t4_exchange():
    print("=== ЗАДАЧА 4: обменное взаимодействие (cos2Δφ) ===")
    c=1.0; R=20
    for dphi_deg in [0,45,90,135,180]:
        dphi=math.radians(dphi_deg)
        A=P.make_fermion((0,0,0),omega=0.2,phase0=0.0)
        B=P.make_fermion((R,0,0),omega=0.2,phase0=dphi)
        F=P.pair_force(A,B,0.0,c,P.TASK_EXCHANGE)
        mag=math.sqrt(sum(x*x for x in F))
        print("  Δφ={:3d}°: |F|={:.5e}  (cos2Δφ={:+.3f})".format(dphi_deg,mag,math.cos(2*dphi)))

def t5_bohr():
    print("=== ЗАДАЧА 5: боровская орбита монополей (резонанс nλ=2πa) ===")
    c=1.0
    # подберём omega так, чтобы для n=1 резидуал ~0: 2πa = λ = 2πc/ω => ω=c/a
    F=P.make_fermion((0,0,0),a=1.0,omega=1.0)  # ω=c/a=1
    for n in [1,2,3]:
        r=P.bohr_residual(F,c,n)
        print("  n={}: резидуал(2πa - nλ)={:+.4f} {}".format(n,r,"<- резонанс" if abs(r)<1e-6 else ""))

if __name__=="__main__":
    for t in (t1_stability,t2_coulomb,t3_gravity,t4_exchange,t5_bohr):
        try: t()
        except OSError as e: print("[{}] .so не найден — соберите gcc ({})".format(t.__name__,e))
        except Exception as e: print("[{}] ошибка: {}".format(t.__name__,e))
        print("")
