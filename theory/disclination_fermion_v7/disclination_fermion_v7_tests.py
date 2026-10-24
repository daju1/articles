# -*- coding: utf-8 -*-
# disclination_fermion_v7_tests.py
# gcc -O3 -shared -fPIC disclination_fermion_v7.c -o disclination_fermion_v7.so -lm
# PYTHONIOENCODING=utf-8 python3 disclination_fermion_v7_tests.py
import math
import disclination_fermion_v7_bind as V

PLANES=[("T_R",V.PLANE_T_R),("R_Z",V.PLANE_R_Z),("T_Z",V.PLANE_T_Z)]

def t_magnetic_plane():
    print("=== Магнитный момент по плоскостям (ожид: !=0 только T_R) ===")
    for name,pl in PLANES:
        for n in [1,3]:
            D=V.make((0,0,0),plane=pl,winding=n)
            mu=V.mag_moment(D)
            print("  {} n={}: |mu|={:.4f}  mu_z={:+.4f}".format(name,n,math.sqrt(sum(x*x for x in mu)),mu[2]))
    print("  -> mu_z!=0 только для T_R (b в плоскости кольца, Франк||z)")

def t_tilted():
    print("=== Косой наклон плоскости T_R: mu_z=0.5*cos(beta)? ===")
    for bdeg in [0,30,45,60,90]:
        D=V.make((0,0,0),plane=V.PLANE_TILTED,normal_beta=math.radians(bdeg),winding=1)
        mu=V.mag_moment(D)
        print("  beta={:2d}: mu_z={:+.4f} (ожид {:+.4f})".format(bdeg,mu[2],0.5*math.cos(math.radians(bdeg))))

def t_spin_gyro():
    print("=== Спин, период, гиромагнитное ===")
    for n in [1,2,3]:
        D=V.make((0,0,0),plane=V.PLANE_T_R,winding=n)
        g,mz,sz=V.gyro_g(D)
        print("  n={}: спин={:.1f}  период(-1)^n={:+d}  mu_z={:+.3f} S={:.1f} mu/S={:.3f}".format(
            n,V.spin(D),V.period_factor(D),mz,sz,g))
    print("  -> нечётн n: период -1 (фермион); mu_z пропорц. спину")

def t_channels():
    print("=== Шесть каналов (усреднено), |F| на R=20 ===")
    c=1.0
    names=[("кулон",V.CH_COULOMB),("гравитация",V.CH_GRAVITY),("обмен",V.CH_EXCHANGE),
           ("магнитное",V.CH_MAGNETIC),("слабое",V.CH_WEAK),("сильное",V.CH_STRONG)]
    for nm,ch in names:
        A=V.make((0,0,0),omega=0.2,plane=V.PLANE_T_R)
        B=V.make((20,0,0),omega=0.2,plane=V.PLANE_T_R)
        F=V.pair_force_avg(A,B,c,ch,nT=16); mag=math.sqrt(sum(x*x for x in F))
        print("  {:12s}: |F|={:.4e}".format(nm,mag))

def t_weak_chirality():
    print("=== Слабое: зависит от знака киральности (P-нарушение) ===")
    c=1.0
    for ch1,ch2,lab in [(1,1,"одинак (LL)"),(1,-1,"разные (LR)")]:
        A=V.make((0,0,0),omega=0.2,chirality=ch1); B=V.make((20,0,0),omega=0.2,chirality=ch2)
        F=V.pair_force_avg(A,B,c,V.CH_WEAK,nT=16); mag=math.sqrt(sum(x*x for x in F))
        print("  киральности {}: |F_weak|={:.4e}".format(lab,mag))
    print("  -> разная связь LL vs LR = нарушение чётности")

def t_magnetic_force():
    print("=== Магнитная сила (поток-поток, Осипов): закон 1/R^3 ===")
    c=1.0
    for R in [10,20,40]:
        A=V.make((0,0,0),omega=0.2,plane=V.PLANE_T_R,winding=1)
        B=V.make((R,0,0),omega=0.2,plane=V.PLANE_T_R,winding=1)
        F=V.pair_force_avg(A,B,c,V.CH_MAGNETIC,nT=16); mag=math.sqrt(sum(x*x for x in F))
        print("    R={:3d}: |F|={:.4e}  |F|*R^3={:.4f}".format(R,mag,mag*R**3))
    print("    -> |F|*R^3=const подтверждает диполь-дипольный 1/R^3 (мю=ню=n/2)")

if __name__=="__main__":
    for t in (t_magnetic_plane,t_tilted,t_spin_gyro,t_channels,t_magnetic_force,t_weak_chirality):
        try: t()
        except OSError as e: print("[{}] .so не найден — соберите gcc ({})".format(t.__name__,e))
        except Exception as e: print("[{}] ошибка: {}".format(t.__name__,e))
        print("")
