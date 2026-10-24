# -*- coding: utf-8 -*-
# disclination_fermion_v6_tests.py
# gcc -O3 -shared -fPIC disclination_fermion_v6.c -o disclination_fermion_v6.so -lm
# python3 disclination_fermion_v6_tests.py
import math
import disclination_fermion_v6_bind as V

def t_multipoles():
    print("=== Твист ≡ пара зарядов: мультиполи ===")
    D=V.make((0,0,0))
    m,d=V.multipoles(D)
    print("  монополь={:.3e} (≈0), |диполь|={:.4f} (ожид 8/3={:.4f})".format(m,d,8/3))
    print("  -> нулевой монополь + ненулевой диполь = эффективная пара ±q_m")

def t_spin():
    print("=== Спин и g-фактор из твиста (период 4π) ===")
    D=V.make((0,0,0))
    print("  S/L_орб={:.3f} (ожид 1/2), g={:.1f} (ожид 2)".format(V.spin_ratio(D),V.g_factor(D)))

def t_channels():
    print("=== Три канала (УСРЕДНЕНО по периоду): закон силы и величина ===")
    c=1.0; Rs=[10,20,40,80]
    for ch,name in [(V.CH_COULOMB,"кулон"),(V.CH_EXCHANGE,"обмен"),(V.CH_GRAVITY,"гравитация")]:
        mags=[]
        for R in Rs:
            A=V.make((0,0,0),omega=0.2); B=V.make((R,0,0),omega=0.2)
            F=V.pair_force_avg(A,B,c,ch,nT=16); mag=math.sqrt(sum(x*x for x in F)); mags.append(mag)
        # оценка показателя F~R^-p по крайним точкам
        p=-math.log(mags[-1]/mags[0])/math.log(Rs[-1]/Rs[0])
        print("  {:10s}: |F|(R=10)={:.3e}  закон F~R^-{:.2f}".format(name,mags[0],p))

if __name__=="__main__":
    for t in (t_multipoles,t_spin,t_channels):
        try: t()
        except OSError as e: print("[{}] .so не найден — соберите gcc ({})".format(t.__name__,e))
        except Exception as e: print("[{}] ошибка: {}".format(t.__name__,e))
        print("")
