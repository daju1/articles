# -*- coding: utf-8 -*-
# unified_channels_tests.py
# gcc -O3 -shared -fPIC unified_channels.c -o unified_channels.so -lm
# PYTHONIOENCODING=utf-8 python3 unified_channels_tests.py
import math
import unified_channels_bind as U

def t_particles():
    print("="*72)
    print("СВОЙСТВА ЧАСТИЦ (заряд, диполь, масса, период, спин)")
    print("="*72)
    print("{:26} {:>7} {:>8} {:>8} {:>8} {:>6}".format("тип","q_m","μ=q·d","масса","период/π","спин"))
    for pt in [U.PARTICLE_GLIDE_DIPOLE,U.PARTICLE_MOBIUS_DIPOLE,U.PARTICLE_SINGLE_RING,U.PARTICLE_QUARK_PARTIAL]:
        p=U.make(pt,a=1.0,omega=0.3)
        print("{:26} {:7.3f} {:8.3f} {:8.4f} {:8.2f} {:6.1f}".format(
            U.PARTICLE_NAMES[pt],U.charge(p),U.magnetic_dipole(p),U.mass(p,1.0),
            U.period(p)/math.pi,U.spin(p)))

def t_channels_all():
    print("\n"+"="*72)
    print("ВСЕ КАНАЛЫ × ВСЕ ПАРЫ ЧАСТИЦ: |F| на R=10 (усреднено по периоду)")
    print("="*72)
    c=1.0; R=10.0
    types=[U.PARTICLE_GLIDE_DIPOLE,U.PARTICLE_MOBIUS_DIPOLE,U.PARTICLE_QUARK_PARTIAL]
    for pt in types:
        p1=U.make(pt,center=(0,0,0),a=1.0,omega=0.3)
        p2=U.make(pt,center=(R,0,0),a=1.0,omega=0.3)
        print("\nпара: {} + {}".format(U.PARTICLE_NAMES[pt],U.PARTICLE_NAMES[pt]))
        for ch in range(6):
            F=U.force(p1,p2,R,c,ch,24)
            print("   {:12}: |F|={:.4e}".format(U.CHANNEL_NAMES[ch],F))

def t_scaling():
    print("\n"+"="*72)
    print("ЗАКОНЫ СПАДАНИЯ: показатель p в F~R^-p (мёбиус-электрон)")
    print("="*72)
    c=1.0; Rs=[10.,20.,40.,80.]
    for ch in range(6):
        mags=[]
        for R in Rs:
            p1=U.make(U.PARTICLE_MOBIUS_DIPOLE,center=(0,0,0),a=1.0,omega=0.3)
            p2=U.make(U.PARTICLE_MOBIUS_DIPOLE,center=(R,0,0),a=1.0,omega=0.3)
            mags.append(U.force(p1,p2,R,c,ch,24))
        if mags[0]>0 and mags[-1]>0:
            p=-math.log(mags[-1]/mags[0])/math.log(Rs[-1]/Rs[0])
            law="~R^-{:.2f}".format(p)
        else:
            law="(нет/короткодейств.)"
        print("   {:12}: {}".format(U.CHANNEL_NAMES[ch],law))

def t_weak_parity():
    print("\n"+"="*72)
    print("СЛАБОЕ: нарушение чётности (одинаковая vs разная киральность)")
    print("="*72)
    c=1.0; R=3.0
    for ch1,ch2,lab in [(1,1,"одинак (LL)"),(1,-1,"разные (LR)")]:
        p1=U.make(U.PARTICLE_MOBIUS_DIPOLE,center=(0,0,0),chirality=ch1)
        p2=U.make(U.PARTICLE_MOBIUS_DIPOLE,center=(R,0,0),chirality=ch2)
        F=U.force(p1,p2,R,c,U.CH_WEAK,24)
        print("   киральности {}: |F_weak|={:.4e}".format(lab,F))

if __name__=="__main__":
    for t in (t_particles,t_channels_all,t_scaling,t_weak_parity):
        try: t()
        except OSError as e: print("[{}] .so не найден — соберите gcc ({})".format(t.__name__,e))
        except Exception as e: print("[{}] ошибка: {}".format(t.__name__,e))
