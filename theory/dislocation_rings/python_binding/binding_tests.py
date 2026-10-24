# -*- coding: utf-8 -*-
# binding_tests.py
# Прогон тестов через НАСТОЯЩИЙ C-код (ctypes). Требует предварительной сборки:
#   ./build.sh    (или вручную gcc ... -o *.so)
# Запуск: python3 binding_tests.py
# Совместимо с Python 3.5+ (без f-строк, через .format()).
# PYTHONIOENCODING=utf-8 python3 binding_tests.py
import math

def test_v4_stability():
    print("=== v4: устойчивость (натяжение) через C ===")
    import binding_v4 as B
    for T in [0.02,0.05,0.1,0.2,0.5]:
        a_eq,st=B.equilibrium_radius(T=T)
        print("  T={:.2f}: a_eq={:.3f} a_eq*sqrt(T)={:.3f} stable={}".format(
            T, a_eq, a_eq*math.sqrt(T), st))

def test_v4_dipole():
    print("=== v4: диполь-диполь F*R^2 (PLAIN/RET/INST) через C ===")
    import binding_v4 as B
    c=1.0
    for R in [5,10,20,40,80,160]:
        A=B.make_dipole((0,0,0)); C=B.make_dipole((R,0,0))
        v=[math.sqrt(sum(x*x for x in B.force_multi(A,C,0.0,c,m)))*R*R
           for m in (B.FERMI_PLAIN,B.FERMI_FILM_RETARDED,B.FERMI_FILM_INSTANT)]
        print("  R={:3d}: {:.5f} {:.5f} {:.5f}".format(R, v[0], v[1], v[2]))

def test_v4_mobius():
    print("=== v4: мёбиус-кольцо, вектор Бюргерса (период 4pi) через C ===")
    import binding_v4 as B
    R=B.make_ring_b((0,0,0),1.0,0.0,1.0,n_elem=1,sign=1,mobius=1)
    for frac in [0,0.5,1.0,1.5,2.0]:
        R.phase0=frac*2*math.pi
        bb=B.burgers_of(R,0,0.0)
        print("  phase={:.1f}*2pi: bx={:+.4f}".format(frac, bb[0]))

def test_v3_climb():
    print("=== v3: осевая сила N/S кольца (FERMI_NONE) через C ===")
    import binding_v3 as B3
    L=B3.v3; c=1.0; a=1.0
    for dfrac in [0.1,0.3,0.5,1.0,1.5,2.0,3.0]:
        d=dfrac*a
        N=L.make_ring_b((0,0,d/2),a,1.0,sign=1)
        S=L.make_ring_b((0,0,-d/2),a,1.0,sign=-1)
        F=L.force([N],[S],0.0,c,B3.FERMI_NONE)
        print("  d/a={:.2f}: Fz={:.6f}".format(dfrac, F[2]))

if __name__=="__main__":
    for t in (test_v3_climb,test_v4_stability,test_v4_dipole,test_v4_mobius):
        try:
            t()
        except OSError as e:
            print("[{}] .so не найден -- сначала ./build.sh ({})".format(t.__name__, e))
        except Exception as e:
            print("[{}] ошибка: {}".format(t.__name__, e))
        print("")
