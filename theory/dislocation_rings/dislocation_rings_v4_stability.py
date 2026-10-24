# dislocation_rings_v4_stability.py
# Добавляет к v4 МЕХАНИЗМ УСТОЙЧИВОСТИ: линейное натяжение дислокации.
# Полная энергия кольца E(a) = E_self_coulomb(a) + T·2πa, T=μb² (натяжение).
# Даёт устойчивый равновесный радиус a_eq ≈ √(kq²/2πT), E''>0.
# Это недостающий "отталкивающий механизм", которого нет в ядре v4 (.c).
#
# Использование: импортировать вместе с dislocation_rings_bind_v4.
import numpy as np

def ring_total_energy(a, q=1.0, r0=1e-3, k=1.0, T=0.05):
    """Полная энергия заряженного дислокационного кольца радиуса a.
    E_coul: кулоновская самоэнергия тонкого кольца (core r0), ~ (kq²/πa)·ln(8a/r0).
    E_tens: энергия линейного натяжения дислокации, T·2πa.
    Возвращает (E_total, E_coul, E_tens)."""
    E_coul = (k*q*q/(np.pi*a))*np.log(8.0*a/r0)
    E_tens = T*2.0*np.pi*a
    return E_coul+E_tens, E_coul, E_tens

def equilibrium_radius(q=1.0, r0=1e-3, k=1.0, T=0.05, amax=50.0, n=20000):
    """Численный поиск равновесного радиуса (минимум полной энергии)."""
    ag=np.linspace(1e-3, amax, n)
    E=np.array([ring_total_energy(a,q,r0,k,T)[0] for a in ag])
    i=int(np.argmin(E))
    stable = (i>0 and i<len(ag)-1 and (E[i+1]-2*E[i]+E[i-1])>0)
    return ag[i], E[i], stable

def analytic_a_eq(q=1.0, k=1.0, T=0.05):
    """Аналитическая оценка a_eq = √(kq²/(2πT)) (без лог-поправки)."""
    return np.sqrt(k*q*q/(2.0*np.pi*T))

if __name__=="__main__":
    a_eq,E_eq,stab = equilibrium_radius()
    print(f"a_eq(числ.)={a_eq:.4f}  E={E_eq:.4f}  устойчиво={stab}")
    print(f"a_eq(аналит.)≈{analytic_a_eq():.4f}")
    for T in [0.02,0.05,0.1,0.2,0.5]:
        a_eq,_,_=equilibrium_radius(T=T)
        print(f"  T={T:4.2f}: a_eq={a_eq:.3f}  a_eq·√T={a_eq*np.sqrt(T):.3f}")
