# Проверки к разбору гл. 6-7 Čápek & Sheehan (2005). Запуск: python ch6_7_checks.py
import numpy as np
R=8.314; kB=1.380649e-23; amu=1.6605e-27; h=6.62607e-34

# ---------- Гл.7: DSPG, таблицы 7.1-7.3, ур. 7.15-7.21, условие Вегшайдера ----------
T=1000.0; nu0=1e13; nuvib=1e13; rA=5e-10; vA=1.1e3; vA2=790.; EAA=240.
tab={1:dict(EA=250,EA2=260,Ed=0, FA=1e-2,FA2=1e3,gd=1e-1,gr=1e-6),   # Табл.7.2, S1
     2:dict(EA=200,EA2=190,Ed=30,FA=1e3, FA2=1.0,gd=1e-9,gr=1e-1)}   # Табл.7.2, S2
def surf(p, gr_scale=1.0):
    tA =p['FA'] /nu0*np.exp(p['EA'] *1e3/(R*T))                       # (7.20)
    tA2=p['FA2']/nu0*np.exp(p['EA2']*1e3/(R*T))
    K  =rA*vA/nuvib*p['gr']*gr_scale/p['gd']*np.exp(p['Ed']*1e3/(R*T)) # (7.21)
    HA,HA2=vA*tA/4, vA2*tA2/4            # n_s/n_g из (7.15),(7.16) (прилипание =1)
    Kg=K*HA**2/HA2                       # K(g), равновесный с поверхностью i (цикл газ-поверхность)
    return dict(tA=tA,tA2=tA2,K=K,HA=HA,HA2=HA2,Kg=Kg,hess=EAA-2*p['EA']+p['EA2'])
s1,s2=surf(tab[1]),surf(tab[2])
for i,s in ((1,s1),(2,s2)):
    print(f"S{i}: tau_A={s['tA']:.3g}s tau_A2={s['tA2']:.3g}s K={s['K']:.3g}m^2 "
          f"Kg={s['Kg']:.3g}m^3  Hess dU={s['hess']} kJ/mol")
print(f"W = Kg(2)/Kg(1) = {s2['Kg']/s1['Kg']:.3g}")
# (7.6) для сравнения
for i,s in ((1,s1),(2,s2)): print(f"S{i}: K(g) по (7.6) = {vA2/(2*nu0)*s['K']:.3g} m^3 против {s['Kg']:.3g} из разд.7.3")
# Табл.7.3: газ задаётся S2, S1 - примесь
nc2=2e16; n22=nc2*vA2*s2['tA2']/4; n2A=np.sqrt(n22/s2['K']); ncA=4*n2A/(vA*s2['tA'])
RT=vA/4*ncA+2*vA2/4*nc2
def s1_ss(s):
    a,b=2*s['K']/s['tA2'],1/s['tA']; x=(-b+np.sqrt(b*b+4*a*RT))/(2*a); y=s['K']*x*x
    return x,y,(y/s['tA2'])/(x/s['tA']), vA2/4*nc2-y/s['tA2']
x,y,a1,J=s1_ss(s1); a2=(n22/s2['tA2'])/(n2A/s2['tA'])
mA=20*amu
dP=lambda a1:(2-np.sqrt(2))*mA*vA*RT*(a2-a1)/((2*a1+1)*(2*a2+1))   # (7.26)
print(f"n(2,A2)={n22:.3g} n(2,A)={n2A:.3g} n(c,A)={ncA:.3g} n(1,A)={x:.3g} n(1,A2)={y:.3g}")
print(f"alpha1={a1:.3g} alpha2={a2:.3g} dP={dP(a1):.3g} Pa  J_cycle(S1)={J:.3g} m^-2 s^-1 "
      f"-> {J*EAA*1e3/6.022e23:.3g} W/m^2")
sW=surf(tab[1],gr_scale=s1['Kg']**-1*s2['Kg'])        # формально восстановить W=1
x,y,a1W,JW=s1_ss(sW); print(f"при W=1: alpha1={a1W:.3g} (=alpha2), dP={dP(a1W):.2g} Pa, J={JW:.2g}")
# статмех K_eq для гипотетического A2 (r=2e-10 м, nu_vib=1e13 Гц) - для ориентира
lam=lambda m: h/np.sqrt(2*np.pi*m*kB*T); I=10*amu*(2e-10)**2
Keq=lam(20*amu)**6/lam(40*amu)**3*(8*np.pi**2*I*kB*T/(2*h**2))/(1-np.exp(-h*nuvib/(kB*T)))*np.exp(EAA*1e3/(R*T))
print(f"K_eq(gas) ~ {Keq:.2g} m^3; Kg(1)/Keq={s1['Kg']/Keq:.1g}, Kg(2)/Keq={s2['Kg']/Keq:.1g}")

# ---------- Гл.6: 1D-модель гравитатора (нормальная компонента, m=kT=1) ----------
v=np.linspace(1e-6,12,240001); dv=v[1]-v[0]; jM=v*np.exp(-v**2/2)
def ICJM(ms,mA=4,Es=2,W=4):                      # (6.1), энергии в kT; пол 1e-3 для 1D связанных орбит
    mu=mA/ms; E=v**2/2; return np.clip(1-(1+mu)**2*E/(4*mu*(E-Es+W)),1e-3,1)
def wall_ratio(P,b,kernel):
    vinc=np.sqrt(2*b); ge=v>=vinc
    J=np.where(ge,v*np.exp(-(v**2-2*b)/2),0.)     # поток со стенки у поверхности
    D=jM if kernel=='sheehan' else P*jM; D=D/np.sum(D*dv)
    Tt=np.sum((P*J)[ge])/np.sum(D[ge])           # стационарный поток захваченных
    O=(1-P)*J+D*Tt; u=np.sqrt(np.clip(v**2-2*b,0,None))
    return np.sum((u*O)[ge])*dv/(np.sum(v*jM)*dv) # поток импульса вверх/вниз на уровне стенки
for b in (0.,1.,1.5):
    for k in ('sheehan','reciprocal'):
        r1,r2=wall_ratio(ICJM(7),b,k),wall_ratio(ICJM(96),b,k)
        print(f"PE={b}kT {k:10s}: S1 {r1:.4f}  S2 {r2:.4f}  P_S1/P_S2(стенка)={(1+r1)/(1+r2):.4f}")
for ms in (7,96):
    P=ICJM(ms); E=v**2/2; D=jM/np.sum(jM*dv)
    print(f"g=0, ядро Sheehan, ms={ms}: поток энергии out/in = "
          f"{(np.sum((1-P)*jM*E)+np.sum(P*jM)*np.sum(E*D)*dv)/np.sum(E*jM):.4f}")
A,B=0.5,2.0   # A*B=1 (PE=kT)
f63=np.sqrt(A/(v**2+B))*v*np.exp(-A*(v**2+B)); fe=np.exp(-A*v**2)
print("<v> по (6.3):",np.sum(v*f63)/np.sum(f63)," <v> полумаксвелл:",np.sum(v*fe)/np.sum(fe))
for gh in (0,1,3):   # Лошмидт: фильтрация по энергии не меняет T
    sel=v>=np.sqrt(2*gh); uu=np.sqrt(v[sel]**2-2*gh); w=jM[sel]
    print(f"mgh={gh}kT: <u^2>_поток={np.sum(w*uu**2)/np.sum(w):.4f} (2 = та же T), доля={np.sum(w)*dv:.4f}")
