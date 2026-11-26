import numpy as np, time
pi=np.pi
trapz=getattr(np,'trapezoid',None) or np.trapz
def check(*a): pass
def K0K1(x):
    """K0(x), K1(x) для x>0 — Абрамовиц–Стиган 9.8.1-9.8.8."""
    x = np.asarray(x, float); K0 = np.empty_like(x); K1 = np.empty_like(x)
    s = x <= 2.0; t = x[s]/2; t2 = (x[s]/3.75)**2
    I0 = 1 + t2*(3.5156229 + t2*(3.0899424 + t2*(1.2067492 + t2*(0.2659732 + t2*(0.0360768 + t2*0.0045813)))))
    I1 = x[s]*(0.5 + t2*(0.87890594 + t2*(0.51498869 + t2*(0.15084934 + t2*(0.02658733 + t2*(0.00301532 + t2*0.00032411))))))
    y = t*t
    K0[s] = -np.log(t)*I0 + (-0.57721566 + y*(0.42278420 + y*(0.23069756 + y*(0.03488590 + y*(0.00262698 + y*(0.00010750 + y*0.0000074))))))
    K1[s] = np.log(t)*I1 + (1/x[s])*(1 + y*(0.15443144 + y*(-0.67278579 + y*(-0.18156897 + y*(-0.01919402 + y*(-0.00110404 + y*(-0.00004686)))))))
    b = ~s; xb = x[b]; z = 2/xb; e = np.exp(-xb)/np.sqrt(xb)
    K0[b] = e*(1.25331414 + z*(-0.07832358 + z*(0.02189568 + z*(-0.01062446 + z*(0.00587872 + z*(-0.00251540 + z*0.00053208))))))
    K1[b] = e*(1.25331414 + z*(0.23498619 + z*(-0.03655620 + z*(0.01504268 + z*(-0.00780353 + z*(0.00325614 + z*(-0.00068245)))))))
    return K0, K1



def k0_selfpanel(kap, h):
    """S_jj: (1/2pi) * интеграл K0(kap|t|) по панели длины h (особенность ln вычтена аналитически)."""
    xg, wg = np.polynomial.legendre.leggauss(16)
    T = h/2; t = (xg+1)/2*T; w = wg/2*T
    K0, _ = K0K1(kap*t)
    reg = K0 + np.log(kap*t/2)
    return (np.sum(w*reg) - (T*np.log(kap*T/2) - T)) / pi

def panels(P0, P1, N):
    """N панелей на отрезке P0-P1, сгущение по косинусу к концам. Возвращает (середины, длины)."""
    u = (1 - np.cos(np.linspace(0, pi, N+1))) / 2
    P0 = np.asarray(P0, float); P1 = np.asarray(P1, float)
    nodes = P0[None, :] + u[:, None]*(P1 - P0)[None, :]
    return (nodes[1:] + nodes[:-1])/2, np.linalg.norm(nodes[1:] - nodes[:-1], axis=1)

def rect_wing(P0, P1, t, N, Nt):
    """Крыло толщины t: прямоугольник вдоль P0->P1; N панелей на длинные стороны, Nt на торцы."""
    P0 = np.asarray(P0, float); P1 = np.asarray(P1, float)
    u = (P1-P0)/np.linalg.norm(P1-P0); n = np.array([-u[1], u[0]])
    A, B, C, D = P0-n*t/2, P1-n*t/2, P1+n*t/2, P0+n*t/2
    segs = [panels(A, B, N), panels(B, C, Nt), panels(C, D, N), panels(D, A, Nt)]
    return np.vstack([s[0] for s in segs]), np.concatenate([s[1] for s in segs])

def Smat(mid, h, kap):
    """Матрица простого слоя S_jb."""
    R = np.linalg.norm(mid[:, None, :] - mid[None, :, :], axis=2); np.fill_diagonal(R, 1.0)
    K0, _ = K0K1(kap*R); S = K0/(2*pi)*h[None, :]
    for j in range(len(h)): S[j, j] = k0_selfpanel(kap, h[j])
    return S

def G0grad(r, s, k):
    """g0(r_m - s_b) и его градиент по r: массивы (M,B) и (M,B,2)."""
    d = r[:, None, :] - s[None, :, :]; R = np.linalg.norm(d, axis=2)
    K0, K1 = K0K1(k*R)
    return K0/(2*pi), (-k*K1/(2*pi*R))[..., None]*d



def kgrid(n=60, s=1.0):
    """Узлы и веса для интеграла по kappa от 0 до бесконечности: kappa = s*tan(pi u/2)."""
    u = (np.arange(n) + 0.5)/n
    return s*np.tan(pi/2*u), s*pi/2/np.cos(pi/2*u)**2/n

def tensor(bodies, pts, nk=60, smin=None):
    """<Theta_ij>(r) в точках pts (M,2), интегрированный по kappa. Возвращает (M,2,2)."""
    mid = np.vstack([b[0] for b in bodies]); h = np.concatenate([b[1] for b in bodies])
    dmin = np.min(np.linalg.norm(pts[:, None, :] - mid[None, :, :], axis=2))
    kk, ww = kgrid(nk, smin if smin else 1.0/dmin)
    Th = np.zeros((len(pts), 2, 2))
    for k, w in zip(kk, ww):
        S = Smat(mid, h, k); G, gr = G0grad(pts, mid, k)     # g0(r-m_b), grad_r g0(r-m_b)
        X0 = np.linalg.solve(S, G.T)                          # S^{-1} g0(m-r)
        Xg = np.stack([np.linalg.solve(S, gr[:, :, j].T) for j in range(2)], -1)
        # grad_{r'} g0(m-r') = +grad_r g0(r-m)|_{r=r'} * (-1) * (-1): g0 зависит от |r-m|, d/dr' g0(m-r') = -d/dm = grad_r g0(r-m)
        gs = -np.einsum('mb,bm->m', G*h[None, :], X0)
        D = -np.einsum('mbi,bmj->mij', gr*h[None, :, None], Xg)
        Th += w*k/(2*pi)*(D - 0.5*np.eye(2)[None]*(k*k*gs + D[:, 0, 0] + D[:, 1, 1])[:, None, None])
    return Th

def traction(bodies, pts, nrm, nk=60, smin=None):
    """Плотность силы на единицу длины контура: sigma_ij n_j = -Theta_ij n_j."""
    return -np.einsum('mij,mj->mi', tensor(bodies, pts, nk, smin), nrm)

def graded(L, eps, nq=8, ratio=2.0):
    """Составная квадратура Гаусса на [0,L] с геометрическим сгущением к обоим концам."""
    edges = [0.0]; x = eps/2
    while x < L/2: edges.append(x); x *= ratio
    half = np.array(edges + [L/2]); allE = np.concatenate([half, L - half[::-1][1:]])
    xg, wg = np.polynomial.legendre.leggauss(nq); s = []; w = []
    for a_, b_ in zip(allE[:-1], allE[1:]):
        s.extend(a_ + (xg+1)/2*(b_-a_)); w.extend(wg/2*(b_-a_))
    return np.array(s), np.array(w)

def hug(P0, P1, t, eps, nq=8, narc=16):
    """Контур на расстоянии eps вокруг крыла P0->P1 толщины t: точки, нормали, веса dl, метки частей."""
    P0 = np.asarray(P0, float); P1 = np.asarray(P1, float); L = np.linalg.norm(P1-P0)
    u = (P1-P0)/L; n = np.array([-u[1], u[0]]); hw = t/2
    P = []; N_ = []; DL = []; LB = []
    def line(A, B, nn, lab):
        Ls = np.linalg.norm(B-A); s, w = graded(Ls, eps, nq)
        P.extend(A + np.outer(s/Ls, B-A)); N_.extend([nn]*len(s)); DL.extend(w); LB.extend([lab]*len(s))
    def arc(C, th0, th1, lab):
        xg, wg = np.polynomial.legendre.leggauss(narc); th = th0 + (xg+1)/2*(th1-th0)
        d = np.cos(th)[:, None]*u + np.sin(th)[:, None]*n
        P.extend(C + eps*d); N_.extend(d); DL.extend(eps*wg/2*(th1-th0)); LB.extend([lab]*narc)
    line(P0-n*(hw+eps), P1-n*(hw+eps), -n, 'side-')
    arc(P1-n*hw, -pi/2, 0, 'tip'); line(P1+u*eps-n*hw, P1+u*eps+n*hw, u, 'tip'); arc(P1+n*hw, 0, pi/2, 'tip')
    line(P1+n*(hw+eps), P0+n*(hw+eps), n, 'side+')
    arc(P0+n*hw, pi/2, pi, 'vertex'); line(P0-u*eps+n*hw, P0-u*eps-n*hw, -u, 'vertex'); arc(P0-n*hw, pi, 1.5*pi, 'vertex')
    return np.array(P), np.array(N_), np.array(DL), np.array(LB)



def poly_body(V, Ns):
    """Многоугольник с вершинами V (против часовой), Ns панелей на сторонах: (середины, длины, внешние нормали)."""
    M = []; H = []; NN = []
    for i in range(len(V)):
        X = np.asarray(V[i], float); Y = np.asarray(V[(i+1) % len(V)], float)
        m, h = panels(X, Y, Ns[i]); tau = (Y-X)/np.linalg.norm(Y-X)
        M.append(m); H.append(h); NN.append(np.tile([tau[1], -tau[0]], (Ns[i], 1)))
    return np.vstack(M), np.concatenate(H), np.vstack(NN)

def rect_wing_n(P0, P1, t, N, Nt):
    """То же крыло, что rect_wing, но с нормалями (обход против часовой)."""
    P0 = np.asarray(P0, float); P1 = np.asarray(P1, float)
    u = (P1-P0)/np.linalg.norm(P1-P0); n = np.array([-u[1], u[0]])
    return poly_body([P0-n*t/2, P1-n*t/2, P1+n*t/2, P0+n*t/2], [N, Nt, N, Nt])

def G0hess(r, s, k):
    d = r[:, None, :] - s[None, :, :]; R = np.linalg.norm(d, axis=2)
    K0, K1 = K0K1(k*R); e = d/R[..., None]
    g1 = -k*K1/(2*pi); g2 = k*k/(2*pi)*(K0 + K1/(k*R))
    return (g2 - g1/R)[..., None, None]*e[..., :, None]*e[..., None, :] + (g1/R)[..., None, None]*np.eye(2)

def Nmat(mid, h, nrm, k):
    """A = -I/2 + K',  K'_jb = h_b n_j . grad g0(m_j - m_b),  K'_jj = 0 (прямая панель)."""
    d = mid[:, None, :] - mid[None, :, :]; R = np.linalg.norm(d, axis=2); np.fill_diagonal(R, 1.0)
    _, K1 = K0K1(k*R); gr = (-k*K1/(2*pi*R))[..., None]*d
    Kp = np.einsum('jbc,jc->jb', gr, nrm)*h[None, :]; np.fill_diagonal(Kp, 0.0)
    return -0.5*np.eye(len(h)) + Kp

def tensorN(bodies, pts, nk=60, smin=None):
    """<Theta_ij> для условия Неймана (ТЕ)."""
    mid = np.vstack([b[0] for b in bodies]); h = np.concatenate([b[1] for b in bodies]); nrm = np.vstack([b[2] for b in bodies])
    dmin = np.min(np.linalg.norm(pts[:, None, :] - mid[None, :, :], axis=2))
    kk, ww = kgrid(nk, smin if smin else 1.0/dmin); Th = np.zeros((len(pts), 2, 2))
    for k, w in zip(kk, ww):
        A = Nmat(mid, h, nrm, k); G, gr = G0grad(pts, mid, k); Hs = G0hess(pts, mid, k)
        Bv = np.einsum('mbi,bi->bm', gr, nrm)          # b_c(r') = -dn_c g0(m_c - r') = n_c . grad_r g0(r - m_c)
        Bd = np.einsum('mbij,bi->bmj', Hs, nrm)        # d/dr'_j  b_c(r')
        X0 = np.linalg.solve(A, Bv); Xd = np.stack([np.linalg.solve(A, Bd[:, :, j]) for j in range(2)], -1)
        gs = np.einsum('mb,bm->m', G*h[None, :], X0)
        D = np.einsum('mbi,bmj->mij', gr*h[None, :, None], Xd)
        Th += w*k/(2*pi)*(D - 0.5*np.eye(2)[None]*(k*k*gs + D[:, 0, 0] + D[:, 1, 1])[:, None, None])
    return Th

def tractionN(bodies, pts, nrm, nk=60, smin=None):
    return -np.einsum('mij,mj->mi', tensorN(bodies, pts, nk, smin), nrm)

def energyN(bodies, nk=60, ks=3.0):
    mid = np.vstack([b[0] for b in bodies]); h = np.concatenate([b[1] for b in bodies]); nrm = np.vstack([b[2] for b in bodies])
    kk, ww = kgrid(nk, ks); E = 0.0
    for k, w in zip(kk, ww):
        E += w*k*(np.linalg.slogdet(Nmat(mid, h, nrm, k))[1] - sum(np.linalg.slogdet(Nmat(b[0], b[1], b[2], k))[1] for b in bodies))
    return E/(4*pi)



def Smat_l(mid, h, k, lam):
    S = Smat(mid, h, k)
    return S if lam is None else S + np.eye(len(h))/lam

def circle_force_l(bodies, center, Rc, lam, M=400, nk=60):
    mid = np.vstack([b[0] for b in bodies]); h = np.concatenate([b[1] for b in bodies])
    th = 2*pi*np.arange(M)/M; n = np.stack([np.cos(th), np.sin(th)], 1); pts = center + Rc*n
    dmin = np.min(np.linalg.norm(pts[:, None, :] - mid[None, :, :], axis=2)); kk, ww = kgrid(nk, 1/dmin); F = np.zeros(2)
    for k, w in zip(kk, ww):
        S = Smat_l(mid, h, k, lam); Gm, gr = G0grad(pts, mid, k)
        X0 = np.linalg.solve(S, Gm.T); Xg = np.stack([np.linalg.solve(S, gr[:, :, j].T) for j in range(2)], -1)
        gs = -np.einsum('mb,bm->m', Gm*h[None, :], X0); D = -np.einsum('mbi,bmj->mij', gr*h[None, :, None], Xg)
        Th = D - 0.5*np.eye(2)[None]*(k*k*gs + D[:, 0, 0] + D[:, 1, 1])[:, None, None]
        F += w*k/(2*pi)*np.einsum('mij,mj->i', Th, n)*(2*pi*Rc/M)
    return -F

def force_l(W1, W2, lam, nk=60, ks=3.0):
    (m1, h1), (m2, h2) = W1, W2; mid = np.vstack([m1, m2]); h = np.concatenate([h1, h2]); N1 = len(h1)
    def gl(mid, h, k):
        Si = np.linalg.inv(Smat_l(mid, h, k, lam))
        d = mid[:, None, :] - mid[None, :, :]; R = np.linalg.norm(d, axis=2); np.fill_diagonal(R, 1.0)
        _, K1 = K0K1(k*R); np.fill_diagonal(K1, 0.0); g = (k*K1/(2*pi*R))[..., None]*d
        return np.einsum('ji,ijc->jc', Si, g)*h[:, None] + np.einsum('ij,ijc,i->jc', Si, g, h)
    kk, ww = kgrid(nk, ks); F = np.zeros((len(h), 2))
    for k, w in zip(kk, ww):
        G_ = gl(mid, h, k); G_[:N1] -= gl(m1, h1, k); G_[N1:] -= gl(m2, h2, k); F += -(1/(4*pi))*w*k*G_
    return F[:N1].sum(0), F[N1:].sum(0)




def tensor_l(bodies, pts, lam, nk=50, smin=None, chunk=150):
    mid=np.vstack([b[0] for b in bodies]); h=np.concatenate([b[1] for b in bodies])
    dmin=np.min(np.linalg.norm(pts[:,None,:]-mid[None,:,:],axis=2)); kk,ww=kgrid(nk,smin if smin else 1/dmin)
    Th=np.zeros((len(pts),2,2))
    for k,w in zip(kk,ww):
        Si=np.linalg.inv(Smat_l(mid,h,k,lam))
        for i0 in range(0,len(pts),chunk):
            q=pts[i0:i0+chunk]; Gm,gr=G0grad(q,mid,k)
            X0=Si@Gm.T; Xg=np.einsum('bc,mcj->bmj',Si,gr)
            gs=-np.einsum('mb,bm->m',Gm*h[None,:],X0); D=-np.einsum('mbi,bmj->mij',gr*h[None,:,None],Xg)
            Th[i0:i0+chunk]+=w*k/(2*pi)*(D-0.5*np.eye(2)[None]*(k*k*gs+D[:,0,0]+D[:,1,1])[:,None,None])
    return Th
def stadium(P0,P1,rho,nq=6,narc=24):
    P0=np.asarray(P0,float);P1=np.asarray(P1,float);L=np.linalg.norm(P1-P0);u=(P1-P0)/L;n=np.array([-u[1],u[0]])
    P=[];N_=[];DL=[];LB=[]
    s,w=graded(L,rho,nq)
    for sg,lab in ((-1,'outer'),(1,'inner')):
        P.extend(P0+sg*rho*n+np.outer(s,u)); N_.extend([sg*n]*len(s)); DL.extend(w); LB.extend([lab]*len(s))
    xg,wg=np.polynomial.legendre.leggauss(narc)
    for C,th0,lab in ((P1,-pi/2,'tip'),(P0,pi/2,'vertex')):
        th=th0+(xg+1)/2*pi; d=np.cos(th)[:,None]*u+np.sin(th)[:,None]*n
        P.extend(C+rho*d); N_.extend(d); DL.extend(rho*wg/2*pi); LB.extend([lab]*narc)
    return np.array(P),np.array(N_),np.array(DL),np.array(LB)
def wing_parts(bodies,P0,P1,rho,lam):
    p,nn,dl,lb=stadium(P0,P1,rho)
    Tr=-np.einsum('mij,mj->mi',tensor_l(bodies,p,lam,smin=1/rho),nn)
    F={L:(Tr[lb==L]*dl[lb==L,None]).sum(0) for L in ('inner','outer','vertex','tip')}; F['wing']=(Tr*dl[:,None]).sum(0)
    return F,len(p)

def wing_parts_int(A1, A2, P0, P1, rho, lam):
    p, nn, dl, lb = stadium(P0, P1, rho)
    Th = tensor_l([A1, A2], p, lam, smin=1/rho) - tensor_l([A1], p, lam, smin=1/rho) - tensor_l([A2], p, lam, smin=1/rho)
    Tr = -np.einsum('mij,mj->mi', Th, nn)
    F = {L: (Tr[lb == L]*dl[lb == L, None]).sum(0) for L in ('inner', 'outer', 'vertex', 'tip')}; F['wing'] = (Tr*dl[:, None]).sum(0)
    return F




def rounded_polygon(V, rad, Nline, Narc=16):
    """Замкнутый CCW-многоугольник V с радиусами скругления rad[i] (0 = острый угол).
    Возвращает тело (mid, h, nrm) и описание сегментов для построения контуров."""
    V = [np.asarray(v, float) for v in V]; n = len(V); segs = []; corners = []
    for i in range(n):
        Pp, P, Pn = V[i-1], V[i], V[(i+1) % n]
        din = (P-Pp)/np.linalg.norm(P-Pp); dout = (Pn-P)/np.linalg.norm(Pn-P)
        phi = np.arctan2(din[0]*dout[1]-din[1]*dout[0], din@dout)      # >0 выпуклый
        rho = rad[i]; L = rho*np.tan(abs(phi)/2)
        Pa, Pb = P - din*L, P + dout*L
        left = np.array([-din[1], din[0]])
        C = Pa + np.sign(phi)*rho*left if rho > 0 else P
        corners.append(dict(P=P, Pa=Pa, Pb=Pb, C=C, rho=rho, phi=phi, din=din, dout=dout))
    M = []; H = []; NN = []; LB = []
    for i in range(n):
        c0, c1 = corners[i], corners[(i+1) % n]
        A, B = c0['Pb'], c1['Pa']; m, h = panels(A, B, Nline[i]); tau = (B-A)/np.linalg.norm(B-A)
        M.append(m); H.append(h); NN.append(np.tile([tau[1], -tau[0]], (len(h), 1))); LB += [('L', i)]*len(h)
        segs.append(('line', i, A, B))
        c = c1
        if c['rho'] > 0:
            a0 = np.arctan2(*(c['Pa']-c['C'])[::-1]); sgn = np.sign(c['phi'])
            th = a0 + sgn*np.linspace(0, abs(c['phi']), Narc+1)
            pts = c['C'] + c['rho']*np.stack([np.cos(th), np.sin(th)], 1)
            mm = (pts[1:]+pts[:-1])/2; hh = np.linalg.norm(pts[1:]-pts[:-1], axis=1)
            tt = (pts[1:]-pts[:-1])/hh[:, None]
            M.append(mm); H.append(hh); NN.append(np.stack([tt[:, 1], -tt[:, 0]], 1)); LB += [('A', (i+1) % n)]*Narc
    return (np.vstack(M), np.concatenate(H), np.vstack(NN)), corners, np.array(LB, dtype=object)

def offset_contour(corners, eps, nq=8, narc=24):
    """Контур на расстоянии eps снаружи тела; метки ('L',i) — грань i (V[i]->V[i+1]), ('A',i) — угол i."""
    n = len(corners); P = []; N_ = []; DL = []; LB = []
    for i in range(n):
        c0, c1 = corners[i], corners[(i+1) % n]
        A, B = c0['Pb'], c1['Pa']; Ls = np.linalg.norm(B-A); tau = (B-A)/Ls; nr = np.array([tau[1], -tau[0]])
        s, w = graded(Ls, eps, nq)
        P.extend(A + nr*eps + np.outer(s, tau)); N_.extend([nr]*len(s)); DL.extend(w); LB += [('L', i)]*len(s)
        c = c1; sgn = np.sign(c['phi'])
        Rr = c['rho'] + sgn*eps if c['rho'] > 0 else eps
        nin = np.array([c['din'][1], -c['din'][0]])
        a0 = np.arctan2(nin[1], nin[0]); xg, wg = np.polynomial.legendre.leggauss(narc)
        th = a0 + sgn*(xg+1)/2*abs(c['phi'])
        d = np.stack([np.cos(th), np.sin(th)], 1)
        if c['rho'] > 0 and sgn < 0:            # вогнутый угол: центр скругления снаружи тела
            P.extend(c['C'] - Rr*d)
        else:
            P.extend(c['C'] + Rr*d if c['rho'] > 0 else c['P'] + eps*d)
        N_.extend(d); DL.extend(Rr*wg/2*abs(c['phi'])); LB += [('A', (i+1) % n)]*narc
    return np.array(P), np.array(N_), np.array(DL), LB

def solid_corner(alpha_deg=15, t=0.08, r1=3.0, rho_in=0.02, rho_out=0.02, rho_tip=0.0, Nl=140, Nt=16, Narc=16):
    a = np.radians(alpha_deg); u1 = np.array([np.sin(a), np.cos(a)]); u2 = np.array([-np.sin(a), np.cos(a)])
    g1 = np.array([-np.cos(a), np.sin(a)]); g2 = np.array([np.cos(a), np.sin(a)])     # в щель от крыла 1 / крыла 2
    rI = t/(2*np.tan(a)); I = rI*u1 + g1*t/2; O = -rI*u1 - g1*t/2
    V = [O, r1*u1 - g1*t/2, r1*u1 + g1*t/2, I, r1*u2 + g2*t/2, r1*u2 - g2*t/2]
    names = ['внешняя грань 1', 'торец 1', 'внутренняя грань 1', 'внутренняя грань 2', 'торец 2', 'внешняя грань 2']
    cnames = ['внешний угол (вершина)', 'кромка торца 1 внеш.', 'кромка торца 1 внутр.', 'внутренний угол (дно)', 'кромка торца 2 внутр.', 'кромка торца 2 внеш.']
    rad = [rho_out, rho_tip, rho_tip, rho_in, rho_tip, rho_tip]
    body, corners, lb = rounded_polygon(V, rad, [Nl, Nt, Nl, Nl, Nt, Nl], Narc)
    return body, corners, names, cnames

def chunked(fun, bodies, pts, chunk=120, **kw):
    return np.vstack([fun(bodies, pts[i:i+chunk], **kw) for i in range(0, len(pts), chunk)])

def circle_flux(body, center, Rc, pol, M=400, nk=60):
    th = 2*pi*np.arange(M)/M; n = np.stack([np.cos(th), np.sin(th)], 1); pts = center + Rc*n
    if pol == 'TM': Th = tensor_l([body[:2]], pts, None, nk=nk)
    else:           Th = chunked(tensorN, [body], pts, nk=nk, smin=None)
    return -(np.einsum('mij,mj->mi', Th, n)*(2*pi*Rc/M)).sum(0)

def parts_flux(body, corners, eps, pol, nk=60):
    p, nn, dl, lb = offset_contour(corners, eps)
    if pol == 'TM': Th = tensor_l([body[:2]], p, None, nk=nk, smin=1/eps)
    else:           Th = chunked(tensorN, [body], p, nk=nk, smin=1/eps)
    Tr = -np.einsum('mij,mj->mi', Th, nn)
    out = {}
    for k in set(lb):
        msk = np.array([l == k for l in lb]); out[k] = (Tr[msk]*dl[msk, None]).sum(0)
    return out, (Tr*dl[:, None]).sum(0)
