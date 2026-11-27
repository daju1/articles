"""corner_coupled.py — связанная ТМ+ТЕ задача Казимира для z-однородных тел из реального металла,
полный ЭМ-тензор из функции Грина в базисе (E_z,h). hbar=c=1. Основа — Casimir_cylinders_Maxwell_general.ipynb,
нормировка — Casimir_P3_Lifshitz_coupled.ipynb §6."""
import numpy as np, time
pi = np.pi
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

def poly_body(V, Ns):
    """Многоугольник с вершинами V (против часовой), Ns панелей на сторонах: (середины, длины, внешние нормали)."""
    M = []; H = []; NN = []
    for i in range(len(V)):
        X = np.asarray(V[i], float); Y = np.asarray(V[(i+1) % len(V)], float)
        m, h = panels(X, Y, Ns[i]); tau = (Y-X)/np.linalg.norm(Y-X)
        M.append(m); H.append(h); NN.append(np.tile([tau[1], -tau[0]], (Ns[i], 1)))
    return np.vstack(M), np.concatenate(H), np.vstack(NN)

def kgrid(n=60, s=1.0):
    """Узлы и веса для интеграла по kappa от 0 до бесконечности: kappa = s*tan(pi u/2)."""
    u = (np.arange(n) + 0.5)/n
    return s*np.tan(pi/2*u), s*pi/2/np.cos(pi/2*u)**2/n

# ---------------- материалы на мнимой частоте xi (eps, mu вещественны и >= 1) ----------------
class Material:
    """eps(xi), mu(xi) на мнимой частоте omega = i*xi. kind='pec' — идеальный проводник."""
    def __init__(self, kind='vacuum', eps=1.0, mu=1.0, wp=None, gamma=0.0, w0=None):
        self.kind, self.e0, self.m0, self.wp, self.g, self.w0 = kind, eps, mu, wp, gamma, w0
    def eps(self, xi):
        if self.kind in ('vacuum',): return 1.0
        if self.kind == 'dielectric': return self.e0
        if self.kind == 'plasma': return 1.0 + self.wp**2/xi**2
        if self.kind == 'drude':  return 1.0 + self.wp**2/(xi*(xi + self.g))
        if self.kind == 'lorentz': return 1.0 + (self.e0 - 1.0)*self.w0**2/(self.w0**2 + xi**2)
        raise ValueError(self.kind)
    def mu(self, xi): return 1.0 if self.kind != 'dielectric' else self.m0

def Kmat(mid, h, nrm, k):
    """Двойной слой: K_jb = h_b * d/dn_b G(m_j - m_b),  K_jj = 0 (прямая панель)."""
    d = mid[:, None, :] - mid[None, :, :]; R = np.linalg.norm(d, axis=2); np.fill_diagonal(R, 1.0)
    _, K1 = K0K1(k*R); gr = (-k*K1/(2*pi*R))[..., None]*d          # grad_x G(x - y)
    Km = -np.einsum('jbc,bc->jb', gr, nrm)*h[None, :]; np.fill_diagonal(Km, 0.0)
    return Km

def Dtau(mid, h, nrm, loops):
    """Производная по дуге d/ds (s вдоль tau = z x n) на замкнутых контурах: центральные разности
    по серединам панелей. loops — список массивов индексов панелей каждого замкнутого контура (по порядку)."""
    N = len(h); D = np.zeros((N, N))
    for L in loops:
        n = len(L)
        for a in range(n):
            j, jp, jm = L[a], L[(a+1) % n], L[(a-1) % n]
            sp = (h[j] + h[jp])/2; sm = (h[j] + h[jm])/2
            # трёхточечная формула на неравномерной сетке
            D[j, jp] += sm/(sp*(sp+sm)); D[j, jm] += -sp/(sm*(sp+sm)); D[j, j] += (sp - sm)/(sp*sm)
    return D



def hess(d, R, k):
    K0, K1 = K0K1(k*R)
    dd = d[..., :, None]*d[..., None, :]
    return (k**2*K0/R**2)[..., None, None]*dd/(2*pi) + (k*K1)[..., None, None]*(2*dd/R[..., None, None]**3 - np.eye(2)/R[..., None, None])/(2*pi)

class Coupled:
    """Связанная ТМ+ТЕ задача для нескольких z-однородных тел (материал mat в вакууме) при (xi,kz)."""
    def __init__(self, bodies, mat, xi, kz):
        self.mid = np.vstack([b[0] for b in bodies]); self.h = np.concatenate([b[1] for b in bodies])
        self.nrm = np.vstack([b[2] for b in bodies]); N = len(self.h); self.N = N
        idx = np.cumsum([0]+[len(b[1]) for b in bodies]); loops = [np.arange(idx[i], idx[i+1]) for i in range(len(bodies))]
        ei = mat.eps(xi); ko = np.sqrt(xi**2+kz**2); ki = np.sqrt(ei*xi**2+kz**2)
        So, Ko = Smat(self.mid, self.h, ko), Kmat(self.mid, self.h, self.nrm, ko)
        Si = np.zeros((N, N)); Ki = np.zeros((N, N))
        for L in loops:
            s = slice(L[0], L[-1]+1)
            Si[s, s] = Smat(self.mid[s], self.h[s], ki); Ki[s, s] = Kmat(self.mid[s], self.h[s], self.nrm[s], ki)
        Dt = Dtau(self.mid, self.h, self.nrm, loops); I = np.eye(N); Z = np.zeros((N, N))
        r2 = ki**2/ko**2; c = (kz/xi)*(1-r2); aE, bE = r2/ei, c/ei; ah, bh = r2, c
        self.A = np.block([[0.5*I-Ko, Z, So, Z], [Z, 0.5*I-Ko, Z, So],
                           [0.5*I+Ki, -Si@(bE*Dt), -Si*aE, Z], [-Si@(bh*Dt), 0.5*I+Ki, Z, -Si*ah]])
        self.lu = np.linalg.inv(self.A); self.xi, self.kz, self.ko = xi, kz, ko

    def incident(self, r):
        """(E_z, h) на контуре от 6 источников в точке r: эл. диполи x,y,z, магн. x,y,z.  -> (N,6),(N,6)"""
        xi, kz, k = self.xi, self.kz, self.ko
        d = self.mid - r[None, :]; R = np.linalg.norm(d, axis=1)
        K0, K1 = K0K1(k*R); g = K0/(2*pi); gx, gy = (-k*K1/(2*pi*R))*d[:, 0], (-k*K1/(2*pi*R))*d[:, 1]
        Ez = np.zeros((self.N, 6), complex); Hz = np.zeros((self.N, 6), complex)
        # электрический: E = (e - grad grad.e / xi^2) g ;  H_z = -curl_z E / xi
        Ez[:, 0] = -1j*kz*gx/xi**2; Ez[:, 1] = -1j*kz*gy/xi**2; Ez[:, 2] = k**2/xi**2*g
        Hz[:, 0] = gy/xi; Hz[:, 1] = -gx/xi
        # магнитный: H = (e - grad grad.e/xi^2) g ;  E_z = curl_z H / xi
        Hz[:, 3] = -1j*kz*gx/xi**2; Hz[:, 4] = -1j*kz*gy/xi**2; Hz[:, 5] = k**2/xi**2*g
        Ez[:, 3] = -gy/xi; Ez[:, 4] = gx/xi
        return Ez, 1j*Hz

    def scattered_EH(self, r):
        """Рассеянные полные векторы E,H в точке r от 6 источников в той же точке -> (3,6),(3,6)"""
        Einc, hinc = self.incident(r)
        x = self.lu @ np.vstack([Einc, hinc, np.zeros((2*self.N, 6))])
        N = self.N; E, hh, qE, qh = x[:N], x[N:2*N], x[2*N:3*N], x[3*N:]
        k = self.ko; d = r[None, :] - self.mid; R = np.linalg.norm(d, axis=1)
        K0, K1 = K0K1(k*R); G = K0/(2*pi); gr = (-k*K1/(2*pi*R))[:, None]*d       # grad_r G
        Hs = hess(d, R, k)                                                          # (B,2,2)
        dnG = -np.einsum('bc,bc->b', gr, self.nrm); grad_dnG = -np.einsum('bij,bj->bi', Hs, self.nrm)
        w = self.h
        val = lambda u, q: (dnG*w) @ u - (G*w) @ q
        grd = lambda u, q: np.einsum('bi,b,bs->is', grad_dnG, w, u) - np.einsum('bi,b,bs->is', gr, w, q)
        Ez, hz = val(E, qE), val(hh, qh); gE, gh = grd(E, qE), grd(hh, qh)
        xi, kz = self.xi, self.kz; k2 = k**2
        Ex = -1j*(kz*gE[0] + xi*gh[1])/k2; Ey = -1j*(kz*gE[1] - xi*gh[0])/k2
        Hx = -(kz*gh[0] + xi*gE[1])/k2; Hy = (-kz*gh[1] + xi*gE[0])/k2
        return np.array([Ex, Ey, Ez]), np.array([Hx, Hy, -1j*hz])

def T_point(bodies, mat, pts, nQ=16, nt=12, scale=1.0):
    """<T_ik>(r) (i,k = x,y,z) в точках pts, hbar=c=1:  T = -(1/2pi^2) Int dxi Int dkz xi^2 [S - 1/2 d trS]."""
    u = (np.arange(nQ)+0.5)/nQ; Q = scale*np.tan(pi/2*u); wQ = scale*pi/2/np.cos(pi/2*u)**2/nQ
    th = (np.arange(nt)+0.5)/nt*pi/2; wt = pi/2/nt
    T = np.zeros((len(pts), 3, 3))
    for Qi, wq in zip(Q, wQ):
        for t in th:
            xi, kz = Qi*np.cos(t), Qi*np.sin(t); sol = Coupled(bodies, mat, xi, kz)
            for m, r in enumerate(pts):
                Es, Hs = sol.scattered_EH(r)
                S = xi**2*(Es[:, :3] + Hs[:, 3:]); S = 0.5*(S+S.T)
                T[m] += 2*(S - 0.5*np.eye(3)*np.trace(S)).real*Qi*wq*wt   # x2: kz<0
    return -T/(2*pi**2)



def strip(y0, y1, L, n_long, n_short):
    V = [(-L/2, y0), (L/2, y0), (L/2, y1), (-L/2, y1)]
    return poly_body(V, [n_long, n_short, n_long, n_short])

def wing(alpha, r0, r1, t, sign, nl, ns):
    """Крыло толщины t от r0 до r1 вдоль луча под углом sign*alpha к оси x (вершина в начале координат),
    утолщение наружу от биссектрисы. Вершины многоугольника против часовой стрелки."""
    e = np.array([np.cos(alpha), sign*np.sin(alpha)]); nv = sign*np.array([-e[1], e[0]])
    P = [r0*e, r1*e, r1*e + t*nv, r0*e + t*nv]
    if sign < 0: P = P[::-1]
    return poly_body(P, [nl, ns, nl, ns])

def flux_circle(bodies, mat, R, n_ang=48, **kw):
    """F_i = -∮ T_ij n_j dl по окружности радиуса R с центром (r1/2,0)-независимо: центр в начале координат."""
    ph = (np.arange(n_ang)+0.5)/n_ang*2*pi; pts = R*np.stack([np.cos(ph), np.sin(ph)], 1)
    T = T_point(bodies, mat, pts, **kw); nrm = pts/R
    F = -np.einsum('mij,mj->i', T[:, :2, :2], nrm)*(2*pi*R/n_ang)
    return F, T


if __name__ == '__main__':
    import sys, time
    wp = float(sys.argv[1]) if len(sys.argv) > 1 else 40.0
    nQ, nt = 14, 10
    mat = Material('plasma', wp=wp)
    t0 = time.time()
    # эталон: бесконечные пластины, зазор 1, та же модель (Лифшиц)
    print('Эталон: пластины, плазма wp=%g: см. Casimir_P3_Lifshitz_coupled (wp=40: eta=0.8805, T_xx=T_zz=-0.644)' % wp)
    # --- тест 1: две полосы
    for L, nl in ((6, 60), (10, 90), (16, 140)):
        T = T_point([strip(-1.5, -0.5, L, nl, 8), strip(0.5, 1.5, L, nl, 8)], mat, np.array([[0, 0.0]]), nQ=nQ, nt=nt)
        print('полосы L=%2d: diag T/(pi^2/240) = %s   %.0f c' % (L, np.round(np.diag(T[0])/(pi**2/240), 4), time.time()-t0), flush=True)
    # --- тест 2: уголок alpha=15, r0=0.3, r1=3, t=0.08 — поток через окружности
    alpha, r0, r1, t = np.radians(15), 0.3, 3.0, 0.08
    for nl in (100, 140):
        bodies = [wing(alpha, r0, r1, t, +1, nl, 6), wing(alpha, r0, r1, t, -1, nl, 6)]
        for R in (3.5, 5.0):
            F, T = flux_circle(bodies, mat, R, n_ang=96, nQ=nQ, nt=nt)
            print('уголок nl=%d R=%.1f: F=(%+.3e, %+.3e), масштаб max|T|*2piR=%.3e   %.0f c'
                  % (nl, R, F[0], F[1], abs(T).max()*2*pi*R, time.time()-t0), flush=True)
