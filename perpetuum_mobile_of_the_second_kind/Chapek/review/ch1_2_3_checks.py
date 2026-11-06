# Численные проверки к разбору глав 1–3 Čápek & Sheehan (2005).
# Только numpy. Единицы: hbar = 1, k_B = 1.
import numpy as np

# ---------- 1. Fish-trap: уравнения (3.57)-(3.61) в точности как напечатаны (стр. 78-80) ----------
def capek_fishtrap(eps=1.0, J=0.1, I=0.1, Gd=0.5, beta=3.0):
    Gu = Gd*np.exp(-beta*eps)
    k = -1j*eps - 0.5*(Gu+Gd); kc = np.conj(k); iI, iJ = 1j*I, 1j*J
    A = np.zeros((9, 9), complex); D = np.zeros((9, 9), complex)
    A[0,0]=-Gd
    A[1,1]=-Gu; A[1,6]=iI; A[1,8]=-iI
    A[2,2]=-Gd; A[2,6]=-iI; A[2,8]=iI
    A[3,3]=k;  A[3,4]=iI
    A[4,3]=iI; A[4,4]=-Gd
    A[5,5]=kc; A[5,7]=-iI
    A[6,1]=iI; A[6,2]=-iI; A[6,6]=kc
    A[7,5]=-iI; A[7,7]=-Gd
    A[8,1]=-iI; A[8,2]=iI; A[8,8]=k
    D[0,0]=-Gu; D[0,3]=iJ; D[0,5]=-iJ
    D[1,1]=-Gd; D[1,3]=-iJ; D[1,5]=iJ
    D[2,2]=-Gu
    D[3,0]=iJ; D[3,1]=-iJ; D[3,3]=kc
    D[4,4]=-Gu; D[4,6]=-iJ
    D[5,0]=-iJ; D[5,1]=iJ; D[5,5]=k
    D[6,4]=-iJ; D[6,6]=k
    D[7,7]=-Gu; D[7,8]=iJ
    D[8,7]=iJ; D[8,8]=kc
    av = (Gu+Gd)/2
    B = np.diag([Gu,Gd,Gu,av,Gu,av,av,Gu,av]); C = np.diag([Gd,Gu,Gd,av,Gd,av,av,Gd,av])
    return np.block([[A, B], [C, D]])

POP = [0, 1, 2, 9, 10, 11]; NAMES = ['-1u','0u','1u','-1d','0d','1d']
def steady(M, popidx=POP):
    w, v = np.linalg.eig(M); x = v[:, np.argmin(abs(w))]
    return x/sum(x[j] for j in popidx)

def HS_fish(eps=1.0, J=0.1, I=0.1):
    H = np.zeros((6, 6))
    for s, m in enumerate([-1, 0, 1]):
        f = -1 if m == 0 else 1
        H[s, s] = eps/2*f; H[3+s, 3+s] = -eps/2*f
    H[3, 4] = H[4, 3] = J; H[1, 2] = H[2, 1] = I
    return H

def gibbs_pops(H, beta):
    w, v = np.linalg.eigh(H); p = np.exp(-beta*(w-w.min())); p /= p.sum()
    return np.diag((v*p) @ v.T)

if __name__ == "__main__":
    print("== 1. Fish-trap: стационар уравнений книги vs Гиббс ==")
    for beta in [0.5, 1, 2, 3, 5]:
        P = np.real(steady(capek_fishtrap(beta=beta))[POP])
        print(f"beta*eps={beta}: P1d={P[5]:.4f}, P1d/P-1d={P[5]/P[3]:.3f}, exp(2*beta*eps)={np.exp(2*beta):.3f}, "
              f"Gibbs P1d={gibbs_pops(HS_fish(), beta)[5]:.4f}")

    print("== 2. Скейлинг ван Хова: J, I, Gamma ~ lam^2, eps фиксировано ==")
    for lam in [0.3, 0.2, 0.1, 0.05]:
        w = sorted(np.linalg.eigvals(capek_fishtrap(J=lam**2, I=lam**2, Gd=lam**2, beta=3.0)), key=abs)
        print(f"lam={lam}: slowest rate/lam^2={abs(w[1].real)/lam**2:.2e}, /lam^6={abs(w[1].real)/lam**6:.3f}")

# ---------- 3. Периодический fish-trap (винт), редуцированная сеть скоростей ----------
def ring_current(Delta, beta=3.0, eps=1.0, X=0.05, Gd=0.5, Gam=1.0, ldb=False):
    Gu = Gd*np.exp(-beta*eps); gc = 2*Gam + 0.5*(Gu+Gd)
    W = lambda d: 2*X**2*gc/(gc**2 + d**2)       # лоренцева (симметричная!) скорость когерентного прыжка
    links = []
    def coh(i, j, dE):
        w = W(dE)
        f, b = (w*2/(1+np.exp(beta*dE)), w*2/(1+np.exp(-beta*dE))) if ldb else (w, w)
        links.append((i, j, f, b))
    # состояния витка: 0:-1d 1:0d 2:0u 3:1u 4:1d 5:-1u
    coh(0, 1, eps+Delta/3); coh(2, 3, eps+Delta/3); coh(4, 0, Delta/3); coh(3, 5, Delta/3)
    links += [(1, 2, Gd, Gu), (3, 4, Gd, Gu), (5, 0, Gd, Gu)]  # переключения ворот (KMS)
    L = np.zeros((6, 6))
    for i, j, f, b in links:
        L[j, i] += f; L[i, i] -= f; L[i, j] += b; L[j, j] -= b
    w, v = np.linalg.eig(L); p = np.real(v[:, np.argmin(abs(w))]); p /= p.sum()
    turns = sum(p[i]*f - p[j]*b for (i, j, f, b) in links if (i, j) in [(4, 0), (3, 5)])
    lnW = sum(np.log(f/b) for (i, j, f, b) in links if (i, j) in [(0,1),(1,2),(2,3),(3,4),(4,0)])
    return turns, lnW

# ---------- 4. Плазменный тепловой насос (3.105)-(3.111): локальное уравнение Линдблада ----------
def lindblad_ss(H, Ls):
    n = H.shape[0]; Id = np.eye(n)
    S = -1j*(np.kron(H, Id) - np.kron(Id, H.T))
    for L in Ls:
        LdL = L.conj().T @ L
        S += np.kron(L, L.conj()) - 0.5*np.kron(LdL, Id) - 0.5*np.kron(Id, LdL.T)
    w, v = np.linalg.eig(S); x = v[:, np.argmin(abs(w))].reshape(n, n)
    return x/np.trace(x)

def plasma(TI, TII, e2=2.0, e3=1.0, J=0.05, K=0.05, G0=0.2, deph0=1.0):
    bII = 1/TII
    Gd = G0*(1+1/(np.exp(bII*e3)-1)); Gu = G0/(np.exp(bII*e3)-1)
    nI = 1/(np.exp(1.0/TI)-1); gdep = deph0*nI*(1+nI)      # «рамановская» дефазировка, delta(w1-w2)
    H = np.diag([0, e2, e3]).astype(complex); H[0,1]=H[1,0]=J; H[1,2]=H[2,1]=K
    a = lambda i, j: np.outer(np.eye(3)[i], np.eye(3)[j]).astype(complex)
    r = np.real(np.diag(lindblad_ss(H, [np.sqrt(Gd)*a(0,2), np.sqrt(Gu)*a(2,0), np.sqrt(4*gdep)*a(1,1)])))
    return e3*(Gd*r[2] - Gu*r[0])

if __name__ == "__main__":
    print("== 3. Винт: число витков в ед. времени и ln W ==")
    for D in [0, 0.5, 1, 3, 6, 12]:
        a = ring_current(D); b = ring_current(D, ldb=True)
        print(f"Delta={D}: local {a[0]: .3e} (lnW={a[1]:.2f}) | LDB {b[0]: .3e} (lnW={b[1]:.2f}), -beta*Delta={-3*D}")
    print("== 4. Плазма: Q(I->II) ==")
    for TI, TII in [(0.5, 2.0), (1.0, 2.0), (2.0, 1.0), (1.0, 1.0), (0.3, 5.0)]:
        print(f"T_I={TI}, T_II={TII}: Q={plasma(TI, TII): .3e}")
