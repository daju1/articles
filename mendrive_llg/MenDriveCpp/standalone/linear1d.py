# -*- coding: utf-8 -*-
"""
linear1d.py -- ЛИНЕЙНЫЙ эталон для сверки нелинейной 1D FDTD-модели МенДрайв (standalone/mendrive)
с линейными (частотными) моделями из ноутбуков MenDrive_qnm*.ipynb.

Конвенции -- КАК В НОУТБУКАХ MenDrive_qnm*.ipynb:
  * временная зависимость  F(x,t) = Re[ F~(x) * exp(-i*omega*t) ]   (exp(-I*omega_t) в ячейках 308/314/586);
  * Гауссова система, уравнения (ячейка 74 «дуальность уравнений Максвелла»)
        rot H = -i*omega/c * eps' * E + 4*pi/c * j_e ,   eps' = eps + 4*pi*i*sigma_e/omega
        rot E = +i*omega/c * mu'  * H - 4*pi/c * j_m ,   mu'  = mu  + 4*pi*i*sigma_m/omega
    здесь c = 1 (единицы FDTD: время измеряется как c*t, omega -- в 1/длина);
  * геометрия: левая стенка x < -a, вакуум -a <= x <= a, правая стенка (феррит) x > a;
  * анзацы в стенках (ячейки 47/48/276/277):  левая ~ exp(-i*K*(x+a)),  правая ~ exp(+i*K*(x-a)) = exp(-i*K*(a-x)),
    Im K > 0 (Зоммерфельд, §20; sign_K_* = 1, ячейка 263);
  * 1D (k_y = k_z = 0): пара TE = (Ey, Hz)  [TE относительно z: Hz != 0, Ez = 0],
                        пара TM = (Ez, Hy)  [TM: Ez != 0, Hz = 0].

Источник FDTD (src_k_phase = k):  j(x,t) = amp * ramp * exp(-xi/delta) * sin(omega*t - k*xi),
xi -- глубина в стенке (xi = x - a справа, xi = -a - x слева). Комплексная амплитуда:
    j~(x) = i*amp * exp(i*K_s*xi),   K_s = k + i/delta,
т.е. в обозначениях ноутбука j0_ext = i*amp, профиль exp(-i*K*(x+a)) слева / exp(+i*K*(x-a)) справа.

Режимы возбуждения (excitation_mode C++ <-> driving_mode ноутбуков):
  'je_z_left'  (TM) <-> driving_mode_0/_2 'je_z, k_y=0, TM' / 'je_z, TM'  (MenDrive_qnm-mode2_je_z_TM.ipynb)
  'je_y_left'  (TE) <-> driving_mode_1    'je_y, k_y=0, TE'               (MenDrive_qnm.ipynb);  legacy 'electric_left'
  'jm_z_right' (TE) <-> driving_mode_3/_4 'jm_z, TE' / 'jm_z, k_y=0, TE'  (MenDrive_qnm-mode3_jm_z_TE.ipynb); legacy 'magnetic_right'
  'jm_y_right' (TM) <-> driving_mode_5    'jm_y, k_y=0, TM'               (пока только в mode3, ячейка 21)
  + симметричные je_y_right, je_z_right, jm_z_left, jm_y_left (в ноутбуках не используются).

Нужен только numpy.
"""
import cmath
import math
import os

import numpy as np

# имя -> (вид тока 'e'/'m', компонента 'y'/'z', стенка 'left'/'right')
EXCITATIONS = {
    'magnetic_right': ('m', 'z', 'right'),   # legacy-имя (= jm_z_right)
    'electric_left':  ('e', 'y', 'left'),    # legacy-имя (= je_y_left, TE -- НЕ je_z!)
    'jm_z_right': ('m', 'z', 'right'),
    'je_y_left':  ('e', 'y', 'left'),
    'je_z_left':  ('e', 'z', 'left'),
    'jm_y_right': ('m', 'y', 'right'),
    'je_y_right': ('e', 'y', 'right'),
    'je_z_right': ('e', 'z', 'right'),
    'jm_z_left':  ('m', 'z', 'left'),
    'jm_y_left':  ('m', 'y', 'left'),
}
DRIVING_MODE_OF = {
    'je_z_left': ["driving_mode_0 = 'je_z, k_y=0, TM'", "driving_mode_2 = 'je_z, TM'"],
    'je_y_left': ["driving_mode_1 = 'je_y, k_y=0, TE'"],
    'jm_z_right': ["driving_mode_3 = 'jm_z, TE'", "driving_mode_4 = 'jm_z, k_y=0, TE'"],
    'jm_y_right': ["driving_mode_5 = 'jm_y, k_y=0, TM'"],
}


def polarization_of(mode):
    """je_y, jm_z -> TE (Ey,Hz);  je_z, jm_y -> TM (Ez,Hy)."""
    kind, comp, side = EXCITATIONS[mode]
    return 'TE' if (kind, comp) in (('e', 'y'), ('m', 'z')) else 'TM'


def sqrt_im_pos(z):
    """Корень с Im > 0 (при Im=0 -- с Re >= 0), как в ноутбуках."""
    r = cmath.sqrt(z)
    if r.imag < 0 or (r.imag == 0 and r.real < 0):
        r = -r
    return r


def wall_K(omega, eps=1.0, sigma_e=0.0, mu=1.0, sigma_m=0.0, k_t2=0.0):
    """K = sqrt(omega^2 eps' mu' - k_t^2), Im K > 0 (c = 1). k_t2 = k_y^2 + k_z^2 (в 1D FDTD = 0)."""
    epsp = eps + 4j * math.pi * sigma_e / omega
    mup = mu + 4j * math.pi * sigma_m / omega
    return sqrt_im_pos(omega ** 2 * epsp * mup - k_t2)


def source_params_from_K(K):
    """Параметры FDTD-источника, воспроизводящие профиль exp(i K xi): skin_depth = 1/Im K, src_k_phase = Re K."""
    return dict(skin_depth=1.0 / K.imag, src_k_phase=K.real)


class Layer:
    def __init__(self, x0, x1, eps=1.0, sigma_e=0.0, mu_yy=1.0, mu_zz=1.0, sigma_m=0.0):
        self.x0, self.x1 = x0, x1
        self.eps, self.sigma_e, self.mu_yy, self.mu_zz, self.sigma_m = eps, sigma_e, mu_yy, mu_zz, sigma_m


def _A_matrix(pol, omega, L):
    """u = [E, H];  du/dx = A u + b(x).  TM: (Ez,Hy), TE: (Ey,Hz)."""
    epsp = L.eps + 4j * math.pi * L.sigma_e / omega
    if pol == 'TM':
        mup = L.mu_yy + 4j * math.pi * L.sigma_m / omega
        # dEz/dx = -i w mu Hy + 4pi jm_y ;  dHy/dx = -i w eps Ez + 4pi je_z
        return np.array([[0, -1j * omega * mup], [-1j * omega * epsp, 0]], dtype=complex), epsp, mup
    mup = L.mu_zz + 4j * math.pi * L.sigma_m / omega
    # dEy/dx = +i w mu Hz - 4pi jm_z ;  dHz/dx = +i w eps Ey - 4pi je_y
    return np.array([[0, 1j * omega * mup], [1j * omega * epsp, 0]], dtype=complex), epsp, mup


def _b_vector(pol, kind):
    if pol == 'TM':
        return np.array([4 * math.pi, 0], dtype=complex) if kind == 'm' else np.array([0, 4 * math.pi], dtype=complex)
    return np.array([-4 * math.pi, 0], dtype=complex) if kind == 'm' else np.array([0, -4 * math.pi], dtype=complex)


def solve_linear(omega, a, h_l, h_r, mode, amp=1.0, skin_depth=0.06, src_k_phase=0.0,
                 sigma_e_left=3.0, sigma_e_ferrite=0.0, eps_left=1.0, eps_ferrite=1.0,
                 mu_ferrite_yy=1.0, mu_ferrite_zz=1.0, sigma_m_ferrite=0.0, outer='pec'):
    """
    Точное решение линейной частотной задачи (кусочно-однородные слои + экспоненциальный источник).
    Возвращает fields(x) -> dict комплексных профилей {'Ey','Ez','Hy','Hz','je_y','je_z','jm_y','jm_z'}.
    outer='pec'       -- E = 0 на x=-(a+h_l), x=a+h_r (как граничные узлы FDTD);
    outer='halfspace' -- полубесконечные стенки, только затухающие волны (анзацы ноутбуков).
    """
    kind, comp, side = EXCITATIONS[mode]
    pol = polarization_of(mode)
    layers = [Layer(-(a + h_l), -a, eps_left, sigma_e_left, 1.0, 1.0, 0.0),
              Layer(-a, a),
              Layer(a, a + h_r, eps_ferrite, sigma_e_ferrite, mu_ferrite_yy, mu_ferrite_zz, sigma_m_ferrite)]
    if outer == 'halfspace':
        layers[0].x0, layers[2].x1 = -np.inf, np.inf
    Ks = src_k_phase + 1j / skin_depth
    j0 = 1j * amp
    src_layer = 0 if side == 'left' else 2
    if side == 'right':
        q, C = 1j * Ks, j0 * cmath.exp(-1j * Ks * a)       # xi = x - a
    else:
        q, C = -1j * Ks, j0 * cmath.exp(-1j * Ks * a)      # xi = -a - x
    b0 = _b_vector(pol, kind) * C
    info = []
    for il, L in enumerate(layers):
        A, epsp, mup = _A_matrix(pol, omega, L)
        K = sqrt_im_pos(omega ** 2 * epsp * mup)
        lam = 1j * K
        vp = np.array([A[0, 1], lam], dtype=complex)
        vm = np.array([A[0, 1], -lam], dtype=complex)
        up = np.linalg.solve(q * np.eye(2) - A, b0) if il == src_layer else None
        info.append(dict(L=L, A=A, K=K, vp=vp, vm=vm, up=up))

    def basis(il, x):
        d = info[il]; L = d['L']; K = d['K']
        cp = d['vp'] * (np.exp(1j * K * (x - L.x0)) if np.isfinite(L.x0) else 0.0)
        cm = d['vm'] * (np.exp(-1j * K * (x - L.x1)) if np.isfinite(L.x1) else 0.0)
        return np.stack([cp, cm], axis=-1)

    def part(il, x):
        d = info[il]
        return np.zeros(2, dtype=complex) if d['up'] is None else d['up'] * np.exp(q * x)

    M = np.zeros((6, 6), dtype=complex); r = np.zeros(6, dtype=complex); row = 0
    if outer == 'pec':
        M[row, 0:2] = basis(0, layers[0].x0)[0]; r[row] = -part(0, layers[0].x0)[0]
    else:
        M[row, 0] = 1.0
    row += 1
    for il, xb in ((0, -a), (1, a)):
        Bl = basis(il, xb); Br = basis(il + 1, xb)
        for ci in (0, 1):
            M[row, 2 * il:2 * il + 2] = Bl[ci]
            M[row, 2 * il + 2:2 * il + 4] = -Br[ci]
            r[row] = part(il + 1, xb)[ci] - part(il, xb)[ci]
            row += 1
    if outer == 'pec':
        M[row, 4:6] = basis(2, layers[2].x1)[0]; r[row] = -part(2, layers[2].x1)[0]
    else:
        M[row, 5] = 1.0
    c = np.linalg.solve(M, r)

    def fields(x):
        x = np.atleast_1d(np.asarray(x, dtype=float))
        E = np.zeros(x.shape, dtype=complex); H = np.zeros(x.shape, dtype=complex); J = np.zeros(x.shape, dtype=complex)
        for k, xv in enumerate(x):
            il = 0 if xv < -a else (1 if xv <= a else 2)
            u = basis(il, xv) @ c[2 * il:2 * il + 2] + part(il, xv)
            E[k], H[k] = u[0], u[1]
            if il == src_layer:
                xi = (xv - a) if side == 'right' else (-a - xv)
                J[k] = j0 * cmath.exp(1j * Ks * xi)
        z = np.zeros(x.shape, dtype=complex)
        out = {k: z.copy() for k in ('Ey', 'Ez', 'Hy', 'Hz', 'je_y', 'je_z', 'jm_y', 'jm_z')}
        if pol == 'TM':
            out['Ez'], out['Hy'] = E, H
        else:
            out['Ey'], out['Hz'] = E, H
        out[('je_' if kind == 'e' else 'jm_') + comp] = J
        return out

    fields.info = dict(pol=pol, K=[d['K'] for d in info], Ks=Ks, j0=j0, coeffs=c,
                       residual=float(np.max(np.abs(M @ c - r))))
    return fields


def ode_residual(fields, omega, a, h_l, h_r, mode, n=4001, sigma_e_left=3.0, sigma_e_ferrite=0.0,
                 eps_left=1.0, eps_ferrite=1.0, mu_ferrite_yy=1.0, mu_ferrite_zz=1.0, sigma_m_ferrite=0.0):
    """Невязка уравнений Максвелла u' = A u + b (центральные разности по x, вдали от границ слоёв),
    отнесённая к max(|E|,|H|). Для точного решения ~ O(dx^2) (около 1e-6 при n=4001)."""
    pol = polarization_of(mode)
    kind, comp, side = EXCITATIONS[mode]
    x = np.linspace(-(a + h_l), a + h_r, n)
    F = fields(x)
    E = F['Ez'] if pol == 'TM' else F['Ey']
    H = F['Hy'] if pol == 'TM' else F['Hz']
    J = F[('je_' if kind == 'e' else 'jm_') + comp]
    dx = x[1] - x[0]
    scale = max(np.max(np.abs(E)), np.max(np.abs(H)), 1e-300)
    layers = [Layer(0, 0, eps_left, sigma_e_left), Layer(0, 0),
              Layer(0, 0, eps_ferrite, sigma_e_ferrite, mu_ferrite_yy, mu_ferrite_zz, sigma_m_ferrite)]
    res = 0.0
    for k in range(1, n - 1):
        xv = x[k]
        if abs(abs(xv) - a) < 2 * dx:
            continue
        il = 0 if xv < -a else (1 if xv <= a else 2)
        A, _, _ = _A_matrix(pol, omega, layers[il])
        b = _b_vector(pol, kind) * J[k]
        du = np.array([(E[k + 1] - E[k - 1]) / (2 * dx), (H[k + 1] - H[k - 1]) / (2 * dx)])
        res = max(res, float(np.max(np.abs(du - (A @ np.array([E[k], H[k]]) + b)))))
    return res / scale


def notebook_source_profile(x, mode, j0, K, a, sign_K=1, fix_right_typo=False):
    """
    Профили источников ровно по формулам ячейки 21 MenDrive_qnm-mode3_jm_z_TE.ipynb (phase_y = 1):
        левая стенка  (je_ext_z, je_ext_y): j0*exp(-I*sign_K*K*(x+a))
        правая стенка (jm_ext_z, jm_ext_y): j0*exp(-I*sign_K*K*(x+a))   <- как НАПИСАНО (вероятная опечатка)
    fix_right_typo=True: правая стенка j0*exp(-I*sign_K*K*(a-x)) -- анзац поля правой стенки (ячейки 48/276/279/330).
    """
    x = np.asarray(x, dtype=float)
    kind, comp, side = EXCITATIONS[mode]
    if side == 'left' or not fix_right_typo:
        return j0 * np.exp(-1j * sign_K * K * (x + a))
    return j0 * np.exp(-1j * sign_K * K * (a - x))


def write_profiles_csv(path, x, prof, x_cm_scale=1.0, n=1, omega=float('nan'), note=''):
    """Формат сравнения (= fields_harm_n<k>.csv программы mendrive):
    n,omega_n,x,x_cm,<F>_re,<F>_im,<F>_abs,<F>_arg  для каждого F из prof."""
    names = list(prof.keys())
    with open(path, 'w', encoding='utf-8') as f:
        for line in (note.splitlines() if note else []):
            f.write('# ' + line + '\n')
        head = ['n', 'omega_n', 'x', 'x_cm']
        for nm in names:
            head += [nm + '_re', nm + '_im', nm + '_abs', nm + '_arg']
        f.write(','.join(head) + '\n')
        for k, xv in enumerate(np.atleast_1d(x)):
            row = [str(n), repr(float(n * omega)), repr(float(xv)), repr(float(xv * x_cm_scale))]
            for nm in names:
                v = complex(np.atleast_1d(prof[nm])[k])
                row += [repr(v.real), repr(v.imag), repr(abs(v)), repr(cmath.phase(v))]
            f.write(','.join(row) + '\n')


def read_profiles_csv(path):
    """Читает CSV сравнения / fields_harm_*.csv / fields_snap.csv; строки '#' пропускаются."""
    cols = None; data = []
    with open(path, encoding='utf-8') as f:
        for line in f:
            if line.startswith('#') or not line.strip():
                continue
            parts = line.strip().split(',')
            if cols is None:
                cols = parts
            else:
                data.append(parts)
    out = {}
    for j, cname in enumerate(cols or []):
        vals = []
        for r in data:
            try:
                vals.append(float(r[j]))
            except (ValueError, IndexError):
                vals.append(float('nan'))
        out[cname] = np.array(vals)
    return out


def complex_column(d, name):
    return d[name + '_re'] + 1j * d[name + '_im']


def compare_profiles(x_lin, f_lin, x_nl, f_nl, norm='max', fit='phase', xlim=None):
    """
    Сравнение профилей: нелинейный интерполируется на x_lin (в общем диапазоне), нормировка по max|F|,
    подбор общей фазы (fit='phase') или общего комплексного множителя (fit='complex').
    rel_L2 = ||alpha*v - u|| / ||u||.  Возвращает (rel_L2, alpha, x, u_lin, alpha*v_nl).
    """
    x_lin = np.asarray(x_lin, dtype=float); x_nl = np.asarray(x_nl, dtype=float)
    o = np.argsort(x_nl); x_nl = x_nl[o]; f_nl = np.asarray(f_nl)[o]
    m = (x_lin >= x_nl[0]) & (x_lin <= x_nl[-1])
    if xlim is not None:
        m &= (x_lin >= xlim[0]) & (x_lin <= xlim[1])
    xs = x_lin[m]; u = np.asarray(f_lin)[m]
    v = np.interp(xs, x_nl, f_nl.real) + 1j * np.interp(xs, x_nl, f_nl.imag)
    if norm == 'max':
        u = u / max(np.max(np.abs(u)), 1e-300); v = v / max(np.max(np.abs(v)), 1e-300)
    s = np.vdot(v, u)
    if fit == 'complex':
        alpha = s / max(np.vdot(v, v).real, 1e-300)
    elif fit == 'phase':
        alpha = np.exp(1j * np.angle(s)) if abs(s) > 0 else 1.0
    else:
        alpha = 1.0
    r = np.linalg.norm(alpha * v - u) / max(np.linalg.norm(u), 1e-300)
    return float(r), alpha, xs, u, alpha * v


# ------------------------------------------------------------------------------------------
# Экспорт профилей ИЗ ноутбуков MenDrive_qnm*.ipynb (SageMath) в CSV сравнения
# ------------------------------------------------------------------------------------------
def export_sage_profiles(path, fields, A, h_l, h_r, y0=0.0, n_per_segment=400, omega=float('nan'),
                         x_cm_scale=1.0, j0=None, note=''):
    """
    Вызывать В КОНЦЕ ноутбука MenDrive_qnm*.ipynb (после ячейки subs_fields(...), где определены
    Ey_left_conductor_d, Ey_vacuum_d, Ey_right_conductor_d, ...). Пример:

        import sys; sys.path.insert(0, '/Users/oleksiy.drozdov/articles/mendrive_llg/MenDriveCpp/standalone')
        from linear1d import export_sage_profiles
        export_sage_profiles('linear_mode3.csv',
            {'Ey': (Ey_left_conductor_d, Ey_vacuum_d, Ey_right_conductor_d),
             'Hz': (Hz_left_conductor_d, Hz_vacuum_d, Hz_right_conductor_d),
             'Ez': (Ez_left_conductor_d, Ez_vacuum_d, Ez_right_conductor_d),
             'Hy': (Hy_left_conductor_d, Hy_vacuum_d, Hy_right_conductor_d)},
            A=A, h_l=h_conductor_l, h_r=h_conductor_r, y0=0, omega=omega_re_value)

    fields: имя -> (f_left, f_vacuum, f_right) -- символьные выражения Sage от x (и, возможно, y)
            или обычные функции Python f(x). Сегменты: [-A-h_l, -A], [-A, A], [A, A+h_r] (как в
            ячейках 521-537 ноутбуков). Значения в точке y=y0 (phase_y). x -- в см ноутбука.
    j0: если задано (например, je0_ext или jm0_ext), профили делятся на j0 (нормировка
        'source', сравнимая с fields_harm_normalize=source программы mendrive).
    """
    def ev(f, xv):
        try:
            if hasattr(f, 'subs'):
                try:
                    v = f.subs(x=xv, y=y0)
                except Exception:
                    v = f.subs(x=xv)
                try:
                    v = v.n()
                except Exception:
                    pass
                return complex(v)
            return complex(f(xv))
        except Exception:
            return complex(float('nan'), float('nan'))

    segs = [(-A - h_l, -A), (-A, A), (A, A + h_r)]
    xs = []
    for k, (x0, x1) in enumerate(segs):
        xx = np.linspace(float(x0), float(x1), n_per_segment)
        xs.append(xx if k == 0 else xx[1:])
    prof = {nm: [] for nm in fields}
    xall = []
    for k, xx in enumerate(xs):
        for xv in xx:
            xall.append(float(xv))
            for nm, fl in fields.items():
                prof[nm].append(ev(fl[k], float(xv)))
    prof = {nm: np.array(v) for nm, v in prof.items()}
    if j0 is not None:
        prof = {nm: v / complex(j0) for nm, v in prof.items()}
    hdr = (f'linear profile (MenDrive_qnm*): a={float(A)!r} h_l={float(h_l)!r} h_r={float(h_r)!r} y0={y0} '
           f'time_convention=minus normalize={"source" if j0 is not None else "none"}')
    write_profiles_csv(path, np.array(xall), prof, x_cm_scale=x_cm_scale, n=1, omega=float(omega),
                       note=hdr + ('\n' + note if note else ''))
    return path


# ------------------------------------------------------------------------------------------
# Точное решение для параметров прогона ./mendrive (out_dir/config_used.txt)
# ------------------------------------------------------------------------------------------
def read_kv(path):
    d = {}
    if os.path.exists(path):
        with open(path, encoding='utf-8') as f:
            for line in f:
                line = line.split('#', 1)[0].strip()
                if '=' in line:
                    k, v = line.split('=', 1)
                    d[k.strip()] = v.strip()
    return d


def analytic_for_run(out_dir, n_points=2001, outer='pec', path=None):
    """
    Точный линейный профиль (первая гармоника) для параметров прогона из out_dir/config_used.txt
    (omega -- omega_res из summary.txt). Феррит: Linear -> mu_yy_lin, mu_zz_lin, sigma_m_leak;
    для JA/Preisach/LLG -- mu = 1 (только TM-режимы со скалярными моделями сравнимы строго,
    для TE -- лишь качественно, при малом amp). Пишет CSV сравнения и возвращает путь.
    """
    c = read_kv(os.path.join(out_dir, 'config_used.txt'))
    s = read_kv(os.path.join(out_dir, 'summary.txt'))
    g = lambda k, dflt: float(c.get(k, dflt))
    omega = float(s.get('omega_res', c.get('omega0', 'nan')))
    lin = c.get('ferrite_model') == 'Linear'
    mode = c.get('excitation_mode', 'magnetic_right')
    a, h_l, h_r = g('a', 0.2), g('h_l', 0.2), g('h_r', 0.2)
    f = solve_linear(omega, a, h_l, h_r, mode, amp=g('amp', 1.0), skin_depth=g('skin_depth', 0.06),
                     src_k_phase=g('src_k_phase', 0.0), sigma_e_left=g('sigma_e_left', 3.0),
                     sigma_e_ferrite=g('sigma_e_ferrite', 0.0), eps_left=g('eps_left', 1.0),
                     eps_ferrite=g('eps_ferrite', 1.0),
                     mu_ferrite_yy=g('mu_yy_lin', 1.0) if c.get('ferrite_model') not in ('LLG', 'Hybrid') else 1.0,
                     mu_ferrite_zz=g('mu_zz_lin', 1.0) if lin else 1.0,
                     sigma_m_ferrite=g('sigma_m_leak', 0.0) if lin else 0.0, outer=outer)
    x = np.linspace(-(a + h_l), a + h_r, n_points)
    F = f(x)
    norm = c.get('fields_harm_normalize', 'none') == 'source'
    if norm:
        F = {k: v / f.info['j0'] for k, v in F.items()}
    path = path or os.path.join(out_dir, 'linear_analytic_n1.csv')
    write_profiles_csv(path, x, F, x_cm_scale=g('length_unit_cm', 1.0), n=1, omega=omega,
                       note=(f'linear1d.solve_linear: a={a!r} h_l={h_l!r} h_r={h_r!r} mode={mode} outer={outer} '
                             f'time_convention=minus normalize={"source" if norm else "none"} '
                             f'K_walls={f.info["K"]} Ks={f.info["Ks"]}'))
    return path


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description='Точный линейный 1D-профиль для прогона ./mendrive (CSV сравнения)')
    ap.add_argument('out_dir', help='папка прогона (config_used.txt, summary.txt)')
    ap.add_argument('--points', type=int, default=2001)
    ap.add_argument('--outer', default='pec', choices=['pec', 'halfspace'])
    ap.add_argument('-o', '--out', default=None)
    a_ = ap.parse_args()
    print('записано:', analytic_for_run(a_.out_dir, a_.points, a_.outer, a_.out))
