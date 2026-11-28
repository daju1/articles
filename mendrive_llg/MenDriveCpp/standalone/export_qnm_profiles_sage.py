# -*- coding: utf-8 -*-
"""
export_qnm_profiles_sage.py -- выгрузка профилей E~_y(x), H~_z(x) линейной модели MenDrive_qnm.ipynb
(driving_mode_1 = 'je_y, k_y=0, TE') в CSV формата linear1d.write_profiles_csv для
    python3 plot.py <out_dir> --compare <csv> --fit none --norm none

Нормировка CSV совпадает с ./mendrive fields_harm_normalize=source при L0 = 1 см (examples/qnm_drive_je_y_left_linear.cfg):
    F_csv = F~ / j0_code,   j0_code = je0_ext * L0 / c   (=> F_csv = c * F~/je0_ext),  x -- в см.

Три способа:

(1) ИЗ НОУТБУКА (ядро SageMath), после ячеек subs_fields/calc_thrust, с ВЫНУЖДЕННЫМ решением M4 X = B_ext*je0:

        import sys; sys.path.insert(0, '/mnt/raid/articles/sagemath_docker_build/work/mendrive_llg/MenDriveCpp/standalone')
        import export_qnm_profiles_sage as E
        F = E.driven_fields_from_notebook(globals(), omega_drive=1.7988142094590323e11, je0=1e6, source_fix=True)
        E.export_notebook_fields(globals(), F, 'out/qnm_drive_je_y_left_linear/nb_profile.csv', je0=1e6)

    source_fix=True делит свободный член B_ext на 2 (правильный источник 2*pi вместо 4*pi, см. NOTE §4.5);
    source_fix=False -- как сейчас в ноутбуке (амплитуды вдвое больше). Поля ВНУТРИ левой стенки в ноутбуке
    без секулярного члена (неточны); в зазоре и в правой стенке при source_fix=True они точные.

(2) ИЗ НОУТБУКА готовыми Ey_*_d / Hz_*_d (любой режим), через linear1d.export_sage_profiles:

        E.export_notebook_fields(globals(), None, 'nb_profile.csv', je0=1e6)

(3) БЕЗ Sage: точное решение с параметрами ноутбука (секулярный источник, полупространства):

        python3 export_qnm_profiles_sage.py --omega 1.7988142094590323e11 --out nb_exact.csv [--kz0] [--notebook-source]

    --kz0: k_z = 0 (как в 1D ./mendrive); по умолчанию k_z = pi/L_z, L_z = 100 см (ноутбук).
    --notebook-source: профиль как в ноутбуке (источник 4*pi, без секулярного члена) -- для сравнения.
Нужен только python3 (cmath); для (1),(2) -- ядро SageMath ноутбука.
"""
import cmath
import math
import os
import sys

C_CGS = 29979245800.0
# параметры MenDrive_qnm.ipynb (ячейки 6, 10, 22; вывод ячейки 420)
NB = dict(a=0.5, Lz=100.0, eps_l=1.0, mu_l=1.0, se_l=5.625e17, sm_l=0.0,
          eps_r=16.0, mu_r=100.0, se_r=9.9e12, sm_r=6.48e7, h_l=0.0001875, h_r=0.00625)


def _sqrt_im_pos(z):
    r = cmath.sqrt(z)
    return -r if (r.imag < 0 or (r.imag == 0 and r.real < 0)) else r


def _solve4(A, b):
    A = [row[:] for row in A]; b = b[:]; n = len(b)
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(A[i][k])); A[k], A[p] = A[p], A[k]; b[k], b[p] = b[p], b[k]
        for i in range(k + 1, n):
            f = A[i][k] / A[k][k]
            for j in range(k, n):
                A[i][j] -= f * A[k][j]
            b[i] -= f * b[k]
    x = [0j] * n
    for i in range(n - 1, -1, -1):
        x[i] = (b[i] - sum(A[i][j] * x[j] for j in range(i + 1, n))) / A[i][i]
    return x


def exact_profiles(omega, je0=1.0, p=NB, kz0=False, notebook_source=False):
    """Точное решение (СГС): u = E_y, u'' + K^2 u = -(4 pi i w mu'/c^2) j, j = je0*exp(-i K_l (x+a)).
    Возвращает функцию x -> (E_y, H_z) [СГС] и словарь с K."""
    c = C_CGS; w = float(omega); a = p['a']
    kz = 0.0 if kz0 else math.pi / p['Lz']
    el = p['eps_l'] + 4j * math.pi * p['se_l'] / w; ml = p['mu_l'] + 4j * math.pi * p['sm_l'] / w
    er = p['eps_r'] + 4j * math.pi * p['se_r'] / w; mr = p['mu_r'] + 4j * math.pi * p['sm_r'] / w
    Kv = cmath.sqrt(w * w / c ** 2 - kz ** 2)
    Kl = _sqrt_im_pos(el * ml * w * w / c ** 2 - kz ** 2); Kr = _sqrt_im_pos(er * mr * w * w / c ** 2 - kz ** 2)
    if notebook_source:
        # ноутбук: H_z однородно, E_y = (c mu w K B + 4 pi i mu w je0) e / D  (без секулярного члена)
        D = c ** 2 * kz ** 2 - el * ml * w * w
        up0 = 4j * math.pi * ml * w * je0 / D              # слагаемое E_y (ячейки 110/232), ~ exp(-iK(x+a))
        dup0 = -1j * Kl * up0
        up = lambda s: up0 * cmath.exp(-1j * Kl * s)
        dup = lambda s: dup0 * cmath.exp(-1j * Kl * s)
    else:
        g = -(4j * math.pi * w * ml / c ** 2) * je0
        al = g / (-2j * Kl)
        up = lambda s: al * s * cmath.exp(-1j * Kl * s)
        dup = lambda s: al * cmath.exp(-1j * Kl * s) * (1 - 1j * Kl * s)
    em, ep = cmath.exp(1j * Kv * a), cmath.exp(-1j * Kv * a)
    A = [[1, -em, -ep, 0], [-1j * Kl / ml, 1j * Kv * em, -1j * Kv * ep, 0],
         [0, ep, em, -1], [0, -1j * Kv * ep, 1j * Kv * em, -1j * Kr / mr]]
    b = [-up(0.0), -dup(0.0) / ml, 0, 0]
    if notebook_source:
        # в ноутбуке сшиваются E_y и H_z; H_z слева однородна: H_z = c u_h'/(i w mu) -> добавка источника в u' не входит
        b = [-up(0.0), 0, 0, 0]
    Cl, Pp, Qq, Cr = _solve4(A, b)

    def F(x):
        if x < -a:
            s = x + a
            u = Cl * cmath.exp(-1j * Kl * s) + up(s)
            du = -1j * Kl * Cl * cmath.exp(-1j * Kl * s) + (0 if notebook_source else dup(s))
            mz = ml
        elif x <= a:
            u = Pp * cmath.exp(-1j * Kv * x) + Qq * cmath.exp(1j * Kv * x)
            du = -1j * Kv * Pp * cmath.exp(-1j * Kv * x) + 1j * Kv * Qq * cmath.exp(1j * Kv * x); mz = 1.0
        else:
            u = Cr * cmath.exp(1j * Kr * (x - a)); du = 1j * Kr * u; mz = mr
        return u, c * du / (1j * w * mz)
    return F, dict(Kv=Kv, Kl=Kl, Kr=Kr, kz=kz)


def write_csv(path, xs, prof, omega_code, note=''):
    """Формат linear1d.write_profiles_csv: n,omega_n,x,x_cm,<F>_re,<F>_im,<F>_abs,<F>_arg."""
    names = list(prof.keys())
    os.makedirs(os.path.dirname(os.path.abspath(path)), exist_ok=True)
    with open(path, 'w', encoding='utf-8') as f:
        for line in note.splitlines():
            f.write('# ' + line + '\n')
        head = ['n', 'omega_n', 'x', 'x_cm']
        for nm in names:
            head += [nm + '_re', nm + '_im', nm + '_abs', nm + '_arg']
        f.write(','.join(head) + '\n')
        for k, xv in enumerate(xs):
            row = ['1', repr(float(omega_code)), repr(float(xv)), repr(float(xv))]
            for nm in names:
                v = complex(prof[nm][k])
                row += [repr(v.real), repr(v.imag), repr(abs(v)), repr(cmath.phase(v))]
            f.write(','.join(row) + '\n')
    return path


def grid(p=NB, n=400):
    a, hl, hr = p['a'], p['h_l'], p['h_r']
    xs = [-(a + hl) + hl * k / (n - 1) for k in range(n)]
    xs += [-a + 2 * a * k / (n - 1) for k in range(1, n)]
    xs += [a + hr * k / (n - 1) for k in range(1, n)]
    return xs


def export_exact(path, omega, je0=1.0, kz0=False, notebook_source=False, n=400):
    F, info = exact_profiles(omega, je0, kz0=kz0, notebook_source=notebook_source)
    xs = grid(n=n); j0_code = je0 / C_CGS                       # L0 = 1 см
    Ey, Hz = [], []
    for x in xs:
        e, h = F(x); Ey.append(e / j0_code); Hz.append(h / j0_code)
    note = (f'MenDrive_qnm.ipynb (driving_mode_1), {"источник как в ноутбуке (4 pi, без секулярного члена)" if notebook_source else "точное решение (секулярный источник)"}; '
            f'omega={omega!r} рад/с, omega_code={omega / C_CGS!r}; k_z={info["kz"]!r}; K_l={info["Kl"]}, K_r={info["Kr"]}\n'
            f'normalize=source: F/j0_code, j0_code = je0*L0/c, L0 = 1 см; time_convention=minus; x в см')
    return write_csv(path, xs, {'Ey': Ey, 'Hz': Hz}, omega / C_CGS, note)


# ---------------------------------------------------------------- Sage-часть (вызывать из ядра ноутбука)
def driven_fields_from_notebook(ns, omega_drive, je0=1e6, source_fix=True):
    """Повторяет шаги 6-9 calc_thrust (ячейка 396), но решает M4 X = B_ext*je0 (а не SVD однородной системы).
    ns = globals() ноутбука. Возвращает кортеж полей subs_fields(...)."""
    g = ns
    om, om_re, om_im, j_e, j_m = g['omega'], g['omega_re'], g['omega_im'], g['je0_ext'], g['jm0_ext']
    sol = [om_re == omega_drive, om_im == 0, om == omega_drive, j_e == je0, j_m == 0]
    dv = g['digit_values']
    K_v_d = [K.subs(g['k_subs']).subs(g['kappa_vacuum_sol']).subs(sol).subs(dv) for K in g['K_vacuum_subs']]
    K_l_d = g['fix_K_sign']([K.subs(g['k_subs']).subs(g['kappa_vacuum_sol']).subs(sol).subs(dv) for K in g['K_left_conductor_subs']])
    K_r_d = g['fix_K_sign']([K.subs(g['k_subs']).subs(g['kappa_vacuum_sol']).subs(sol).subs(dv) for K in g['K_right_conductor_subs']])
    def num(E):
        return E.subs(K_v_d).subs(K_l_d).subs(K_r_d).subs(g['kappa_vacuum_sol']).subs(g['epsilon_mu_subs']) \
                .subs(g['k_subs']).subs(dv).subs(sol).subs(g['sign_K_1_subs'])
    from sage.all import ComplexField, matrix, vector
    CF = ComplexField(128)
    M = matrix(CF, num(g['M4_K']).n(prec=128))
    B = vector(CF, [CF(num(b).n(prec=128)) for b in g['B4_K']])
    if source_fix:
        B = B / 2
    X = M.solve_right(B)
    res = (M * X - B).norm()
    print(f'[driven] |M X - B| = {float(res):.3e},  X = {[complex(v) for v in X]}')
    coeffs = [g['vars4'][i] == X[i] for i in range(len(g['vars4']))]
    return g['subs_fields'](sol, K_l_d, K_v_d, K_r_d, coeffs, dv)


def export_notebook_fields(ns, fields, path, je0=1e6, n_per_segment=400):
    """CSV из полей ноутбука. fields = кортеж subs_fields(...) или None (взять Ey_*_d, Hz_*_d из ns)."""
    here = os.path.dirname(os.path.abspath(__file__))
    if here not in sys.path:
        sys.path.insert(0, here)
    from linear1d import export_sage_profiles
    if fields is not None:
        (Ex_l, Ex_r, Ex_v, Ey_l, Ey_r, Ey_v, Ez_l, Ez_r, Ez_v, Dx_l, Dx_r, Dx_v, Dy_l, Dy_r, Dy_v, Dz_l, Dz_r, Dz_v,
         Hx_l, Hx_r, Hx_v, Hy_l, Hy_r, Hy_v, Hz_l, Hz_r, Hz_v) = fields[:27]
    else:
        Ey_l, Ey_v, Ey_r = ns['Ey_left_conductor_d'], ns['Ey_vacuum_d'], ns['Ey_right_conductor_d']
        Hz_l, Hz_v, Hz_r = ns['Hz_left_conductor_d'], ns['Hz_vacuum_d'], ns['Hz_right_conductor_d']
    A = float(ns.get('A', NB['a'])); hl = float(ns.get('h_conductor_l', NB['h_l'])); hr = float(ns.get('h_conductor_r', NB['h_r']))
    return export_sage_profiles(path, {'Ey': (Ey_l, Ey_v, Ey_r), 'Hz': (Hz_l, Hz_v, Hz_r)},
                                A=A, h_l=hl, h_r=hr, y0=0, n_per_segment=n_per_segment,
                                j0=float(je0) / C_CGS,
                                note='MenDrive_qnm.ipynb, je_y left, TE; нормировка F/(je0/c) = fields_harm_normalize=source при L0=1 см')


if __name__ == '__main__':
    import argparse
    ap = argparse.ArgumentParser(description='Профили Ey, Hz линейной модели MenDrive_qnm.ipynb -> CSV для plot.py --compare')
    ap.add_argument('--omega', type=float, default=1.7988142094590323e11, help='рад/с (по умолчанию Re omega_QNM)')
    ap.add_argument('--out', default='nb_exact_profile.csv')
    ap.add_argument('--kz0', action='store_true', help='k_z = 0 (как в 1D ./mendrive)')
    ap.add_argument('--notebook-source', action='store_true', help='источник как в ноутбуке (4 pi, без секулярного члена)')
    ap.add_argument('--n', type=int, default=400)
    a_ = ap.parse_args()
    print('записано:', export_exact(a_.out, a_.omega, kz0=a_.kz0, notebook_source=a_.notebook_source, n=a_.n))
