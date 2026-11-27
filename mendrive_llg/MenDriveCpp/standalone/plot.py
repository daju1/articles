#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot.py -- НЕОБЯЗАТЕЛЬНЫЙ просмотр результатов автономной программы mendrive.
Только читает CSV из папки результатов (расчёт от Python не зависит).

    python3 plot.py out_dir            # сохраняет PNG в out_dir/plots/
    python3 plot.py out_dir --show     # и показывает окна

FIELDS (README §9): если в out_dir есть fields_harm_n*.csv / fields_snap.csv, дополнительно строятся
профили полей по x в стиле ноутбуков MenDrive_qnm*.ipynb (Re, Im, Abs, Arg/pi тангенциальных
компонент E и H по трём участкам: левая стенка, вакуум, правая стенка; TM -- зелёным, как в ноутбуках)
и снимки по фазам периода. Сравнение с линейной моделью:

    python3 plot.py out_dir --compare linear.csv [--harm 1] [--fit phase|complex|none] [--norm max|none]
                            [--xcol x|x_cm|xa] [--xlim lo,hi] [--comp Ey,Hz]

linear.csv -- CSV того же формата (linear1d.write_profiles_csv / export_sage_profiles / analytic_for_run):
столбцы x (и/или x_cm), <F>_re, <F>_im. Кривые нормируются по max|F|, подбирается общая фаза
(или комплексный множитель), печатается относительная невязка L2 = ||alpha*F_cpp - F_lin|| / ||F_lin||.

Нужен matplotlib (numpy не обязателен; для --compare и профилей нужен numpy).
"""
import csv
import math
import os
import sys

import matplotlib
if '--show' not in sys.argv:
    matplotlib.use('Agg')
import matplotlib.pyplot as plt


def read_csv(path):
    if not os.path.exists(path):
        return None
    with open(path, newline='') as f:
        rows = list(csv.reader(f))
    if len(rows) < 2:
        return None
    head = rows[0]
    cols = {h: [] for h in head}
    for r in rows[1:]:
        for h, v in zip(head, r):
            try:
                cols[h].append(float(v))
            except ValueError:
                cols[h].append(math.nan)
    return cols


def read_csv_hash(path):
    """CSV со строками-комментариями '#' (fields_*.csv, linear*.csv). Возвращает (cols, comments)."""
    if not os.path.exists(path):
        return None, []
    comments, rows = [], []
    with open(path, newline='', encoding='utf-8') as f:
        for line in f:
            if line.startswith('#'):
                comments.append(line[1:].strip()); continue
            if line.strip():
                rows.append(line.rstrip('\n').split(','))
    if len(rows) < 2:
        return None, comments
    head = rows[0]
    cols = {h: [] for h in head}
    for r in rows[1:]:
        for h, v in zip(head, r):
            try:
                cols[h].append(float(v))
            except ValueError:
                cols[h].append(math.nan)
    return cols, comments


def header_value(comments, key, default=None):
    for c in comments:
        for tok in c.replace(';', ' ').split():
            if tok.startswith(key + '='):
                try:
                    return float(tok.split('=', 1)[1])
                except ValueError:
                    return tok.split('=', 1)[1]
    return default


# --- стиль графиков полей как в ноутбуках MenDrive_qnm*.ipynb (ячейки 521-537) ---
E_TANG = ('Ey', 'Ez')
H_TANG = ('Hy', 'Hz')
TITLES = {
    'E': 'График %s тангенциальной компоненты электрического поля',
    'H': 'График %s тангенциальной компоненты напряженности магнитного поля H',
    'B': 'График %s компоненты магнитной индукции B',
    'M': 'График %s компоненты намагниченности M',
    'j': 'График %s тока источника',
}
PARTS = (('Re', lambda z: z.real), ('Im', lambda z: z.imag), ('Abs', abs), ('Arg', None))


def comp_color(name):
    """TM-пара (Ez, Hy, By, My, je_z, jm_y) -- зелёный, как color="green" в ноутбуках; TE -- синий."""
    return 'green' if name in ('Ez', 'Hy', 'By', 'My', 'je_z', 'jm_y') else 'tab:blue'


def _segments(x, region):
    """Разбиение на участки левая стенка / вакуум / правая стенка (как три plot() в ноутбуке)."""
    segs = []
    for r in (-1, 0, 1):
        idx = [k for k, v in enumerate(region) if v == r]
        if idx:
            lo, hi = idx[0], idx[-1]
            # стыкуем участки по общему узлу, чтобы кривая была непрерывной
            if r == 0 and lo > 0: lo -= 1
            if r == 1 and lo > 0: lo -= 1
            segs.append((lo, hi + 1))
    return segs or [(0, len(x))]


def plot_harm_file(path, pdir, title0, saved):
    import numpy as np
    d, comm = read_csv_hash(path)
    if not d:
        return
    n = int(d['n'][0]) if 'n' in d else 1
    x = np.array(d['x']); region = d.get('region', [0] * len(x))
    a = header_value(comm, 'a')
    groups = (('E', E_TANG), ('H', H_TANG), ('B', ('By', 'Bz')), ('M', ('Mx', 'My', 'Mz')),
              ('j', ('je_y', 'je_z', 'jm_y', 'jm_z')))
    vals = {}
    for g, names in groups:
        for nm in names:
            if nm + '_re' in d:
                z = np.array(d[nm + '_re']) + 1j * np.array(d[nm + '_im'])
                if np.nanmax(np.abs(z)) > 0:
                    vals[nm] = z
    for g, names in groups:
        present = [nm for nm in names if nm in vals]
        if not present:
            continue
        fig, ax = plt.subplots(2, 2, figsize=(13, 8.5))
        for (pname, fn), axk in zip(PARTS, ax.flat):
            for nm in present:
                z = vals[nm]
                y = (np.angle(z) / np.pi) if fn is None else fn(z)
                for (lo, hi) in _segments(x, region):
                    axk.plot(x[lo:hi], y[lo:hi], color=comp_color(nm), lw=1.2,
                             label=nm if (lo, hi) == _segments(x, region)[0] else None)
            if a is not None:
                for xv in (-a, a):
                    axk.axvline(xv, color='gray', ls=':', lw=0.8)
            axk.set_title(TITLES[g] % (pname + ('/pi' if fn is None else '')), fontsize=10)
            axk.set_xlabel('x'); axk.grid(True); axk.legend(fontsize=8)
        fig.suptitle(f'{title0}: гармоника n={n} ({", ".join(present)})')
        fig.tight_layout()
        p = os.path.join(pdir, f'fields_harm_n{n}_{g}.png'); fig.savefig(p, dpi=120); plt.close(fig); saved.append(p)
    # delta Arg H и E (ячейка 537 ноутбуков)
    pairs = [(h, e) for h, e in (('Hz', 'Ey'), ('Hy', 'Ez')) if h in vals and e in vals]
    if pairs:
        plt.figure(figsize=(10, 4.5))
        for h, e in pairs:
            dphi = (np.angle(vals[h]) - np.angle(vals[e])) / np.pi
            for (lo, hi) in _segments(x, region):
                plt.plot(x[lo:hi], dphi[lo:hi], color=comp_color(h), lw=1.2)
        if a is not None:
            for xv in (-a, a):
                plt.axvline(xv, color='gray', ls=':', lw=0.8)
        plt.title('График delta Arg тангенциальных компонент H и E (/pi)'); plt.xlabel('x'); plt.grid(True)
        plt.tight_layout()
        p = os.path.join(pdir, f'fields_harm_n{n}_dArg.png'); plt.savefig(p, dpi=120); plt.close(); saved.append(p)
    # энергия / Пойнтинг / Txx (n=0 -- средние за период)
    q = [nm for nm in ('w', 'Sx', 'Txx') if nm + '_re' in d and max(abs(v) for v in d[nm + '_abs']) > 0]
    if q:
        plt.figure(figsize=(10, 4.5))
        for nm in q:
            yv = np.array(d[nm + '_re']) if n == 0 else np.array(d[nm + '_abs'])
            plt.plot(x, yv, lw=1.2, label=nm + (' (среднее)' if n == 0 else ' |.|'))
        if a is not None:
            for xv in (-a, a):
                plt.axvline(xv, color='gray', ls=':', lw=0.8)
        plt.title(f'{title0}: w, S_x, T_xx, гармоника n={n}'); plt.xlabel('x'); plt.grid(True); plt.legend()
        plt.tight_layout()
        p = os.path.join(pdir, f'fields_harm_n{n}_energy.png'); plt.savefig(p, dpi=120); plt.close(); saved.append(p)


def plot_snapshots(path, pdir, title0, saved):
    import numpy as np
    d, comm = read_csv_hash(path)
    if not d:
        return
    a = header_value(comm, 'a')
    per = np.array(d['period']); last = np.nanmax(per)
    for nm in ('Ey', 'Ez', 'Hy', 'Hz', 'Bz', 'By', 'Sx', 'Txx'):
        if nm not in d or np.nanmax(np.abs(d[nm])) == 0:
            continue
        plt.figure(figsize=(10, 4.8))
        for ph in sorted(set(np.array(d['phase_target'])[per == last])):
            m = (per == last) & (np.array(d['phase_target']) == ph)
            plt.plot(np.array(d['x'])[m], np.array(d[nm])[m], lw=1.2, label=f'фаза {ph:g}')
        if a is not None:
            for xv in (-a, a):
                plt.axvline(xv, color='gray', ls=':', lw=0.8)
        plt.title(f'{title0}: снимки {nm}(x) по фазам периода {int(last)}'); plt.xlabel('x'); plt.grid(True)
        plt.legend(fontsize=8); plt.tight_layout()
        p = os.path.join(pdir, f'fields_snap_{nm}.png'); plt.savefig(p, dpi=120); plt.close(); saved.append(p)


def compare_profiles(x_lin, f_lin, x_nl, f_nl, norm='max', fit='phase', xlim=None):
    """Та же функция, что linear1d.compare_profiles (копия, чтобы plot.py не зависел от модуля)."""
    import numpy as np
    x_lin = np.asarray(x_lin, float); x_nl = np.asarray(x_nl, float)
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


def opt(name, default=None):
    if name in sys.argv:
        k = sys.argv.index(name)
        if k + 1 < len(sys.argv):
            return sys.argv[k + 1]
    return default


def run_compare(out, pdir, title0, saved):
    import numpy as np
    lin_path = opt('--compare'); n = int(opt('--harm', '1'))
    fit = opt('--fit', 'phase'); norm = opt('--norm', 'max'); xcol = opt('--xcol', 'x')
    xl = opt('--xlim'); xlim = tuple(float(v) for v in xl.split(',')) if xl else None
    comps = opt('--comp')
    nl, comm_nl = read_csv_hash(os.path.join(out, f'fields_harm_n{n}.csv'))
    li, comm_li = read_csv_hash(lin_path)
    if not nl or not li:
        print('сравнение: нет данных (fields_harm_n%d.csv или %s)' % (n, lin_path)); return
    def xs(d, comm):
        if xcol == 'xa':
            a = header_value(comm, 'a')
            if a is None:
                raise SystemExit('--xcol xa: в заголовке CSV нет a=...')
            return np.array(d['x']) / a
        return np.array(d[xcol] if xcol in d else d['x'])
    x_nl, x_li = xs(nl, comm_nl), xs(li, comm_li)
    names = comps.split(',') if comps else [k[:-3] for k in li if k.endswith('_re') and k[:-3] + '_re' in nl]
    lines = []
    for nm in names:
        if nm + '_re' not in nl or nm + '_re' not in li:
            continue
        zl = np.array(li[nm + '_re']) + 1j * np.array(li[nm + '_im'])
        zn = np.array(nl[nm + '_re']) + 1j * np.array(nl[nm + '_im'])
        if np.nanmax(np.abs(zl)) == 0 or np.nanmax(np.abs(zn)) == 0:
            continue
        r, alpha, xx, u, v = compare_profiles(x_li, zl, x_nl, zn, norm=norm, fit=fit, xlim=xlim)
        line = (f'{nm:5s}: rel_L2 = {r:.4e}   alpha = {abs(alpha):.6g}*exp(i*{np.angle(alpha) / np.pi:+.6f}*pi)'
                f'   ({"OK" if r < 0.05 else "> 5%"})')
        print('сравнение ' + line); lines.append(line)
        fig, ax = plt.subplots(1, 3, figsize=(15, 4.5))
        for axk, (pname, fn) in zip(ax, (('Abs', np.abs), ('Re', np.real), ('Im', np.imag))):
            axk.plot(xx, fn(u), color='k', lw=2.0, alpha=0.6, label='линейная модель')
            axk.plot(xx, fn(v), color=comp_color(nm), lw=1.1, ls='--', label='C++ (alpha*F)')
            axk.set_title(f'{pname} {nm}, n={n}'); axk.grid(True); axk.set_xlabel(xcol); axk.legend(fontsize=8)
        fig.suptitle(f'{title0}: {nm} -- линейная vs нелинейная, rel_L2={r:.3e} (norm={norm}, fit={fit})')
        fig.tight_layout()
        p = os.path.join(pdir, f'compare_n{n}_{nm}.png'); fig.savefig(p, dpi=120); plt.close(fig); saved.append(p)
    with open(os.path.join(out, f'compare_n{n}.txt'), 'w', encoding='utf-8') as f:
        f.write(f'# linear={lin_path} norm={norm} fit={fit} xcol={xcol} xlim={xlim}\n')
        f.write('\n'.join(lines) + '\n')


def read_summary(path):
    d = {}
    if os.path.exists(path):
        for line in open(path, encoding='utf-8'):
            if '=' in line and not line.startswith('#'):
                k, v = line.strip().split('=', 1)
                d[k] = v
    return d


def main():
    vopts = ('--compare', '--harm', '--fit', '--norm', '--xcol', '--xlim', '--comp')
    args = [a for k, a in enumerate(sys.argv[1:], 1) if not a.startswith('--') and sys.argv[k - 1] not in vopts]
    out = args[0] if args else 'out'
    pdir = os.path.join(out, 'plots')
    os.makedirs(pdir, exist_ok=True)
    s = read_summary(os.path.join(out, 'summary.txt'))
    title = f"{s.get('ferrite_model', '?')} {s.get('excitation_mode', '?')} amp={s.get('amp', '?')}"
    saved = []

    # 1) скан резонанса
    sw, sf = read_csv(os.path.join(out, 'scan_wide.csv')), read_csv(os.path.join(out, 'scan_fine.csv'))
    if sw or sf:
        plt.figure(figsize=(11, 5))
        for c, lab, st in ((sw, 'грубый скан', 'o--'), (sf, 'уточняющий скан', 'o-')):
            if c:
                pts = sorted(zip(c['omega'], c['resp']))
                plt.plot([p[0] for p in pts], [p[1] for p in pts], st, ms=3, label=lab)
        if 'omega_res' in s:
            plt.axvline(float(s['omega_res']), color='red', ls='--', label=f"omega0={float(s['omega_res']):.4g}")
        plt.xlabel('omega0'); plt.ylabel('амплитуда отклика Hz'); plt.title(title + ': скан резонанса')
        plt.grid(True); plt.legend(); plt.tight_layout()
        p = os.path.join(pdir, 'scan.png'); plt.savefig(p, dpi=120); saved.append(p)

    # 2) сходимость и баланс по периодам
    per = read_csv(os.path.join(out, 'periods.csv'))
    if per:
        b = per['block']
        fig, ax = plt.subplots(3, 1, figsize=(11, 10), sharex=True)
        ax[0].plot(b, per['cum_dTxx_over_P'], lw=1)
        ax[0].set_ylabel('накопл. <dTxx>/<P>'); ax[0].grid(True); ax[0].set_title(title + ': сходимость и баланс')
        ax[1].semilogy(b, [abs(x) for x in per['P_src']], label='|P_src| (источник)')
        ax[1].semilogy(b, [abs(x) for x in per['Pdiss']], label='|Pdiss|')
        ax[1].semilogy(b, [abs(x) for x in per['P_num']], label='|P_num| (невязка схемы)')
        ax[1].set_ylabel('мощность за период'); ax[1].grid(True); ax[1].legend()
        ax[2].plot(b, per['Hn_max'], label='max Hn'); ax[2].plot(b, per['Hn_min'], label='min Hn')
        ax[2].plot(b, per['dTxx'], label='<dTxx> за период')
        ax[2].set_xlabel('период записи'); ax[2].grid(True); ax[2].legend()
        plt.tight_layout()
        p = os.path.join(pdir, 'periods.png'); plt.savefig(p, dpi=120); saved.append(p)

    # 3) петля гистерезиса (последний период)
    lp = read_csv(os.path.join(out, 'loop.csv'))
    if lp:
        plt.figure(figsize=(8, 7))
        plt.plot(lp['H'], lp['B'], 'o-', ms=2, lw=1, color='tab:purple')
        plt.xlabel('H (Э)'); plt.ylabel('B (Гс)')
        plt.title(f"{title}: петля (площадь={s.get('hysteresis_area', '?')})")
        plt.grid(True); plt.tight_layout()
        p = os.path.join(pdir, 'loop.png'); plt.savefig(p, dpi=120); saved.append(p)

    # 4) спектр dTxx
    sp = read_csv(os.path.join(out, 'spectrum.csv'))
    if sp:
        plt.figure(figsize=(10, 5))
        plt.semilogy(sp['omega'], [max(x, 1e-300) for x in sp['abs']], lw=1)
        if 'omega_res' in s:
            plt.axvline(float(s['omega_res']), color='red', ls='--', alpha=0.6)
        plt.xlabel('omega'); plt.ylabel('|ДПФ(dTxx·окно Ханна)|'); plt.title(title + ': спектр dTxx')
        plt.grid(True); plt.tight_layout()
        p = os.path.join(pdir, 'spectrum.png'); plt.savefig(p, dpi=120); saved.append(p)

    # 5) хвост временного ряда
    ts = read_csv(os.path.join(out, 'timeseries.csv'))
    if ts:
        fig, ax = plt.subplots(3, 1, figsize=(11, 9), sharex=True)
        ax[0].plot(ts['t'], ts['Hn'], lw=0.8); ax[0].set_ylabel('Hn'); ax[0].grid(True)
        ax[1].plot(ts['t'], ts['dTxx'], lw=0.8); ax[1].set_ylabel('dTxx'); ax[1].grid(True)
        ax[2].plot(ts['t'], ts['P_src'], lw=0.8, label='P_src')
        ax[2].plot(ts['t'], ts['Pdiss'], lw=0.8, label='Pdiss')
        ax[2].set_xlabel('t'); ax[2].grid(True); ax[2].legend()
        ax[0].set_title(title + ': хвост временного ряда')
        plt.tight_layout()
        p = os.path.join(pdir, 'timeseries.png'); plt.savefig(p, dpi=120); saved.append(p)

    # 6) FIELDS: профили полей (гармоники), снимки по фазам, сравнение с линейной моделью
    import glob
    for hp in sorted(glob.glob(os.path.join(out, 'fields_harm_n*.csv'))):
        plot_harm_file(hp, pdir, title, saved)
    if os.path.exists(os.path.join(out, 'fields_snap.csv')):
        plot_snapshots(os.path.join(out, 'fields_snap.csv'), pdir, title, saved)
    if '--compare' in sys.argv:
        run_compare(out, pdir, title, saved)

    for p in saved:
        print('сохранено:', p)
    if '--show' in sys.argv:
        plt.show()


if __name__ == '__main__':
    main()
