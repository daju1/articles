#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Проверка автономной версии против исходного ctypes-пути (запускать ПОСЛЕ `make test`
или вручную: make && make lib && make reflib && python3 tests/compare_ctypes.py).

Что сравнивается для каждого тестового случая (короткие прогоны, малое N):
  1) ИСХОДНОЕ ядро ../mendrive_core.cpp (libmendrive_ref.*) через ctypes -- ровно
     как MenDriveCpp.run() + force_per_power() + check_energy_balance() +
     shoelace_area() из ноутбука;
  2) НОВОЕ ядро standalone/mendrive_core.cpp (libmendrive_new.*) через тот же
     старый extern "C" API -- поэлементное сравнение всех массивов записи;
  3) автономная программа ./mendrive (потоковые средние) -- сравнение итоговых
     чисел из out_dir/summary.txt;
  4) MODE: НОВЫЙ конструктор mendrive_create_mode (polarization=auto, т.е. обновляется только
     нужная пара TE/TM) -- для старых режимов обязан давать ТЕ ЖЕ массивы, что исходное ядро
     (допуск 0 / 1e-12), а новые имена jm_z_right / je_y_left в ./mendrive -- те же итоги,
     что magnetic_right / electric_left.
Нужен только numpy. Python здесь используется лишь для проверки.
Новые режимы (je_z_left, jm_y_right, ...) и экспорт полей проверяет tests/test_fields_linear.py.
"""
import ctypes
import os
import platform
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
EXT = 'dylib' if platform.system() == 'Darwin' else 'so'
BIAS_CODE = {'none': 0, 'x': 1, 'y': 2, 'z': 3, 'parallel': 4}
FERRITE_CODE = {'JA': 0, 'LLG': 1, 'Hybrid': 2, 'Preisach': 3, 'Preisach2D': 4}
EXC_CODE = {'magnetic_right': 0, 'electric_left': 1}
DBLP = ctypes.POINTER(ctypes.c_double)

DEFAULTS = dict(a=0.2, h_l=0.2, h_r=0.2, dt_frac=0.1, excitation_mode='magnetic_right',
                ferrite_model='JA', bias_orientation='x', H0_bias=4.0, sigma_m_leak=3.0,
                sigma_e_left=3.0, sigma_e_ferrite=0.0, Ms=1.0, a_JA=0.3, alpha_JA=0.001,
                k_JA=0.15, c_JA=0.15, n_sub=2, n_fp=2, n_newton=4, Ms_llg=0.3, gamma_llg=1.0,
                alpha_llg=0.1, skin_depth=0.06, n_hyst_pr=201, Hc_mean_pr=1.0, Hc_sigma_pr=0.3,
                Hb_sigma_pr=0.5, sigma_m_leak_pr=3.0)
ORDER = ['N', 'a', 'h_l', 'h_r', 'dt_frac', 'excitation_mode', 'ferrite_model', 'bias_orientation',
         'H0_bias', 'sigma_m_leak', 'sigma_e_left', 'sigma_e_ferrite', 'Ms', 'a_JA', 'alpha_JA', 'k_JA',
         'c_JA', 'n_sub', 'n_fp', 'n_newton', 'Ms_llg', 'gamma_llg', 'alpha_llg', 'skin_depth',
         'n_hyst_pr', 'Hc_mean_pr', 'Hc_sigma_pr', 'Hb_sigma_pr', 'sigma_m_leak_pr']


def load(path):
    lib = ctypes.CDLL(path)
    lib.mendrive_create.restype = ctypes.c_void_p
    lib.mendrive_create.argtypes = [
        ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_int, ctypes.c_int, ctypes.c_int, ctypes.c_double,
        ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_int, ctypes.c_int, ctypes.c_int,
        ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double,
        ctypes.c_int, ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]
    lib.mendrive_destroy.argtypes = [ctypes.c_void_p]
    lib.mendrive_get_n_fer.argtypes = [ctypes.c_void_p]
    lib.mendrive_get_n_fer.restype = ctypes.c_int
    lib.mendrive_get_dt.argtypes = [ctypes.c_void_p]
    lib.mendrive_get_dt.restype = ctypes.c_double
    lib.mendrive_run.argtypes = [
        ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
        ctypes.c_double, ctypes.c_double, ctypes.c_int,
        DBLP, DBLP, DBLP, DBLP, DBLP, DBLP, DBLP, DBLP, DBLP, DBLP, ctypes.POINTER(ctypes.c_int)]
    lib.mendrive_run.restype = ctypes.c_int
    if hasattr(lib, 'mendrive_create_mode'):   # MODE: только в новом ядре
        lib.mendrive_create_mode.restype = ctypes.c_void_p
        lib.mendrive_create_mode.argtypes = list(lib.mendrive_create.argtypes) + [
            ctypes.c_int, ctypes.c_int,                         # pr2d_algo, n_threads
            ctypes.c_int, ctypes.c_double, ctypes.c_int,        # polarization, src_k_phase, src_sigma_fix
            ctypes.c_double, ctypes.c_double, ctypes.c_double, ctypes.c_double]  # eps_left, eps_ferrite, mu_yy, mu_zz
        lib.mendrive_get_polarization.argtypes = [ctypes.c_void_p]
        lib.mendrive_get_polarization.restype = ctypes.c_int
    return lib


def run_ctypes(lib, p, omega0, n_periods, record_from_period, amp, ramp_periods=2.0, probe_idx=0,
               create_mode=False, polarization=0):
    """Точная копия MenDriveCpp.__init__ + MenDriveCpp.run из mendrive_ctypes.py.
    create_mode=True -- MODE: тот же прогон через новый mendrive_create_mode (polarization: 0 auto, 1 TE, 2 TM, 3 both)."""
    args = []
    for k in ORDER:
        v = p[k]
        if k == 'excitation_mode': v = EXC_CODE[v]
        elif k == 'ferrite_model': v = FERRITE_CODE[v]
        elif k == 'bias_orientation': v = BIAS_CODE[v]
        args.append(v)
    if create_mode:
        h = lib.mendrive_create_mode(*args, 0, 1, polarization, 0.0, 0, 1.0, 1.0, 1.0, 1.0)
    else:
        h = lib.mendrive_create(*args)
    assert h, 'mendrive_create вернул NULL'
    dt = lib.mendrive_get_dt(h)
    T = 2*np.pi/omega0
    nsteps = int(n_periods*T/dt)
    record_start = int(record_from_period*T/dt)
    max_rec = max(0, nsteps - record_start)
    outs = [np.zeros(max_rec) for _ in range(10)]
    blew = ctypes.c_int(0)
    n_rec = lib.mendrive_run(h, omega0, nsteps, record_start, amp, ramp_periods, probe_idx,
                             *[o.ctypes.data_as(DBLP) for o in outs], ctypes.byref(blew))
    lib.mendrive_destroy(h)
    names = ['t', 'dTxx', 'Hn', 'P', 'Mx', 'My', 'Mz', 'HzJA', 'MzJA', 'Pdiss']
    res = {k: o[:n_rec] for k, o in zip(names, outs)}
    res.update(T=T, dt=dt, blew_up=bool(blew.value))
    return res


def metrics(res, model, last_frac=0.5):
    """force_per_power / check_energy_balance / shoelace_area -- как в ноутбуке."""
    t, dt, T = res['t'], res['dt'], res['T']
    spp = max(1, int(round(T/dt)))
    n_periods = max(1, len(t)//spp)
    n_use = min(max(1, int(n_periods*last_frac))*spp, len(t))
    dTxx = float(np.mean(res['dTxx'][-n_use:]))
    P = float(np.mean(res['P'][-n_use:]))
    Pd = float(np.mean(res['Pdiss'][-n_use:]))
    out = dict(dTxx_avg=dTxx, P_avg=P, Pdiss_avg=Pd, force_per_power_code=dTxx/P,
               rel_dif=abs(P-Pd)/max(abs(P), abs(Pd), 1e-300))
    if model in ('JA', 'Preisach', 'Preisach2D'):
        H = res['Hn'][-spp:]; B = H + 4*np.pi*res['My'][-spp:]
        out['hysteresis_area'] = 0.5*abs(np.sum(H*np.roll(B, -1) - np.roll(H, -1)*B))
    elif model == 'Hybrid':
        H = res['HzJA'][-spp:]; B = H + 4*np.pi*res['MzJA'][-spp:]
        out['hysteresis_area'] = 0.5*abs(np.sum(H*np.roll(B, -1) - np.roll(H, -1)*B))
    return out


def run_standalone(p, omega0, n_periods, record_from_period, amp, out_dir, pr2d_algo):
    exe = os.path.join(ROOT, 'mendrive')
    cmd = [exe] + [f'{k}={p[k]}' for k in ORDER] + [
        f'omega0={omega0!r}', f'n_periods_total={n_periods}', f'record_from_period={record_from_period}',
        f'amp={amp}', 'ramp_periods=2.0', 'probe_idx=0', 'last_frac=0.5', 'ts_every=0',
        'fft_bins=0', 'fft_harmonics=2', f'pr2d_algo={pr2d_algo}', f'out_dir={out_dir}']
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout); print(r.stderr)
        raise RuntimeError('mendrive завершился с ошибкой')
    vals = {}
    with open(os.path.join(out_dir, 'summary.txt')) as f:
        for line in f:
            if '=' in line and not line.startswith('#'):
                k, v = line.strip().split('=', 1)
                try: vals[k] = float(v)
                except ValueError: vals[k] = v
    return vals


def reldiff(a, b):
    return abs(a-b)/max(abs(a), abs(b), 1e-300)


CASES = [
    # имя, параметры, omega0, периоды, запись с, amp, алгоритм 2D, допуск
    ('JA_magnetic', dict(N=40, ferrite_model='JA'), 3.0, 12, 4, 1.0, 'rows', 1e-9),
    ('Preisach_magnetic', dict(N=40, ferrite_model='Preisach', n_hyst_pr=201), 3.0, 12, 4, 1.0, 'rows', 1e-9),
    ('Preisach_electric', dict(N=40, ferrite_model='Preisach', n_hyst_pr=201,
                               excitation_mode='electric_left'), 16.9, 30, 10, 3.0, 'rows', 1e-9),
    ('Preisach2D_bits', dict(N=40, ferrite_model='Preisach2D', n_hyst_pr=21), 3.0, 12, 4, 5.0, 'bits', 1e-9),
    ('Preisach2D_rows', dict(N=40, ferrite_model='Preisach2D', n_hyst_pr=21), 3.0, 12, 4, 5.0, 'rows', 1e-6),
    ('LLG_magnetic', dict(N=40, ferrite_model='LLG'), 3.0, 12, 4, 1.0, 'rows', 1e-9),
    ('Hybrid_magnetic', dict(N=40, ferrite_model='Hybrid'), 3.0, 12, 4, 1.0, 'rows', 1e-9),
]


def main():
    ref = load(os.path.join(ROOT, f'libmendrive_ref.{EXT}'))
    new = load(os.path.join(ROOT, f'libmendrive_new.{EXT}'))
    ok_all = True
    for name, over, w, nper, rec, amp, algo, tol in CASES:
        p = dict(DEFAULTS); p.update(over)
        r_ref = run_ctypes(ref, p, w, nper, rec, amp)
        r_new = run_ctypes(new, p, w, nper, rec, amp)
        dmax = max(float(np.max(np.abs(r_ref[k] - r_new[k]))) if len(r_ref[k]) else 0.0
                   for k in ('t', 'dTxx', 'Hn', 'P', 'My', 'Pdiss'))
        m_ref = metrics(r_ref, p['ferrite_model'])
        out_dir = os.path.join(HERE, 'out_' + name)
        s = run_standalone(p, w, nper, rec, amp, out_dir, algo)
        print(f'== {name}: точек={len(r_ref["t"])}, blew_up={r_ref["blew_up"]}, '
              f'max|старое ядро - новое ядро| (ctypes) = {dmax:.3e}')
        ok = True
        for k in ('dTxx_avg', 'P_avg', 'Pdiss_avg', 'force_per_power_code', 'rel_dif', 'hysteresis_area'):
            if k not in m_ref: continue
            d = reldiff(m_ref[k], s[k])
            flag = 'OK' if d <= tol else 'РАСХОЖДЕНИЕ'
            ok &= d <= tol
            print(f'   {k:22s} ноутбук={m_ref[k]: .12e}  mendrive={s[k]: .12e}  отн.разн.={d:.2e}  {flag}')
        print(f'   P_src_avg={s["P_src_avg"]:.6e}  P_num={s["P_num_residual"]:.6e}  '
              f'rel(P_src+P_num, Pdiss)={s["rel_dif_src_plus_num"]:.3e}  пик RSS={s["peak_rss_mb"]:.1f} МБ')
        # MODE: новый конструктор с polarization=auto (обновляется только нужная пара) -- те же массивы
        if hasattr(new, 'mendrive_create_mode'):
            r_mode = run_ctypes(new, p, w, nper, rec, amp, create_mode=True, polarization=0)
            dm = max(float(np.max(np.abs(r_ref[k] - r_mode[k]))) if len(r_ref[k]) else 0.0
                     for k in ('t', 'dTxx', 'Hn', 'P', 'Mx', 'My', 'Mz', 'Pdiss'))
            okm = dm <= 1e-12 * max(1.0, float(np.max(np.abs(r_ref['P']))) if len(r_ref['P']) else 1.0)
            print(f'   MODE: mendrive_create_mode(polarization=auto): max|исходное - новое| = {dm:.3e}  '
                  f'{"OK" if okm else "РАСХОЖДЕНИЕ"}')
            ok &= okm
        # MODE: новые имена режимов в ./mendrive -> те же итоги
        alias = {'magnetic_right': 'jm_z_right', 'electric_left': 'je_y_left'}[p['excitation_mode']]
        p2 = dict(p); p2['excitation_mode'] = alias
        s2 = run_standalone(p2, w, nper, rec, amp, out_dir + '_alias', algo)
        da = max(reldiff(s[k], s2[k]) for k in ('dTxx_avg', 'P_avg', 'Pdiss_avg', 'P_src_avg'))
        print(f'   MODE: excitation_mode={alias} vs {p["excitation_mode"]}: макс. отн. разн. итогов = {da:.2e}  '
              f'{"OK" if da == 0 else "РАСХОЖДЕНИЕ"}')
        ok &= (da == 0)
        ok_all &= ok
    print('\nИТОГ:', 'все совпадения в пределах допуска' if ok_all else 'ЕСТЬ РАСХОЖДЕНИЯ')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
