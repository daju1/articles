#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
tests/test_fields_linear.py -- проверка НОВЫХ режимов возбуждения (// MODE:) и экспорта полей (// FIELDS:)
программы ./mendrive. Запуск из папки standalone:  python3 tests/test_fields_linear.py   (или make test-fields)

Для каждого из 8 режимов (je_y_left, je_z_left, jm_z_right, jm_y_right, je_y_right, je_z_right, jm_z_left,
jm_y_left) с ЛИНЕЙНЫМ ферритом (ferrite_model=Linear) проверяется:
  1) ./mendrive  vs  numpy-копия схемы fdtd1d_ref.py         -- одна и та же разностная схема, допуск 1e-8;
  2) ./mendrive  vs  точное решение линейной задачи linear1d.py -- ошибка дискретизации O(dx), допуск 10% при N=80;
  3) неактивная пара (TE<->TM) тождественно равна нулю (polarization=auto);
  4) дуальность TE<->TM: (Ez, Hy)[je_z_left] = (Ey, -Hz)[je_y_left] и т.п. (допуск 1e-10);
  5) снимок по фазе (fields_snap.csv) = Re[F~ exp(-i w t)] из fields_harm_n1.csv (конвенция времени, допуск 1e-3);
  6) старые имена magnetic_right/electric_left == новые jm_z_right/je_y_left (JA, итоги бит-в-бит).
Нужен только numpy.
"""
import math
import os
import subprocess
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
sys.path.insert(0, ROOT)
import linear1d as L          # noqa: E402
import fdtd1d_ref as FD       # noqa: E402

EXE = os.path.join(ROOT, 'mendrive')
MODES = ['je_y_left', 'je_z_left', 'jm_z_right', 'jm_y_right', 'je_y_right', 'je_z_right', 'jm_z_left', 'jm_y_left']
MAT = dict(sigma_e_left=3.0, sigma_e_ferrite=0.2, eps_left=1.0, eps_ferrite=2.0)
LIN = dict(mu_yy_lin=3.0, mu_zz_lin=3.0, sigma_m_leak=0.3)
GEOM = dict(N=80, a=0.2, h_l=0.2, h_r=0.2, dt_frac=0.2)
SRC = dict(skin_depth=0.06, src_k_phase=5.0, src_sigma_fix=1)
RUN = dict(omega0=3.0, amp=1.0, ramp_periods=2.0, n_periods_total=40, record_from_period=10)
HARM_PERIODS = 10


def run_mendrive(mode, out_dir, extra=()):
    kv = dict(GEOM); kv.update(MAT); kv.update(LIN); kv.update(SRC); kv.update(RUN)
    kv.update(excitation_mode=mode, ferrite_model='Linear', polarization='auto', fields_harm='1',
              fields_harm_periods=HARM_PERIODS, ts_every=0, fft_bins=0, fft_harmonics=1, out_dir=out_dir)
    cmd = [EXE] + [f'{k}={v}' for k, v in kv.items()] + list(extra)
    r = subprocess.run(cmd, capture_output=True, text=True)
    if r.returncode != 0:
        print(r.stdout); print(r.stderr)
        raise RuntimeError('mendrive завершился с ошибкой')
    return r.stdout


def rel(a, b):
    return float(np.linalg.norm(a - b) / max(np.linalg.norm(b), 1e-300))


def flag(ok):
    return 'OK' if ok else 'FAIL'


def main():
    if not os.path.exists(EXE):
        print('сначала соберите: make'); return 1
    ok_all = True
    harm = {}
    for mode in MODES:
        out = os.path.join(HERE, 'out_fields_' + mode)
        txt = run_mendrive(mode, out)
        mode_line = [l for l in txt.splitlines() if l.startswith('[mode]')]
        d = L.read_profiles_csv(os.path.join(out, 'fields_harm_n1.csv'))
        harm[mode] = d
        pol = L.polarization_of(mode)
        En, Hn = ('Ez', 'Hy') if pol == 'TM' else ('Ey', 'Hz')
        Eo, Ho = ('Ey', 'Hz') if pol == 'TM' else ('Ez', 'Hy')
        E = L.complex_column(d, En); H = L.complex_column(d, Hn)
        # 1) numpy-копия схемы
        sim = FD.FDTD1D(excitation_mode=mode, polarization='auto', **GEOM, **MAT, **LIN, **SRC)
        ref = FD.run_harmonics(sim, RUN['omega0'], RUN['n_periods_total'], HARM_PERIODS, amp=RUN['amp'],
                               ramp_periods=RUN['ramp_periods'])[1]
        r1 = max(rel(E, ref[En]), rel(H, ref[Hn + '_A']))
        # 2) точное решение
        f = L.solve_linear(RUN['omega0'], GEOM['a'], GEOM['h_l'], GEOM['h_r'], mode, amp=RUN['amp'],
                           skin_depth=SRC['skin_depth'], src_k_phase=SRC['src_k_phase'],
                           mu_ferrite_yy=LIN['mu_yy_lin'], mu_ferrite_zz=LIN['mu_zz_lin'],
                           sigma_m_ferrite=LIN['sigma_m_leak'], **MAT)
        F = f(d['x'])
        r2 = max(rel(E, F[En]), rel(H, F[Hn]))
        # 3) неактивная пара
        r3 = float(np.max(np.abs(L.complex_column(d, Eo))) + np.max(np.abs(L.complex_column(d, Ho))))
        ok = (r1 <= 1e-8) and (r2 <= 0.10) and (r3 == 0.0)
        ok_all &= ok
        print(f'== {mode:11s} ({pol}): vs numpy-схема {r1:.2e} [{flag(r1 <= 1e-8)}], vs точное {r2:.2e} '
              f'[{flag(r2 <= 0.10)}], неактивная пара {r3:.1e} [{flag(r3 == 0.0)}]')
        if mode_line:
            print('   ' + mode_line[0])
    # 4) дуальность: (Ez,Hy)_TM = (Ey,-Hz)_TE  для je_*,  (Ez,Hy)_TM(jm_y) = (-Ey, Hz)_TE(jm_z)
    for te, tm, sE, sH in (('je_y_left', 'je_z_left', 1, -1), ('je_y_right', 'je_z_right', 1, -1),
                           ('jm_z_right', 'jm_y_right', -1, 1), ('jm_z_left', 'jm_y_left', -1, 1)):
        a, b = harm[te], harm[tm]
        rE = rel(L.complex_column(b, 'Ez'), sE * L.complex_column(a, 'Ey'))
        rH = rel(L.complex_column(b, 'Hy'), sH * L.complex_column(a, 'Hz'))
        okd = max(rE, rH) <= 1e-10
        ok_all &= okd
        print(f'== дуальность {te} <-> {tm}: {max(rE, rH):.2e} [{flag(okd)}]')
    # 5) снимки по фазам vs гармоника (конвенция exp(-i w t))
    out = os.path.join(HERE, 'out_fields_snap')
    run_mendrive('je_z_left', out, extra=['fields_snap_phases=0,0.25,0.5,0.75', 'fields_snap_every_periods=0'])
    d = L.read_profiles_csv(os.path.join(out, 'fields_harm_n1.csv'))
    sn = L.read_profiles_csv(os.path.join(out, 'fields_snap.csv'))
    w = RUN['omega0']; dt = None
    with open(os.path.join(out, 'fields_snap.csv'), encoding='utf-8') as fh:
        head = fh.readline()
    dt = float(head.split('dt=')[1].split()[0])
    worst = 0.0
    for ph in sorted(set(sn['phase_target'])):
        m = sn['phase_target'] == ph
        t = sn['t'][m][0]
        Ez_s = sn['Ez'][m]; Hy_s = sn['Hy'][m]
        Ez_h = (L.complex_column(d, 'Ez') * np.exp(-1j * w * t)).real
        Hy_h = (L.complex_column(d, 'Hy') * np.exp(-1j * w * (t + 0.5 * dt))).real   # H в момент t+dt/2
        e = max(np.max(np.abs(Ez_s - Ez_h)) / np.max(np.abs(Ez_h)), np.max(np.abs(Hy_s - Hy_h)) / np.max(np.abs(Hy_h)))
        worst = max(worst, e)
    oks = worst <= 1e-3
    ok_all &= oks
    print(f'== снимки fields_snap.csv vs Re[F~ exp(-i w t)]: макс. отн. откл. {worst:.2e} [{flag(oks)}]')
    # 6) старые имена == новые
    vals = {}
    for name in ('magnetic_right', 'jm_z_right', 'electric_left', 'je_y_left'):
        o = os.path.join(HERE, 'out_alias_' + name)
        r = subprocess.run([EXE, 'N=40', 'ferrite_model=JA', f'excitation_mode={name}', 'omega0=3.0',
                            'n_periods_total=6', 'record_from_period=2', 'ts_every=0', 'fft_bins=0',
                            f'out_dir={o}'], capture_output=True, text=True)
        if r.returncode != 0:
            print(r.stderr); return 1
        with open(os.path.join(o, 'summary.txt'), encoding='utf-8') as fh:
            vals[name] = [l for l in fh if l.split('=')[0] in ('dTxx_avg', 'P_src_avg', 'Pdiss_avg', 'hysteresis_area')]
    oka = vals['magnetic_right'] == vals['jm_z_right'] and vals['electric_left'] == vals['je_y_left']
    ok_all &= oka
    print(f'== legacy-имена == новые имена: [{flag(oka)}]')
    print('\nИТОГ:', 'все проверки OK' if ok_all else 'ЕСТЬ FAIL')
    return 0 if ok_all else 1


if __name__ == '__main__':
    sys.exit(main())
