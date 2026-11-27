#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
plot.py -- НЕОБЯЗАТЕЛЬНЫЙ просмотр результатов автономной программы mendrive.
Только читает CSV из папки результатов (расчёт от Python не зависит).

    python3 plot.py out_dir            # сохраняет PNG в out_dir/plots/
    python3 plot.py out_dir --show     # и показывает окна

Нужен matplotlib (numpy не обязателен).
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


def read_summary(path):
    d = {}
    if os.path.exists(path):
        for line in open(path, encoding='utf-8'):
            if '=' in line and not line.startswith('#'):
                k, v = line.strip().split('=', 1)
                d[k] = v
    return d


def main():
    args = [a for a in sys.argv[1:] if not a.startswith('--')]
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

    for p in saved:
        print('сохранено:', p)
    if '--show' in sys.argv:
        plt.show()


if __name__ == '__main__':
    main()
