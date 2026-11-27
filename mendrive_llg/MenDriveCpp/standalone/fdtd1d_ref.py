# -*- coding: utf-8 -*-
"""
fdtd1d_ref.py -- построчная numpy-копия ЛИНЕЙНОЙ части шага MenDriveSim::step_once (standalone/mendrive_core.cpp,
пометки // MODE:) и накопления гармоник main.cpp (// FIELDS:). Нужна, чтобы проверять схему, знаки источников,
конвенцию фазы и нормировку БЕЗ компилятора, а после сборки -- сравнивать с ./mendrive бит-в-бит (до ~1e-12).

Поддерживаются: все excitation_mode (legacy и новые), polarization TE/TM/both, феррит 'Linear'
(mu_yy_lin, mu_zz_lin, sigma_m_leak) и 'none'-подобный (mu=1), eps_left/eps_ferrite, src_k_phase, src_sigma_fix.
Нелинейные модели (JA/LLG/Preisach) здесь НЕ эмулируются (для них есть tests/compare_ctypes.py).
"""
import math

import numpy as np

from linear1d import EXCITATIONS, polarization_of


class FDTD1D:
    def __init__(self, N=160, a=0.2, h_l=0.2, h_r=0.2, dt_frac=0.1, excitation_mode='jm_z_right',
                 polarization='auto', sigma_e_left=3.0, sigma_e_ferrite=0.0, eps_left=1.0, eps_ferrite=1.0,
                 mu_yy_lin=1.0, mu_zz_lin=1.0, sigma_m_leak=0.0, skin_depth=0.06, src_k_phase=0.0,
                 src_sigma_fix=0):
        self.N = N; self.a = a
        x_min, x_max = -(a + h_l), (a + h_r)
        self.dx = (x_max - x_min) / N; self.dt = dt_frac * self.dx
        dx, dt = self.dx, self.dt
        self.xA = x_min + dx * np.arange(N + 1); self.xB = x_min + dx * (np.arange(N) + 0.5)
        leftA, rightA = self.xA < -a, self.xA > a
        leftB, rightB = self.xB < -a, self.xB > a
        self.rightB = rightB
        kind, comp, side = EXCITATIONS[excitation_mode]
        self.kind, self.comp, self.side = kind, comp, side
        pol = polarization_of(excitation_mode) if polarization == 'auto' else polarization
        self.upd_TE = pol in ('TE', 'both'); self.upd_TM = pol in ('TM', 'both')
        sig = np.where(leftA, sigma_e_left, np.where(rightA, sigma_e_ferrite, 0.0))
        eps = np.where(leftA, eps_left, np.where(rightA, eps_ferrite, 1.0))
        self.Ca_e = (eps - 2 * math.pi * sig * dt) / (eps + 2 * math.pi * sig * dt)
        self.Cb_e = (dt / dx) / (eps + 2 * math.pi * sig * dt)
        self.Cj_e = dt * 4 * math.pi / ((eps + 2 * math.pi * sig * dt) if src_sigma_fix else 1.0)
        # магнитные коэффициенты (феррит 'Linear'); вне феррита mu=1, sigma_m=0
        mu_y = np.where(rightB, mu_yy_lin, 1.0); mu_z = np.where(rightB, mu_zz_lin, 1.0)
        sm = np.where(rightB, sigma_m_leak, 0.0)
        self.Ca_hy = (mu_y - 2 * math.pi * sm * dt) / (mu_y + 2 * math.pi * sm * dt)
        self.Cb_hy = dt / (mu_y + 2 * math.pi * sm * dt)
        self.Ca_hz = (mu_z - 2 * math.pi * sm * dt) / (mu_z + 2 * math.pi * sm * dt)
        self.Cb_hz = dt / (mu_z + 2 * math.pi * sm * dt)
        # источник: огибающая exp(-xi/delta), фаза k*xi (xi -- глубина в стенке)
        on_A = (kind == 'e'); xg = self.xA if on_A else self.xB
        inwall = (xg > a) if side == 'right' else (xg < -a)
        xi = (xg - a) if side == 'right' else (-a - xg)
        env = np.where(inwall, np.exp(-xi / skin_depth), 0.0)
        self.src_c = env * np.cos(src_k_phase * xi); self.src_s = env * np.sin(src_k_phase * xi)
        self.k_phase = src_k_phase
        self.Ey = np.zeros(N + 1); self.Ez = np.zeros(N + 1); self.Hy = np.zeros(N); self.Hz = np.zeros(N)
        self.Hy_prev = np.zeros(N); self.Hz_prev = np.zeros(N)
        self.j = np.zeros(N + 1 if on_A else N)
        self.t = 0.0

    def step(self, omega0, amp, ramp_periods, T):
        dt, dx, N = self.dt, self.dx, self.N
        ramp = min(1.0, self.t / (ramp_periods * T))
        if self.k_phase == 0.0:
            s = amp * ramp * math.sin(omega0 * self.t)
            self.j = s * self.src_c
        else:
            sa = amp * ramp
            self.j = sa * (math.sin(omega0 * self.t) * self.src_c - math.cos(omega0 * self.t) * self.src_s)
        je_y = self.j if (self.kind, self.comp) == ('e', 'y') else 0.0
        je_z = self.j if (self.kind, self.comp) == ('e', 'z') else 0.0
        jm_y = self.j if (self.kind, self.comp) == ('m', 'y') else 0.0
        jm_z = self.j if (self.kind, self.comp) == ('m', 'z') else 0.0
        if self.upd_TE:
            Ey_new = np.zeros(N + 1)
            Ey_new[1:N] = (self.Ca_e[1:N] * self.Ey[1:N] - self.Cb_e[1:N] * (self.Hz[1:] - self.Hz[:-1])
                           - self.Cj_e[1:N] * (je_y[1:N] if np.ndim(je_y) else 0.0))
            self.Ey = Ey_new
        if self.upd_TM:
            Ez_new = np.zeros(N + 1)
            Ez_new[1:N] = (self.Ca_e[1:N] * self.Ez[1:N] + self.Cb_e[1:N] * (self.Hy[1:] - self.Hy[:-1])
                           - self.Cj_e[1:N] * (je_z[1:N] if np.ndim(je_z) else 0.0))
            self.Ez = Ez_new
        if self.upd_TE:
            rotE_z = (self.Ey[1:] - self.Ey[:-1]) / dx
            self.Hz_prev = self.Hz
            self.Hz = self.Ca_hz * self.Hz - self.Cb_hz * (rotE_z + 4 * math.pi * jm_z)
        if self.upd_TM:
            rotE_y = (self.Ez[1:] - self.Ez[:-1]) / dx
            self.Hy_prev = self.Hy
            self.Hy = self.Ca_hy * self.Hy + self.Cb_hy * (rotE_y - 4 * math.pi * jm_y)
        self.t += dt


def run_harmonics(sim, omega0, n_periods_total, harm_periods, amp=1.0, ramp_periods=2.0, harmonics=(1,),
                  time_convention='minus'):
    """
    Прогон + накопление гармоник на лету (как // FIELDS: в main.cpp):
        F~_n(x) = (c_n / T_rec) * sum_k F(x, tau_k) * exp(+s*i*n*omega*tau_k) * dt,   c_0 = 1, c_n = 2,
    s = +1 для time_convention='minus' (F = Re[F~ e^{-i w t}], как в ноутбуках), s = -1 для 'plus'.
    tau_k -- истинное время отсчёта: E в t, H в t + dt/2 (схема Йи), j в t - dt.
    Возвращает dict {n: {'Ey','Ez','Hy','Hz','j'}} на A-сетке (H интерполируется в узлы xA).
    """
    T = 2 * math.pi / omega0; dt = sim.dt
    nsteps = int(n_periods_total * T / dt)
    spp = max(1, int(round(T / dt)))
    n_acc = harm_periods * spp
    start = nsteps - n_acc
    s = 1.0 if time_convention == 'minus' else -1.0
    acc = {n: {k: 0j for k in ()} for n in harmonics}
    for n in harmonics:
        acc[n] = dict(Ey=np.zeros(sim.N + 1, complex), Ez=np.zeros(sim.N + 1, complex),
                      Hy=np.zeros(sim.N, complex), Hz=np.zeros(sim.N, complex), j=np.zeros(sim.j.size, complex))
    for step in range(nsteps):
        sim.step(omega0, amp, ramp_periods, T)
        if step >= start:
            t = sim.t
            for n in harmonics:
                zE = np.exp(1j * s * n * omega0 * t)
                zH = np.exp(1j * s * n * omega0 * (t + 0.5 * dt))
                zj = np.exp(1j * s * n * omega0 * (t - dt))
                A = acc[n]
                A['Ey'] += sim.Ey * zE; A['Ez'] += sim.Ez * zE
                A['Hy'] += sim.Hy * zH; A['Hz'] += sim.Hz * zH
                A['j'] += sim.j * zj
    T_rec = n_acc * dt
    out = {}
    for n in harmonics:
        c = (1.0 if n == 0 else 2.0) * dt / T_rec
        A = acc[n]
        res = {k: c * v for k, v in A.items()}
        for k in ('Hy', 'Hz'):
            h = res[k]; hA = np.empty(sim.N + 1, complex)
            hA[1:sim.N] = 0.5 * (h[1:] + h[:-1]); hA[0] = h[0]; hA[sim.N] = h[-1]
            res[k + '_A'] = hA
        out[n] = res
    return out
