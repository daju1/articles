"""
Обёртка ctypes для libmendrive.so -- загружается в любой Jupyter Notebook
без установки pybind11/cmake/pip-пакетов.

Использование:
    import ctypes, numpy as np
    lib = ctypes.CDLL('./libmendrive.so')
    # ... (см. настройку argtypes ниже, либо импортировать этот файл целиком) ...
    sim = MenDriveCpp(N=20, ferrite_model='Hybrid', bias_orientation='x', H0_bias=4.0)
    res = sim.run(omega0=7.45, n_periods=80, record_from_period=15)
"""
import ctypes
import numpy as np

lib = ctypes.CDLL('./libmendrive.so')  # укажите актуальный путь к .so

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

DBLP = ctypes.POINTER(ctypes.c_double)
lib.mendrive_run.argtypes = [
    ctypes.c_void_p, ctypes.c_double, ctypes.c_int, ctypes.c_int,
    ctypes.c_double, ctypes.c_double, ctypes.c_int,
    DBLP, DBLP, DBLP, DBLP, DBLP, DBLP, DBLP, DBLP, DBLP, DBLP,
    ctypes.POINTER(ctypes.c_int)]
lib.mendrive_run.restype = ctypes.c_int

BIAS_CODE = {'none': 0, 'x': 1, 'y': 2, 'z': 3, 'parallel': 4}
FERRITE_CODE = {'JA': 0, 'LLG': 1, 'Hybrid': 2, 'Preisach': 3}
EXC_CODE = {'magnetic_right': 0, 'electric_left': 1}


class MenDriveCpp:
    """Python-обёртка над C++ ядром резонатора МенДрайв. Интерфейс совпадает
    с оригинальным MenDriveFDTD -- run() возвращает тот же dict с ключами
    t, dTxx, Hn, Mn, HzJA, MzJA, P, Pdiss, T, dt, blew_up."""

    def __init__(self, N, a=0.2, h_l=0.2, h_r=0.2, dt_frac=0.1,
                 excitation_mode='magnetic_right', ferrite_model='JA',
                 bias_orientation='none', H0_bias=0.0,
                 sigma_m_leak=3.0, sigma_e_left=3.0, sigma_e_ferrite=0.0,
                 Ms=1.0, a_JA=0.3, alpha_JA=0.001, k_JA=0.15, c_JA=0.15,
                 n_sub=4, n_fp=3, n_newton=6,
                 Ms_llg=0.3, gamma_llg=1.0, alpha_llg=0.1, skin_depth=0.06,
                 n_hyst_pr=201, Hc_mean_pr=1.0, Hc_sigma_pr=0.3, Hb_sigma_pr=0.5,
                 sigma_m_leak_pr=3.0):
        self.handle = lib.mendrive_create(
            N, a, h_l, h_r, dt_frac,
            EXC_CODE[excitation_mode], FERRITE_CODE[ferrite_model],
            BIAS_CODE[bias_orientation], H0_bias,
            sigma_m_leak, sigma_e_left, sigma_e_ferrite,
            Ms, a_JA, alpha_JA, k_JA, c_JA,
            n_sub, n_fp, n_newton,
            Ms_llg, gamma_llg, alpha_llg, skin_depth,
            n_hyst_pr, Hc_mean_pr, Hc_sigma_pr, Hb_sigma_pr, sigma_m_leak_pr)
        self.n_fer = lib.mendrive_get_n_fer(self.handle)
        self.dt = lib.mendrive_get_dt(self.handle)

    def __del__(self):
        if getattr(self, 'handle', None):
            lib.mendrive_destroy(self.handle)

    def run(self, omega0, n_periods, record_from_period, amp=1.0, ramp_periods=2, probe_idx=0):
        T = 2*np.pi/omega0
        nsteps = int(n_periods*T/self.dt)
        record_start = int(record_from_period*T/self.dt)
        max_rec = max(0, nsteps - record_start)
        if max_rec <= 0:
            return dict(t=np.array([]), dTxx=np.array([]), Hn=np.array([]), Mn=np.zeros((0, 3)),
                        HzJA=np.array([]), MzJA=np.array([]), P=np.array([]), Pdiss=np.array([]),
                        T=T, dt=self.dt, blew_up=False)
        out_t = np.zeros(max_rec); out_dTxx = np.zeros(max_rec)
        out_Hn = np.zeros(max_rec); out_P = np.zeros(max_rec)
        out_Mx = np.zeros(max_rec); out_My = np.zeros(max_rec); out_Mz = np.zeros(max_rec)
        out_HzJA = np.zeros(max_rec); out_MzJA = np.zeros(max_rec)
        out_Pdiss = np.zeros(max_rec)
        blew_up = ctypes.c_int(0)
        n_rec = lib.mendrive_run(
            self.handle, omega0, nsteps, record_start, amp, ramp_periods, probe_idx,
            out_t.ctypes.data_as(DBLP), out_dTxx.ctypes.data_as(DBLP),
            out_Hn.ctypes.data_as(DBLP), out_P.ctypes.data_as(DBLP),
            out_Mx.ctypes.data_as(DBLP), out_My.ctypes.data_as(DBLP), out_Mz.ctypes.data_as(DBLP),
            out_HzJA.ctypes.data_as(DBLP), out_MzJA.ctypes.data_as(DBLP),
            out_Pdiss.ctypes.data_as(DBLP), ctypes.byref(blew_up))
        Mn = np.stack([out_Mx[:n_rec], out_My[:n_rec], out_Mz[:n_rec]], axis=1)
        return dict(t=out_t[:n_rec], dTxx=out_dTxx[:n_rec], Hn=out_Hn[:n_rec], Mn=Mn,
                    HzJA=out_HzJA[:n_rec], MzJA=out_MzJA[:n_rec], P=out_P[:n_rec],
                    Pdiss=out_Pdiss[:n_rec],
                    T=T, dt=self.dt, blew_up=bool(blew_up.value))

    def force_per_power(self, res, last_frac=0.5,
                         length_unit_cm=1.0, area_cm2=1.0, time_unit_s=1.0):
        """Отношение силы (тензор Максвелла, dTxx) к затраченной мощности (P) в
        Н/кВт, усреднённое по последней доле last_frac записанной траектории.

        ВНИМАНИЕ (калибровка СГС -> СИ): dTxx возвращается ядром в тех же
        единицах, что и энергия поля (E^2+H^2)/8*pi, то есть как плотность
        энергии в коде (код-единицы длины/времени/поля). Чтобы получить
        физическую силу в динах, эту плотность нужно умножить на площадь
        поперечного сечения резонатора в см^2 (area_cm2); чтобы получить
        физическую мощность в эрг/с, P нужно домножить на характерный масштаб
        мощности вашей нормировки (time_unit_s, length_unit_cm).
        По умолчанию все три калибровочных множителя равны 1.0, то есть метод
        возвращает "код-единицы Н/кВт" -- ПОДСТАВЬТЕ фактические масштабы из
        вашей нормировки (см. MenDrive_qnm_theory.ipynb), прежде чем
        интерпретировать force_per_kW как физическую величину.
        """
        t, dTxx, P, dt, T = res['t'], res['dTxx'], res['P'], res['dt'], res['T']
        if len(t) == 0:
            return np.nan, np.nan
        spp = max(1, int(round(T / dt)))
        n_periods = max(1, len(t) // spp)
        n_use = max(1, int(n_periods * last_frac)) * spp
        n_use = min(n_use, len(t))
        dTxx_avg = float(np.mean(dTxx[-n_use:]))
        P_avg = float(np.mean(P[-n_use:]))
        force_per_power_code = dTxx_avg / P_avg if abs(P_avg) > 1e-300 else np.nan

        force_dyn = dTxx_avg * area_cm2                        # эрг/см^3 * см^2 = дин
        power_erg_s = P_avg * (length_unit_cm ** 2) / time_unit_s
        force_N = force_dyn * 1e-5                             # дин -> Н
        power_kW = power_erg_s * 1e-7 * 1e-3                   # эрг/с -> Вт -> кВт
        force_per_kW = force_N / power_kW if abs(power_kW) > 1e-300 else np.nan
        return force_per_power_code, force_per_kW

    def check_energy_balance(self, res, last_frac=0.5):
        """Сверка энергобаланса: сравнение мощности источника (P, из -j*E) с
        ПОЛНЫМ временным интегралом закона сохранения энергии Пойнтинга:
        Pdiss = sigma_e*E^2 (обе стенки, точно) + sigma_m_leak*H^2 (вязкий
        канал JA/Hybrid/Preisach, точно) + H*dM/dt (энергия, передаваемая
        намагниченности -- площадь петли гистерезиса для JA/Preisach/Hybrid,
        работа демпфирования Гильберта для LLG).
        ВАЖНО: при коротких/грубых прогонах внутри полного резонатора P и
        Pdiss могут не совпадать даже с верной формулой -- если резонатор ещё
        не вышел на установившийся периодический режим, либо если динамика M
        заходит в нефизичные минорные петли. Перед выводами по балансу
        убедитесь в сходимости на достаточно длинном прогоне.
        Возвращает (P_avg, Pdiss_avg, относительное расхождение)."""
        t, dt, T = res['t'], res['dt'], res['T']
        if len(t) == 0:
            return np.nan, np.nan, np.nan
        spp = max(1, int(round(T / dt)))
        n_periods = max(1, len(t) // spp)
        n_use = max(1, int(round(n_periods * last_frac))) * spp
        n_use = min(n_use, len(t))
        P_avg = float(np.mean(res['P'][-n_use:]))
        Pdiss_avg = float(np.mean(res['Pdiss'][-n_use:]))
        denom = max(abs(P_avg), abs(Pdiss_avg), 1e-300)
        rel_diff = abs(P_avg - Pdiss_avg) / denom
        return P_avg, Pdiss_avg, rel_diff
