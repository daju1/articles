#!/usr/bin/env bash
# qnm_sweep.sh -- сетка по N и по частоте для конфигурации MenDrive_qnm.ipynb (je_y_left, TE, Linear),
# затем сводная таблица: dTxx, P_src, P_diss, P_num, Н/кВт (суммарно и на правую стенку) против точного 1D-решения.
#
#   ./qnm_sweep.sh                                   # по умолчанию N_LIST и W_LIST ниже
#   N_LIST="16384 65536 131072" W_LIST="6.000198342077812" JOBS=3 ./qnm_sweep.sh
#   SUMMARY_ONLY=1 ./qnm_sweep.sh                    # только пересобрать таблицу по готовым прогонам
#   DRY=1 ./qnm_sweep.sh                             # показать команды и оценку стоимости, не запускать
#
# Нужно: собранный ./mendrive (make CXX=g++), python3 (только стандартная библиотека).
# Частоты (omega_code = omega_СГС/c, L0 = 1 см) -- корни QNM ноутбука (ячейка 441) и +-10% около основного:
#   2.944600605222343 (Re omega = 8.8277e10), 5.400178507870031 (0.9*1.7988e11), 6.000198342077812 (1.7988e11, ячейка 465),
#   6.600218176285594 (1.1*1.7988e11), 9.073479328679792 (2.7202e11).
# Частоты скана ячейки 479 (omega = 1e6..1e7 рад/с -> omega_code = 3.3e-5..3.3e-4) во FDTD недостижимы:
#   ~2*pi/(omega_code*dt) ~ 1e10..1e11 шагов на период.
# N = 80/160/320 -- контрольные: dx > h_l, источник je_y вне сетки (P_src = 0 -> nan). Скин-слой источника
#   delta_l = 3.76e-5 см разрешается при N >= ~65536 (delta_l/dx >= 2.5); основной прогон N = 131072.
# Стоимость ~ N^2: N=131072 -- ~3.6e10 обновлений ячеек за период (60 периодов ~1-2.5 ч/ядро), N=65536 -- в 4 раза меньше.
set -euo pipefail
cd "$(dirname "$0")"
MENDRIVE=${MENDRIVE:-./mendrive}
CFG=${CFG:-examples/qnm_drive_je_y_left_linear.cfg}
OUT=${OUT:-out/qnm_sweep}
PY=${PY:-python3}
JOBS=${JOBS:-1}
N_LIST=${N_LIST:-"80 160 320 16384 65536 131072"}
W_LIST=${W_LIST:-"2.944600605222343 5.400178507870031 6.000198342077812 6.600218176285594 9.073479328679792"}
EXTRA=${EXTRA:-}
mkdir -p "$OUT"

# K_l(omega) левой стенки (eps=1, mu=1, sigma_e из конфига) -> skin_depth = 1/Im K_l, src_k_phase = Re K_l
src_params() {
  "$PY" - "$1" "$CFG" <<'PYEOF'
import cmath, math, sys
w = float(sys.argv[1]); d = {}
for line in open(sys.argv[2], encoding='utf-8'):
    line = line.split('#', 1)[0].strip()
    if '=' in line:
        k, v = line.split('=', 1); d[k.strip()] = v.strip()
eps = float(d.get('eps_left', 1)); sig = float(d.get('sigma_e_left', 0))
K = cmath.sqrt(w * w * (eps + 4j * math.pi * sig / w))
if K.imag < 0: K = -K
print(repr(1.0 / K.imag), repr(K.real))
PYEOF
}

if [ "${SUMMARY_ONLY:-0}" != "1" ]; then
  cmds="$OUT/commands.txt"; : > "$cmds"
  for N in $N_LIST; do
    for W in $W_LIST; do
      read -r SD KP < <(src_params "$W")
      od="$OUT/N${N}_w${W}"
      if [ -f "$od/summary.txt" ]; then echo "[skip] $od (уже есть summary.txt)"; continue; fi
      est=$("$PY" -c "import math;N=$N;w=$W;dx=(1.0+1.875e-4+6.25e-3)/N;dt=0.5*dx;print('%.2e'%(60*2*math.pi/w/dt*N))")
      echo "[plan] N=$N omega0=$W skin_depth=$SD src_k_phase=$KP  (~$est обновлений ячеек)"
      echo "$MENDRIVE -c $CFG N=$N omega0=$W skin_depth=$SD src_k_phase=$KP out_dir=$od $EXTRA > $od.log 2>&1" >> "$cmds"
    done
  done
  if [ "${DRY:-0}" = "1" ]; then cat "$cmds"; exit 0; fi
  [ -x "$MENDRIVE" ] || { echo "нет $MENDRIVE: соберите  make CXX=g++"; exit 1; }
  if [ -s "$cmds" ]; then
    if [ "$JOBS" -gt 1 ]; then
      xargs -P "$JOBS" -I{} sh -c '{}' < "$cmds" || true
    else
      while IFS= read -r c; do echo "[run] $c"; sh -c "$c" || echo "[fail] $c"; done < "$cmds"
    fi
  fi
fi

"$PY" - "$OUT" <<'PYEOF'
import cmath, csv, glob, math, os, sys
C_CGS = 29979245800.0
def sqrt_im_pos(z):
    r = cmath.sqrt(z)
    return -r if (r.imag < 0 or (r.imag == 0 and r.real < 0)) else r
def read_kv(path):
    d = {}
    if os.path.exists(path):
        for line in open(path, encoding='utf-8'):
            line = line.split('#', 1)[0].strip()
            if '=' in line:
                k, v = line.split('=', 1); d[k.strip()] = v.strip()
    return d
def solve4(A, b):
    A = [row[:] for row in A]; b = b[:]; n = 4
    for k in range(n):
        p = max(range(k, n), key=lambda i: abs(A[i][k])); A[k], A[p] = A[p], A[k]; b[k], b[p] = b[p], b[k]
        for i in range(k + 1, n):
            f = A[i][k] / A[k][k]
            for j in range(k, n): A[i][j] -= f * A[k][j]
            b[i] -= f * b[k]
    x = [0j] * n
    for i in range(n - 1, -1, -1):
        x[i] = (b[i] - sum(A[i][j] * x[j] for j in range(i + 1, n))) / A[i][i]
    return x
def exact_TE(cfg):
    """Точное 1D TE-решение (k_z=0, c=1, полубесконечные стенки, источник exp(i Ks xi) в левой стенке, секулярный член при Ks=K_l)."""
    g = lambda k, d: float(cfg.get(k, d))
    w = g('omega0', 'nan'); a = g('a', 0.5); amp = g('amp', 1.0)
    el = g('eps_left', 1) + 4j * math.pi * g('sigma_e_left', 0) / w; ml = 1.0
    er = g('eps_ferrite', 1) + 4j * math.pi * g('sigma_e_ferrite', 0) / w
    mr = g('mu_zz_lin', 1) + 4j * math.pi * g('sigma_m_leak', 0) / w
    Kl = sqrt_im_pos(w * w * el * ml); Kr = sqrt_im_pos(w * w * er * mr); Kv = w
    Ks = g('src_k_phase', 0) + 1j / g('skin_depth', 0.06)
    j0 = 1j * amp
    gsrc = -(4j * math.pi * w * ml) * j0                       # u'' + K^2 u = gsrc*exp(-i Ks s)
    if abs(Ks - Kl) <= 1e-9 * abs(Kl):
        al = gsrc / (-2j * Kl); up0 = 0; dup0 = al                # u_p = al*s*exp(-iKl s)
        Pint = lambda Cl: -0.5 * (j0.conjugate() * (Cl / (2 * Kl.imag) - al / (4 * Kl.imag ** 2))).real
    else:
        be = gsrc / (Kl ** 2 - Ks ** 2); up0 = be; dup0 = -1j * Ks * be
        q = 1j * (Ks.conjugate() - Kl)                            # для интеграла P_src
        Pint = lambda Cl: -0.5 * (j0.conjugate() * (Cl / (-1j * (Kl - Ks.conjugate())) + be / (2 * Ks.imag))).real
    em, ep = cmath.exp(1j * Kv * a), cmath.exp(-1j * Kv * a)
    A = [[1, -em, -ep, 0], [-1j * Kl / ml, 1j * Kv * em, -1j * Kv * ep, 0], [0, ep, em, -1], [0, -1j * Kv * ep, 1j * Kv * em, -1j * Kr / mr]]
    b = [-up0, -dup0 / ml, 0, 0]
    Cl, Pp, Qq, Cr = solve4(A, b)
    def vac(x):
        u = Pp * cmath.exp(-1j * Kv * x) + Qq * cmath.exp(1j * Kv * x)
        du = -1j * Kv * Pp * cmath.exp(-1j * Kv * x) + 1j * Kv * Qq * cmath.exp(1j * Kv * x)
        return u, du / (1j * w)
    Ea, Ha = vac(a); Em, Hm = vac(-a)
    p_r = (abs(Ea) ** 2 + abs(Ha) ** 2) / (16 * math.pi); p_l = (abs(Em) ** 2 + abs(Hm) ** 2) / (16 * math.pi)
    return dict(P_src=Pint(Cl), p_right=p_r, p_left=p_l, Kl=Kl, Kr=Kr)
def n0_pressure(path, a):
    if not os.path.exists(path): return float('nan'), float('nan')
    rows = [r for r in csv.reader(l for l in open(path, encoding='utf-8') if not l.startswith('#'))]
    h = rows[0]; ix, ir, it = h.index('x'), h.index('region'), h.index('Txx_re')
    vac = [(float(r[ix]), float(r[it])) for r in rows[1:] if r and float(r[ir]) == 0]
    if not vac: return float('nan'), float('nan')
    return max(vac)[1], min(vac)[1]
def main(base):
    out = []
    for d in sorted(glob.glob(os.path.join(base, 'N*_w*'))):
        cfg = read_kv(os.path.join(d, 'config_used.txt')); s = read_kv(os.path.join(d, 'summary.txt'))
        if not s: continue
        f = lambda k: float(s.get(k, 'nan'))
        N = int(float(cfg.get('N', 'nan'))); w = float(cfg.get('omega0', 'nan')); a = float(cfg.get('a', 0.5))
        dT, Ps, Pd, Pn = f('dTxx_avg'), f('P_src_avg'), f('Pdiss_avg'), f('P_num_residual')
        pr, pl = n0_pressure(os.path.join(d, 'fields_harm_n0.csv'), a)
        ex = exact_TE(cfg)
        nkw = lambda F, P: 1e5 * F / (C_CGS * P) if P not in (0.0,) and P == P else float('nan')
        out.append(dict(N=N, omega_code=w, omega_cgs=w * C_CGS, dTxx=dT, P_src=Ps, P_diss=Pd, P_num=Pn,
                        rel_src=f('rel_dif_src'), rel_src_num=f('rel_dif_src_plus_num'),
                        p_right=pr, p_left=pl, NkW_thrust=nkw(dT, Ps), NkW_right=nkw(pr, Ps),
                        P_src_exact=ex['P_src'], NkW_right_exact=nkw(ex['p_right'], ex['P_src']),
                        P_src_ratio=Ps / ex['P_src'] if ex['P_src'] else float('nan'), blew_up=s.get('blew_up', '?')))
    if not out:
        print('нет готовых прогонов в', base); return
    keys = list(out[0].keys())
    with open(os.path.join(base, 'summary_table.csv'), 'w', encoding='utf-8') as fo:
        fo.write(','.join(keys) + '\n')
        for r in out: fo.write(','.join(repr(r[k]) for k in keys) + '\n')
    hdr = ['N', 'omega_code', 'dTxx', 'P_src', 'P_diss', 'P_num', 'rel_src', 'NkW_thrust', 'NkW_right', 'NkW_right_exact', 'P_src_ratio']
    line = '| ' + ' | '.join(hdr) + ' |'
    md = [line, '|' + '---|' * len(hdr)]
    for r in out:
        md.append('| ' + ' | '.join(f'{r[k]:.4g}' if isinstance(r[k], float) else str(r[k]) for k in hdr) + ' |')
    open(os.path.join(base, 'summary_table.md'), 'w', encoding='utf-8').write('\n'.join(md) + '\n')
    print('\n'.join(md))
    print('\nНазвания: dTxx = <p(a) - p(-a)> (код); NkW_thrust = 1e5*dTxx/(c*P_src) -- суммарная тяга (точно 0 для линейной задачи);')
    print('NkW_right = 1e5*<p(a)>/(c*P_src) -- сила на правую стенку; *_exact -- точное 1D-решение (k_z=0, полупространства); P_src_ratio = P_src/P_src_exact.')
    print('таблицы:', os.path.join(base, 'summary_table.csv'), os.path.join(base, 'summary_table.md'))
if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else 'out/qnm_sweep')
PYEOF
