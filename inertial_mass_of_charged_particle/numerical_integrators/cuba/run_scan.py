#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
run_scan.py: серия запусков calc_transverse_sphere_mass_film с разными параметрами.
Совместим с Python 3.5+ (без f-строк, без capture_output/text).

Пример:
  python3 run_scan.py --rho0 0.92,2,100 --delay 0,1 --fast 0,1 --jobs 4

python3 run_scan.py --rho0 0.92,2,100 --modes none,fermi_P,general_P --delay 0,1 --fast 1 --rot 0,1 --jobs 4 --out scan_rot01.csv

python3 run_scan.py --rho0 100 --modes fermi_O,film,film_inst --delay 1 --fast 1 --rot 1 \
    --v-c 0.01 --film-K 100 --out scan_film_calib.csv --logdir logs_film_calib

python3 run_scan.py --rho0 0.92,1.5,2,100 --modes none,fermi_O,fermi_P,film,film_inst,pde \
    --delay 1 --fast 1 --rot 1 --v-c 0.01 --film-K 100 --pde-n 24 --jobs 4 \
    --out scan_film_small.csv --logdir logs_film_small

python3 run_scan.py --rho0 0.92,1.5,2,3,5,10,100 --modes none,fermi_P \
    --delay 0 --fast 1 --rot 1 --v-c 0.01 --out scan_vc001.csv

for K in 25 400; do
  python3 run_scan.py --rho0 2 --modes film_inst --delay 1 --fast 1 --rot 1 \
      --lorentz lorentz --v-c 0.5 --film-K $K --out scan_K$K.csv --logdir logs_K$K
done

for vc in 0.4 0.5; do
  python3 run_scan.py --rho0 2,100 --modes fermi_O,fermi_P,film_inst --delay 1 --fast 1 --rot 1 \
      --lorentz lorentz --v-c $vc --film-K 400 --jobs 4 \
      --seed 12345,777 --maxeval 3000000 \
      --out scan_w1_vc$vc.csv --logdir logs_w1_vc$vc
done

for vc in 0.3 0.4 0.5; do
  python3 run_scan.py --rho0 2,100 --modes fermi_P,film_inst --delay 1 --fast 1 --rot 1 \
      --lorentz lorentz --v-c $vc --film-K 100 --jobs 4 \
      --seed 12345,777 --maxeval 3000000 \
      --out scan_fi_s_vc$vc.csv --logdir logs_fi_s_vc$vc
done

for vc in 0.3 0.4 0.5; do
  python3 run_scan.py --rho0 2,100 --modes fermi_P,film_inst --delay 1 --fast 1 --rot 1 \
      --lorentz lorentz --v-c $vc --film-K 100 --jobs 4 \
      --out scan_fi_fix_vc$vc.csv --logdir logs_fi_fix_vc$vc
done

for vc in 0.2 0.3 0.4 0.5; do
  python3 run_scan.py --rho0 2,100 --modes fermi_O,fermi_P,film,film_inst \
      --delay 1 --fast 1 --rot 1 --lorentz lorentz --v-c $vc --film-K 100 --jobs 4 \
      --out scan_filmrel_vc$vc.csv --logdir logs_filmrel_vc$vc
done

for vc in 0.1 0.2 0.3 0.4 0.5; do
  python3 run_scan.py --rho0 2,100 --modes none,fermi_P --delay 0,1 --fast 1 --rot 1 \
      --lorentz off,lorentz --v-c $vc --jobs 4 \
      --out scan_lor_vc$vc.csv --logdir logs_lor_vc$vc
done


for rho0 in 2,3,4,5,6,7; do
  python3 run_scan.py --rho0 $rho0 --modes none,fermi_P,fermi_O,general_P,general_O --delay 0,1 --fast 1 --rot 0,1 \
      --lorentz off,lorentz --v-c 0.01 --jobs 4 \
      --out scan_lor_rho0_$rho0.csv --logdir logs_lor_rho0_$rho0
done

"""
import argparse
import csv
import io
import itertools
import os
import re
import subprocess
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed

# Не зависеть от локали терминала (LANG=C на старых системах)
sys.stdout = io.TextIOWrapper(sys.stdout.buffer, encoding="utf-8",
                              errors="replace", line_buffering=True)
sys.stderr = io.TextIOWrapper(sys.stderr.buffer, encoding="utf-8",
                              errors="replace", line_buffering=True)

RESULT_NAMES = [
    "observed_ratio_1_rho", "observed_ratio_1_phi", "observed_ratio_1_z",
    "observed_ratio_2_rho", "observed_ratio_2_phi", "observed_ratio_2_z",
    "observed_ratio_rho",   "observed_ratio_phi",   "observed_ratio_z",
]

# режим поправки Ферми -> ключи командной строки
MODES = {
    "none":      [],
    "fermi_O":   ["--fermi-O=1"],
    "fermi_P":   ["--fermi=1"],
    "general_P": ["--fermi-general=1"],
    "general_O": ["--fermi-general-O=1"],
    "film":      ["--film=1"],
    "film_inst": ["--film=1", "--film-instant=1"],
    "pde":       ["--pde=1"],
}
NEEDS_DELAY = set(["film", "film_inst", "pde"])

LORENTZ = {
    "off":     [],
    "lorentz": ["--lorentz=1"],
    "general": ["--lorentz-general=1"],
}

RESULT_RE = re.compile(r"^RESULT:\s+(\S+)\s+(\S+)\s+\+-\s+(\S+)\s+p\s*=\s*(\S+)", re.M)
STATUS_RE = re.compile(r"^STATUS\s+fail=(-?\d+)", re.M)
PARAMS_RE = re.compile(r"^PARAMS .*?omega=(\S+) v_c=(\S+) v_c_max=(\S+)", re.M)


def csv_list(text, cast=str):
    return [cast(x.strip()) for x in text.split(",") if x.strip()]


def make_cases(a):
    cases, skipped = [], []
    for rho0, mode, delay, fast, lor, integ, rot, seed in itertools.product(
            a.rho0, a.modes, a.delay, a.fast, a.lorentz, a.integrator, a.rot, a.seed):
        if mode in NEEDS_DELAY and delay == 0:
            skipped.append("rho0={:g} mode={} delay=0".format(rho0, mode))
            continue
        cases.append(dict(rho0=rho0, mode=mode, delay=delay, fast=fast,
                          lorentz=lor, integrator=integ, rot=rot, seed=seed))
    return cases, skipped


def case_tag(a, c):
    tag = "rho0={:g}_{}_delay{}_fast{}_rot{}_{}_{}_R0{:g}_m{}_s{}_ev{}".format(
        c["rho0"], c["mode"], c["delay"], c["fast"], c["rot"], c["lorentz"],
        c["integrator"], a.R0, a.mass, c["seed"], a.maxeval)
    return re.sub(r"[^\w.=+-]", "-", tag)


def build_cmd(a, c):
    cmd = [a.exe,
           "--integrator={}".format(c["integrator"]),
           "--R0={}".format(a.R0),
           "--rho0={}".format(c["rho0"]),
           "--c={}".format(a.c),
           "--delay={}".format(c["delay"]),
           "--fast={}".format(c["fast"]),
           "--rot={}".format(c["rot"])]
    cmd.append("--seed={}".format(c["seed"]))
    cmd.append("--maxeval={}".format(a.maxeval))
    if a.omega is not None:
        cmd.append("--omega={}".format(a.omega))
    elif a.v_c is not None:
        cmd.append("--v-c={}".format(a.v_c))
    else:
        cmd.append("--mass={}".format(a.mass))
    cmd += LORENTZ[c["lorentz"]]
    cmd += MODES[c["mode"]]
    if c["mode"] in ("film", "film_inst") and a.film_K:
        cmd.append("--film-K={}".format(a.film_K))
    if c["mode"] == "pde" and a.pde_n:
        cmd.append("--pde-n={}".format(a.pde_n))
    return cmd


def to_text(x):
    if x is None:
        return ""
    if isinstance(x, bytes):
        return x.decode("utf-8", "replace")
    return x


def parse_output(text):
    res = dict((name, (float("nan"),) * 3) for name in RESULT_NAMES)
    for m in RESULT_RE.finditer(text):
        name, val, err, p = m.groups()
        if name in res:
            res[name] = (float(val), float(err), float(p))
    st = STATUS_RE.search(text)
    pm = PARAMS_RE.search(text)
    return (res,
            int(st.group(1)) if st else None,
            tuple(float(x) for x in pm.groups()) if pm else (None, None, None))


def run_case(a, c, logdir):
    cmd = build_cmd(a, c)
    env = os.environ.copy()
    # Cuba по умолчанию запускает форки на всех ядрах: при --jobs > 1 это
    # приводит к перегрузке. 0 = без внутреннего распараллеливания.
    env["CUBACORES"] = str(a.cuba_cores)
    t0 = time.time()
    note = ""
    try:
        p = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                           timeout=a.timeout, env=env)
        out, err, rc = to_text(p.stdout), to_text(p.stderr), p.returncode
    except subprocess.TimeoutExpired as e:
        out, err, rc = to_text(e.stdout), to_text(e.stderr), "timeout"
    except OSError as e:
        out, err, rc = "", str(e), "oserror"
    seconds = time.time() - t0

    first_err = ""
    for ln in err.splitlines():
        if ln.strip():
            first_err = ln.strip()
            break
    if rc != 0:
        note = first_err

    tag = case_tag(a, c)
    with io.open(os.path.join(logdir, tag + ".log"), "w", encoding="utf-8") as f:
        f.write(u"CMD: {}\nRC: {}\n--- stdout ---\n{}\n--- stderr ---\n{}\n".format(
            " ".join(cmd), rc, out, err))

    res, fail, (omega, v_c, v_c_max) = parse_output(out)
    row = dict(c, maxeval=a.maxeval, R0=a.R0, mass=a.mass, omega=omega, v_c=v_c, v_c_max=v_c_max,
               fail=fail, rc=rc, seconds=round(seconds, 1), note=note)
    for name in RESULT_NAMES:
        val, e, pr = res[name]
        row[name], row[name + "_err"], row[name + "_p"] = val, e, pr
    return row


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--exe", default="./calc_transverse_sphere_mass_film")
    ap.add_argument("--rho0", type=lambda s: csv_list(s, float), default=[0.92, 2.0, 100.0])
    ap.add_argument("--modes", type=csv_list,
                    default=["none", "fermi_O", "fermi_P", "general_P", "general_O"],
                    help="из: " + ",".join(sorted(MODES)))
    ap.add_argument("--delay", type=lambda s: csv_list(s, int), default=[0, 1])
    ap.add_argument("--fast", type=lambda s: csv_list(s, int), default=[1])
    ap.add_argument("--lorentz", type=csv_list, default=["off"],
                    help="из: " + ",".join(sorted(LORENTZ)))
    ap.add_argument("--integrator", type=csv_list, default=["vegas"])
    ap.add_argument("--R0", type=float, default=1.0)
    ap.add_argument("--c", type=float, default=1.0)
    ap.add_argument("--mass", default="1", help="1 или 4/3 (по умолчанию 1)")
    ap.add_argument("--omega", type=float, default=None, help="фиксировать omega для всех запусков")
    ap.add_argument("--v-c", dest="v_c", type=float, default=None, help="фиксировать v/c центра")
    ap.add_argument("--film-K", dest="film_K", type=int, default=0)
    ap.add_argument("--pde-n", dest="pde_n", type=int, default=0)
    ap.add_argument("--jobs", type=int, default=1)
    ap.add_argument("--timeout", type=float, default=3600.0, help="секунд на один запуск")
    ap.add_argument("--cuba-cores", dest="cuba_cores", type=int, default=0)
    ap.add_argument("--out", default="scan_results.csv")
    ap.add_argument("--logdir", default="scan_logs")
    ap.add_argument("--dry-run", action="store_true")
    ap.add_argument("--rot", type=lambda s: csv_list(s, int), default=[1],
                    help="повернуть поля в лаб. оси: 1, 0 или 0,1 для сравнения")
    ap.add_argument("--seed", type=lambda s: csv_list(s, int), default=[38304],
                    help="seed Cuba; список через запятую даёт независимые прогоны")
    ap.add_argument("--maxeval", type=int, default=500000)
    a = ap.parse_args()

    for m in a.modes:
        if m not in MODES:
            ap.error("неизвестный режим '{}'".format(m))
    for l in a.lorentz:
        if l not in LORENTZ:
            ap.error("неизвестный --lorentz '{}'".format(l))

    cases, skipped = make_cases(a)
    for s in skipped:
        print("SKIP (film/pde требуют delay=1): {}".format(s))
    print("Запусков: {}".format(len(cases)))
    if a.dry_run:
        for c in cases:
            print(" ".join(build_cmd(a, c)))
        return 0
    if not cases:
        return 0

    logdir = a.logdir
    if not os.path.isdir(logdir):
        os.makedirs(logdir)

    fields = ["rho0", "mode", "delay", "fast", "rot", "lorentz", "seed", "maxeval", "integrator", "R0", "mass",
              "omega", "v_c", "v_c_max", "fail", "rc", "seconds", "note"]
    for name in RESULT_NAMES:
        fields += [name, name + "_err", name + "_p"]

    rows = []
    with io.open(a.out, "w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fields)
        w.writeheader()
        with ThreadPoolExecutor(max_workers=max(1, a.jobs)) as ex:
            futs = dict((ex.submit(run_case, a, c, logdir), c) for c in cases)
            for i, fu in enumerate(as_completed(futs), 1):
                row = fu.result()
                rows.append(row)
                w.writerow(row)
                f.flush()
                print("[{}/{}] rho0={:g} {} delay={} fast={} rot={} {} seed={} rc={} fail={} {}s".format(
                    i, len(cases), row["rho0"], row["mode"], row["delay"], row["fast"],
                    row["rot"], row["lorentz"], row["seed"], row["rc"], row["fail"], row["seconds"]))

    rows.sort(key=lambda r: (r["rho0"], r["mode"], r["delay"], r["fast"], r["rot"], r["lorentz"], r["seed"]))
    print("\n rho0   mode       delay fast  rot lorentz  seed   ratio_rho (итог)            ratio_phi    fail  rc")
    for r in rows:
        line = "{:6g}  {:<10} {:>3} {:>4} {:>3} {:<8} {:>6} {:>12.6f} +- {:<10.2g}{:>12.6f}  {:>4}  {}".format(
            r["rho0"], r["mode"], r["delay"], r["fast"], r["rot"], r["lorentz"], r["seed"],
            r["observed_ratio_rho"], r["observed_ratio_rho_err"],
            r["observed_ratio_phi"], str(r["fail"]), r["rc"])
        if r["note"]:
            line += "   !! " + r["note"]
        print(line)
    print("\nCSV: {}, логи: {}/".format(a.out, logdir))
    return 0 if all(r["rc"] == 0 for r in rows) else 1


if __name__ == "__main__":
    sys.exit(main())
