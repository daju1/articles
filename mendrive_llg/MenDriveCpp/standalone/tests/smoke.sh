#!/bin/sh
# Короткий smoke-тест автономной программы: Preisach2D на малой и на «боевой»
# сетке, замер пикового потребления памяти и проверка --resume.
# Запуск из папки standalone:  sh tests/smoke.sh      (или make smoke)
set -e
cd "$(dirname "$0")/.."
[ -x ./mendrive ] || { echo "сначала соберите: make"; exit 1; }

if [ "$(uname -s)" = "Darwin" ]; then TIMECMD="/usr/bin/time -l"; else TIMECMD="/usr/bin/time -v"; fi
COMMON="N=160 omega0=3.0 bias_orientation=x H0_bias=4 n_sub=2 n_fp=2 n_newton=4 sigma_m_leak=3 sigma_e_left=3 ts_every=0 fft_bins=50"

echo "=== 1) Preisach2D rows, n_grid=41, 3 периода"
$TIMECMD ./mendrive $COMMON ferrite_model=Preisach2D pr2d_algo=rows n_grid=41 \
    n_periods_total=3 record_from_period=1 out_dir=tests/out_smoke_rows41 2>tests/time_rows41.txt
grep -E "maximum resident|Maximum resident" tests/time_rows41.txt || true

echo "=== 2) Preisach2D bits, n_grid=41 (должно совпасть с rows до ~1e-12)"
./mendrive $COMMON ferrite_model=Preisach2D pr2d_algo=bits n_grid=41 \
    n_periods_total=3 record_from_period=1 out_dir=tests/out_smoke_bits41 >/dev/null
grep -E "^(dTxx_avg|P_src_avg|Pdiss_avg|hysteresis_area)=" tests/out_smoke_rows41/summary.txt
grep -E "^(dTxx_avg|P_src_avg|Pdiss_avg|hysteresis_area)=" tests/out_smoke_bits41/summary.txt

echo "=== 3) Preisach2D rows, n_grid=2001 (как n_hyst_pr=2001 в ноутбуке), 1 период: память и скорость"
$TIMECMD ./mendrive $COMMON ferrite_model=Preisach2D pr2d_algo=rows n_grid=2001 \
    n_periods_total=1.2 record_from_period=0.2 progress_every_periods=1 \
    out_dir=tests/out_smoke_rows2001 2>tests/time_rows2001.txt
grep -E "maximum resident|Maximum resident" tests/time_rows2001.txt || true
grep -E "^(wall_time_s|peak_rss_mb)=" tests/out_smoke_rows2001/summary.txt

echo "=== 4) --resume: прерываем прогон и продолжаем; итог должен совпасть с непрерывным"
./mendrive $COMMON ferrite_model=Preisach n_hyst_pr=201 n_periods_total=40 record_from_period=10 \
    checkpoint_every_periods=1 out_dir=tests/out_smoke_full >/dev/null
./mendrive $COMMON ferrite_model=Preisach n_hyst_pr=201 n_periods_total=40 record_from_period=10 \
    checkpoint_every_periods=1 out_dir=tests/out_smoke_resume >/dev/null &
PID=$!
sleep 3
kill -INT $PID 2>/dev/null || true
wait $PID || true
./mendrive $COMMON ferrite_model=Preisach n_hyst_pr=201 n_periods_total=40 record_from_period=10 \
    checkpoint_every_periods=1 out_dir=tests/out_smoke_resume --resume >/dev/null
for f in tests/out_smoke_full tests/out_smoke_resume; do
    grep -vE "^(wall_time_s|peak_rss_mb)=" $f/summary.txt > $f/summary_cmp.txt
done
if cmp -s tests/out_smoke_full/summary_cmp.txt tests/out_smoke_resume/summary_cmp.txt \
   && cmp -s tests/out_smoke_full/periods.csv tests/out_smoke_resume/periods.csv; then
    echo "resume: OK (итоги и periods.csv совпали бит-в-бит)"
else
    echo "resume: РАСХОЖДЕНИЕ"; diff tests/out_smoke_full/summary_cmp.txt tests/out_smoke_resume/summary_cmp.txt || true
fi
