# MenDrive — автономная C++-версия (без Python / ctypes / Jupyter)

Папка `standalone/` рядом с исходными файлами `../mendrive_core.cpp`, `../mendrive_ctypes.py`,
`../MenDriveCpp_run.ipynb` (их я **не менял**). Здесь:

| файл | назначение |
|---|---|
| `mendrive_core.cpp` | копия ядра с минимальными исправлениями (`// FIX:`), оптимизациями без изменения арифметики (`// PERF:`) и новой диагностикой (`// DIAG:`). extern "C" API совместим со старой `mendrive_ctypes.py` |
| `main.cpp` | порт `main_cpp` / `scan_resonance_cpp` / `force_per_power` / `check_energy_balance` / `shoelace_area` + потоковая запись, чекпойнты, `--resume`, прогресс |
| `Makefile` | сборка (`make`, `make omp`, `make lib`, `make reflib`, `make test`, `make smoke`) |
| `examples/*.cfg` | готовые конфиги: повтор ячейки ноутбука, Preisach2D, быстрая проверка |
| `tests/compare_ctypes.py` | сравнение с исходным ctypes-путём (старое ядро против нового ядра против `./mendrive`) |
| `tests/smoke.sh` | smoke-тест Preisach2D, замер пиковой памяти, проверка `--resume` |
| `plot.py` | необязательные графики по готовым CSV (только чтение) |
| `linear1d.py` | **MODE/FIELDS:** точное линейное 1D-решение (эталон), чтение/запись CSV профилей, экспорт из ноутбуков `MenDrive_qnm*` (`export_sage_profiles`), CLI `python3 linear1d.py out_dir` |
| `fdtd1d_ref.py` | numpy-копия линейной части нового `step_once` и накопления гармоник (проверка C++ без компилятора) |
| `compare_linear_vs_cpp.ipynb` | ноутбук-помощник: проверки «формулы → OK/FAIL», графики в стиле `MenDrive_qnm*`, сверка с C++ |
| `tests/test_fields_linear.py` | проверка новых мод TE/TM и экспорта полей (`make test-fields`) |
| `backup_before_modes/` | версии файлов до этих изменений |

> **Важно о проверке.** На этом Mac (M1 Pro, 16 ГБ) сейчас нет компилятора: `/usr/bin/c++`
> отвечает `xcode-select: note: No developer tools were found`. В облачной песочнице, где я
> работал, запуск процессов запрещён. Поэтому **собрать и запустить C++-код я не смог**. Код
> вычитан вручную, а ключевые алгоритмические изменения проверены отдельно на Python/numpy
> (см. раздел «Что проверено»). Первым делом выполните `make test`: он сравнит новую программу
> с исходным ядром до машинной точности.

---

## 1. Почему падало ядро Python в режиме `ferrite_model='Preisach2D'`

### 1.1. `n_hyst_pr` в Preisach2D — это `n_grid` по **каждой** оси, т.е. $n_{grid}^2$ пар

В `MenDriveSim` для Preisach2D вызывается `PreisachEngine2D(..., n_hyst_pr)`, и этот параметр
становится `n_grid`. Сетка строится как $n_{grid}\times n_{grid}$. При типичном для 1D-Прейзаха
`n_hyst_pr=2001` получается $2001^2 = 4\,004\,001$ пар **на каждую ячейку феррита**.

Память исходного кода (C++-часть), $n_{fer}$ — число ячеек феррита ($n_{fer}=40$ при $N=160$):

$$
M_{C++} \approx \underbrace{n_{fer}\,n_{grid}^2}_{\text{int8-состояния}}
 + \underbrace{3\,n_{grid}^2\,s_{ld}}_{\alpha_k,\beta_k,w_k}
 + \underbrace{n_{grid}^2}_{\text{копия trial на итерацию}},
\qquad s_{ld}=\texttt{sizeof(long double)}=8\ (\text{arm64}),\ 16\ (\text{x86-64}).
$$

| $n_{grid}$ | пар | int8-состояние, 40 яч. | $\alpha,\beta,w$ (arm64) | bits (новое) | rows (новое) |
|---:|---:|---:|---:|---:|---:|
| 41 | 1 681 | 64 КБ | 39 КБ | 8 КБ | 6 КБ |
| 201 | 40 401 | 1.5 МБ | 0.9 МБ | 0.2 МБ | 31 КБ |
| 2001 | 4 004 001 | **153 МБ** | **92 МБ** (184 МБ на x86) | 19 МБ | 0.3 МБ |
| 6001 | 3.6·10⁷ | **1.37 ГБ** | **0.82 ГБ** | 172 МБ | 0.9 МБ |
| 20001 | 4·10⁸ | **15 ГБ** | **9 ГБ** | 1.9 ГБ | 3 МБ |
| 46341 | 2.15·10⁹ | переполнение `int` (UB) | — | 10 ГБ | 7 МБ |

### 1.2. Время: исходный Preisach2D физически не досчитывается

Стоимость шага исходного Preisach2D:

$$
C_{step} \propto n_{fer}\,(n_{newton}+1)\,n_{grid}^2\quad(\text{копия состояния + проход переключения + проход суммы}).
$$

Калибровка по ноутбуку: 1D-Прейзах с 2001 гистероном считал $33.5\cdot10^6$ шагов за 27 385 с,
т.е. $\approx 0.82$ мс/шаг ($\approx 2$ нс на одно «реле × итерацию»). Для Preisach2D с
$n_{grid}=2001$ работы в $n_{grid}^2/n_{hyst}=2001$ раз больше: $\approx 1.6$ с/шаг, весь
прогон на 8000 периодов $\approx 5\cdot10^7$ с $\approx$ **1.7 года**. Даже скан резонанса
не заканчивается.

### 1.3. Память Python-стороны растёт пропорционально числу записанных шагов

`MenDriveCpp.run()` заранее выделяет 10 массивов float64 длиной
$N_{rec} = \lfloor n_{tot}T/\Delta t\rfloor - \lfloor n_{from}T/\Delta t\rfloor$, где
$\Delta t = 0.1\,\Delta x$, $\Delta x = (2a+h_l+h_r)/N$. Затем `main_cpp` добавляет `Mn` (3 столбца),
окно Ханна, произведение, `rfft` (комплексный массив $N_{rec}/2$) и `rfftfreq`:

$$
M_{Py} \approx 8\,\text{Б}\cdot N_{rec}\,(10 + 3 + 2) + 16\,\text{Б}\cdot\tfrac{N_{rec}}{2} + 2\cdot 8\,\text{Б}\cdot\tfrac{N_{rec}}{2}
\approx 136\ \text{Б}\cdot N_{rec}.
$$

Для $N=160$, $\omega_0=3$, 8000/1000 периодов: $N_{rec}=29\,321\,531$, то есть **≈3.7 ГБ на
один вызов `main_cpp`** (2.2 ГБ буферов записи, 0.66 ГБ `Mn`, 0.87 ГБ FFT), плюс отрисовка
14.6 млн точек спектра в matplotlib. Пока создаётся новый `result`, старый ещё жив, поэтому
пик в цепочке ячеек ≈ **7–8 ГБ**. Вместе с состоянием Preisach2D (см. 1.1) и самим Jupyter это
выходит за пределы 16 ГБ.

### 1.4. Прочие дефекты, которые могут уронить процесс

* исключения (`std::bad_alloc` при нехватке памяти и др.) вылетали из конструктора через
  границу `extern "C"`, что ведёт к `std::terminate`, т.е. к падению процесса Python и ядра Jupyter
  (**FIX**: `try/catch`, функция возвращает `nullptr`);
* `probe_idx >= n_fer` давал чтение за границами `M_field/pr_M/...` (**FIX**: проверка);
* `n_grid*n_grid` в `int` переполнялся при $n_{grid}>46340$ (**FIX**: `int64_t`);
* `sigma<=0` давал деление на ноль в весах, `n_grid<1` — деление на `wsum=0` (**FIX**: проверки);
* на каждой итерации Ньютона создавалась копия вектора состояния (лишние аллокации; **FIX**);
* рекурсии и утечек памяти нет (всё в `std::vector`/`unique_ptr`).

**Вывод.** Самая вероятная причина падения — исчерпание памяти. Её съедают две вещи:
буферы записи на Python-стороне, растущие как $O(N_{rec})$ (≈3.7 ГБ на вызов, до 7–8 ГБ в
цепочке ячеек), и состояние Preisach2D, растущее как $O(n_{fer}n_{grid}^2)$, потому что
`n_hyst_pr` здесь — это `n_grid`. Вместе это превышает 16 ГБ, процесс убивает ОС или
`bad_alloc` → `terminate`. Но даже без падения исходный алгоритм с $n_{grid}=2001$ работал бы
годами.

---

## 2. Что изменено в ядре

1. **Preisach2D: два точных представления состояния** (`pr2d_algo`):
   * `bits` — 1 бит на пару (бит-упаковка, память /8). Порядок суммирования $M$ тот же, что в
     оригинале, поэтому числа совпадают с исходным кодом. Сложность $O(n_{grid}^2)$.
   * `rows` (по умолчанию в `main.cpp`) — для фиксированной строки $hc_i$ пороги
     $\alpha_{ij}=hi_j+hc_i$ и $\beta_{ij}=hi_j-hc_i$ монотонны по $j$. Рост поля включает
     префикс строки, спад выключает суффикс. Отсюда состояние строки всегда имеет вид
     $(+1,\dots,+1,-1,\dots,-1)$ и задаётся одним целым $K_i$:
     $$K_i \leftarrow \max\bigl(K_i,\ \#\{j: H\ge hi_j+hc_i\}\bigr)\ \text{(рост)},\qquad
       K_i \leftarrow \min\bigl(K_i,\ \min\{j: H\le hi_j-hc_i\}\bigr)\ \text{(спад)},$$
     $$M=\frac{M_s}{W}\sum_i w^{hc}_i\bigl(2S_{K_i}-S_{n}\bigr),\qquad S_K=\sum_{j<K}w^{hi}_j .$$
     Индекс находится за $O(1)$: оценка по равномерной сетке плюс точная коррекция теми же
     сравнениями, что в оригинале. Память $O(n_{fer}n_{grid})$, время $O(n_{grid})$ вместо
     $O(n_{grid}^2)$. Состояния реле совпадают с полным перебором **точно**, $M$ отличается
     только округлением ($\le 2\cdot10^{-15}$, проверено, см. §6).
2. Массивы $\alpha_k,\beta_k,w_k$ не хранятся: они вычисляются из сеток тем же выражением.
3. Пробный шаг Ньютона идёт без копирования состояния (1D и 2D Прейзах).
4. Рабочие буферы шага переиспользуются, копирование полей заменено на `swap`. Значения те же.
5. Многопоточность по ячейкам феррита (JA, Preisach, Preisach2D): `threads=K`. Пул на
   `std::thread` работает без внешних библиотек. Каждая ячейка пишет только своё, суммы
   собираются последовательно в фиксированном порядке, поэтому **результат не зависит от числа
   потоков**. Сборка `make omp` вместо пула использует OpenMP.
6. **DIAG:** невязка неявного шага $P_{num}$ и составляющие $P_{diss}$ (см. §7).
7. Чекпойнт полного динамического состояния (`save_state/load_state`).

`mendrive_run()` и сигнатура `mendrive_create()` не изменились, поэтому старую обёртку можно
подключить к новому ядру: `make lib` → `libmendrive_new.*`.

---

## 3. Сборка

### macOS (Apple Silicon / Intel)

```sh
xcode-select --install        # один раз: Command Line Tools (clang, make)
cd ~/articles/mendrive_llg/MenDriveCpp/standalone
make                          # -> ./mendrive   (флаги: -O3 -mcpu=native|-march=native -std=c++17)
```

На macOS `g++` — это Apple clang. Он не поддерживает `-fopenmp` без libomp. Многопоточность
`threads=K` работает и без OpenMP (через `std::thread`). OpenMP-вариант собирается так:

```sh
brew install libomp           # нужен Homebrew
make omp                      # -> ./mendrive_omp (Makefile сам найдёт brew --prefix libomp)
```

Если компилятор не понимает `-mcpu=native`/`-march=native`, используйте `make ARCH=`.

### Linux

```sh
make CXX=g++                  # g++ -O3 -march=native -std=c++17 main.cpp -o mendrive -pthread
make omp CXX=g++              # g++ ... -fopenmp
```

---

## 4. Запуск

```sh
./mendrive --help                                   # все параметры и значения по умолчанию
./mendrive -c examples/smoke_ja.cfg                 # секунды: скан + короткий прогон JA
./mendrive -c examples/notebook_preisach_magnetic.cfg            # повтор ячейки 4 ноутбука
./mendrive -c examples/notebook_preisach_magnetic.cfg amp=5 out_dir=out/pr_mag_amp5
./mendrive -c examples/preisach2d_magnetic.cfg                   # Preisach2D, n_grid=201, rows
./mendrive -c examples/preisach2d_magnetic.cfg n_grid=2001 out_dir=out/p2d_2001

# долгий прогон в фоне, Mac не засыпает, лог в файл:
nohup caffeinate -i ./mendrive -c examples/preisach2d_magnetic.cfg > out/p2d.log 2>&1 &
tail -f out/p2d.log

# остановка: Ctrl+C или kill -INT <pid>  -> пишется чекпойнт; продолжение:
./mendrive -c examples/preisach2d_magnetic.cfg --resume
```

Параметры задаются в файле (`ключ = значение`, `#` — комментарий) и/или в командной строке
`ключ=значение`. Они применяются по порядку, поэтому параметр в командной строке после `-c`
переопределяет значение из файла. Имена совпадают с аргументами `main_cpp`/`MenDriveCpp`. Новые
ключи: `n_grid`, `pr2d_algo`, `omega0`, `scan_*`, `fine_*`, `ts_*`, `fft_*`, `threads`,
`checkpoint_*`, `progress_every_periods`, `force_resume`.

Отличия от ноутбука, которые нужно знать:

* `omega0>0` — частота задана, скан пропускается. Иначе выполняется скан как в `main_cpp`:
  `arange(scan_lo, scan_hi, scan_step)`, затем уточнение $\pm$`fine_halfwidth` с шагом
  `fine_step`; в отклике проба 0, `amp=0.3`, 3 периода. Сетка по умолчанию из кода
  `arange(0.005, 2.005, 0.005)` стоит $\sum_\omega 3\cdot 2\pi/(\omega\Delta t)\approx
  4.95\cdot10^7$ шагов, это **больше длинного прогона** (низкие частоты). В сохранённых выводах
  ноутбука пик «3.00» лежит вне этой сетки, значит, при тех запусках сетка была другой (судя по
  времени скана — целые частоты 1..20). Скан параллелится по частотам (`threads`) и
  дописывается в `scan_*.csv` после каждой точки, поэтому прерванный скан продолжается по `--resume`.
* `sigma_m_leak_pr`: `main_cpp` передаёт в модель только `sigma_m_leak`, а для Прейзаха всегда
  остаётся значение по умолчанию 3.0. Здесь это отдельный ключ (по умолчанию тоже 3.0).
* Для Hybrid не повторяется вспомогательный LLG-прогон ноутбука (он нужен только для рисунка
  траектории).
* Полный FFT заменён ДПФ с тем же окном Ханна на гармониках $k\omega_0$ и на равномерной сетке
  `fft_bins` точек в $[0,\ \texttt{fft\_wmax\_factor}\cdot\omega_0]$. По определению это
  `np.fft.rfft(dTxx*np.hanning(N_rec))` в этих частотах.

---

## 5. Выходные файлы (`out_dir/`)

| файл | содержимое |
|---|---|
| `config_used.txt` | все параметры прогона |
| `scan_wide.csv`, `scan_fine.csv` | `idx,omega,resp` (порядок строк — по завершении) |
| `scan_result.txt` | `omega_coarse`, `omega_res` |
| `periods.csv` | по одной строке на период записи: `block,t_end,dTxx,P,P_src,Pdiss,P_num,Hn_min,Hn_max,W_field,cum_dTxx_over_P,cum_dTxx_over_Psrc` (средние за период; `cum_*` = `convergence_ratio_cum` ноутбука) |
| `timeseries.csv` / `.bin` | прореженный хвост: каждый `ts_every`-й шаг последних `ts_last_periods` периодов; столбцы `t,dTxx,Hn,P,P_src,Pdiss,Pe,Pm_wall,Pmag,P_num,Mx,My,Mz,HzJA,MzJA,W_field` (для `.bin` — float64 подряд, имена в `ts_columns.txt`) |
| `loop.csv` | последний период: `H,M,B` ($B=H+4\pi M$), по нему `hysteresis_area` |
| `harmonics.csv`, `spectrum.csv` | $\lvert\mathrm{DFT}\rvert$ (и фаза) $dT_{xx}$ |
| `summary.txt` | итоги `ключ=значение`: `dTxx_avg, P_avg, P_src_avg, Pdiss_avg` и его части, `P_num_residual`, `force_per_power_code`, `force_per_kW` (как в ноутбуке) и `*_src` (с правильной мощностью), `rel_dif`, `rel_dif_src`, `rel_dif_src_plus_num`, `hysteresis_area`, `wall_time_s`, `peak_rss_mb` |
| `checkpoint.bin` | состояние для `--resume`. Пишется атомарно, раз в `checkpoint_every_periods` периодов и не реже `checkpoint_every_sec`, а также по Ctrl+C/SIGTERM |

`P` — та же величина, что `res['P']` в ноутбуке (с ошибкой масштаба $1/\Delta t$, см. §7).
`P_src` — физическая мощность источника.

Графики (необязательно): `python3 plot.py out/pr_mag_amp1` → `out/pr_mag_amp1/plots/*.png`.

---

## 6. Память, время и что проверено

**Память программы** не зависит от длины прогона:

$$
M \approx M_{state} + O(N) + 2\cdot 8\,\text{Б}\cdot\mathrm{spp} + 5\cdot 8\,\text{Б}\cdot n_{bins},
\qquad M_{state}^{rows}=4\,n_{fer}n_{grid}\ \text{Б},\quad M_{state}^{bits}=n_{fer}\lceil n_{grid}^2/64\rceil\cdot 8\ \text{Б},
$$

где $\mathrm{spp}=\mathrm{round}(T/\Delta t)$ (4189 при $\omega_0=3$). Для $N=160$, Preisach2D,
$n_{grid}=2001$, rows это **меньше 2 МБ** данных (весь процесс — единицы МБ) против ≈4 ГБ + 0.25 ГБ
в ноутбуке. Программа печатает оценку при старте, а `peak_rss_mb` записывает в `summary.txt`.

**Время** (оценка по калибровке ноутбука ≈2 нс на «реле × итерацию»; точные цифры даст `make smoke`):

| режим ($N=160$, 8000 периодов, $3.35\cdot10^7$ шагов) | на шаг | прогон, 1 поток |
|---|---:|---:|
| Preisach 1D, 2001 гистерон (ноутбук) | 0.82 мс | 7.6 ч (факт) |
| Preisach2D исходный / `bits`, $n_{grid}=2001$ | ~1.6 с | ~1.7 года |
| Preisach2D `bits`, $n_{grid}=41$ | ~0.7 мс | ~6 ч |
| Preisach2D `rows`, $n_{grid}=2001$ | ~2–4 мс | ~20–35 ч; при `threads=6` ~4–7 ч |
| Preisach2D `rows`, $n_{grid}=201$ | ~0.2–0.4 мс | ~2–4 ч |

Время одного шага растёт линейно по $n_{grid}$ и почти линейно делится на `threads`
(40 независимых ячеек). Сначала сделайте короткий прогон (`n_periods_total=20`): по строке
прогресса `шаг/с` и `осталось ~` видно реальное время.

**Что проверено (без компилятора C++, на Python/numpy):**

* Алгоритм `rows` против полного перебора реле (как в оригинальном `step`) на 3500 случайных
  траекториях поля и при $n_{grid}=1,2,7,41$: **0 несовпадений состояний**,
  $\max\lvert M_{bits}-M_{rows}\rvert \le 2\cdot10^{-15}$.
* Независимая numpy-реализация той же схемы (N=160, $\omega_0=3$, amp=1, σ=3, n_newton=4)
  замыкает баланс энергии через $P_{num}$ (см. §7).

**Что проверить у себя после сборки:** `make test`. Тест прогоняет JA, Preisach (обе
накачки), Preisach2D bits/rows, LLG и Hybrid на N=40. Для всех режимов, кроме rows, ожидается
совпадение `dTxx_avg, P_avg, Pdiss_avg, force_per_power_code, rel_dif, hysteresis_area` с исходным
ctypes-путём до ~1e-12; для rows — до ~1e-9. Затем `make smoke`: Preisach2D с $n_{grid}=41$ и
$2001$, пиковая RSS через `/usr/bin/time -l`, и проверка, что прерванный и продолженный по
`--resume` прогон совпадает с непрерывным бит-в-бит. Если что-то не собирается или не
совпадает, пришлите вывод.

Замечание о воспроизводимости 1D Preisach: `std::normal_distribution` реализован по-разному в
libstdc++ (Linux g++) и libc++ (macOS clang). Поэтому гистероны 1D-Прейзаха на Mac и на Linux
**разные** (статистически эквивалентные), и бит-в-бит сравнивать можно только на одной платформе.
Preisach2D (детерминированная сетка) от этого не зависит.

---

## 7. Замечания по балансу энергии (важно для вопроса о диссипативной тяге)

### 7.1. Что должно выполняться

В гауссовых единицах (в коде $c=1$) теорема Пойнтинга для этой схемы:

$$
\frac{\partial}{\partial t}\frac{E^2+H^2}{8\pi}+\nabla\!\cdot\!\mathbf S=
-\mathbf j_e\!\cdot\!\mathbf E-\mathbf j_m\!\cdot\!\mathbf H-\sigma_eE^2-\sigma_mH^2-\mathbf H\!\cdot\!\frac{\partial\mathbf M}{\partial t}.
$$

Концы расчётной области — идеальные стенки ($E=0$), поэтому в установившемся периодическом
режиме среднее по периоду равно

$$
\langle P_{src}\rangle=\langle P_{diss}\rangle,\qquad
P_{src}=-\!\int\!(j_eE_y+j_mH_z)\,dx,\qquad
P_{diss}=\!\int\!\sigma_eE^2dx+\!\int_{fer}\!\sigma_mH^2dx+\!\int_{fer}\!H\,\dot M\,dx .
$$

### 7.2. Ошибка нормировки мощности в `mendrive_run`: $P^{code}=P_{src}/\Delta t$

`step_once()` уже возвращает **мощность** $P_{src}$: в сумме $-\sum j\,\bar E\,\Delta x$ нет
множителя $\Delta t$. `mendrive_run()` записывает `out_P = P_step / sim->dt`, то есть
$P^{code}=P_{src}/\Delta t$. При $N=160$: $\Delta t=0.1\cdot 0.8/160=5\cdot10^{-4}$, и мощность
завышена **в 2000 раз**.

Подтверждение — выводы самого ноутбука (electric_left):

| amp | 1 | 2 | 3 | 4 | 5 | 6 | 7 | 8 | 9 | 10 |
|---|---|---|---|---|---|---|---|---|---|---|
| $P^{code}/P_{diss}$ | 1981 | 1978 | 1975 | 1964 | 1950 | 1924 | 1879 | 1826 | 1767 | 1707 |
| $P_{src}/P_{diss}=P^{code}\Delta t/P_{diss}$ | 0.991 | 0.989 | 0.987 | 0.982 | 0.975 | 0.962 | 0.939 | 0.913 | 0.884 | 0.853 |

При правильной нормировке баланс для electric_left сходится до 1% при малой амплитуде. Остаток
растёт с amp по мере включения переключений реле: это численная невязка неявного шага (§7.4).
Константа $P_{diss}$ записывается правильно. Ошибка именно в $P$, и значит, `force_per_kW`
в ноутбуке занижена в 2000 раз при любом amp.

В `main.cpp` выводятся обе величины: `P`/`force_per_kW`/`rel_dif` (как в ноутбуке) и
`P_src`/`force_per_kW_src`/`rel_dif_src` (правильно). Исходный `mendrive_run` оставлен как был
(с пометкой `NOTE`), чтобы старый ноутбук давал прежние числа.

### 7.3. magnetic_right: $P_{diss}\approx190.8$ и $\langle dT_{xx}\rangle\approx1.97$ не зависят от amp

После пересчёта $P_{src}=P^{code}\Delta t$ (по выводам ноутбука, Preisach 1D, magnetic_right):

| amp | 1 | 2 | 3 | 5 | 10 |
|---|---|---|---|---|---|
| $P_{src}$ | 1.72·10⁻³ | 6.23·10⁻³ | 1.51·10⁻² | 4.07·10⁻² | 0.163 |
| $P_{diss}$ | 190.83 | 190.80 | 190.79 | 190.81 | 190.77 |
| $\langle dT_{xx}\rangle$ = force_per_power_code·P_avg | 1.972 | 1.972 | 1.972 | 1.974 | 1.972 |
| hysteresis_area | 1.70·10⁴ | 1.58·10⁴ | 1.65·10⁴ | 1.55·10⁴ | 1.66·10⁴ |

Отсюда следует:

1. **Рассеивается в $10^5$ раз больше энергии, чем подаёт источник** ($190.8$ против
   $1.7\cdot10^{-3}$ при amp=1). В установившемся режиме это запрещено теоремой Пойнтинга.
   Значит, в расчёте есть источник энергии, которого нет в физике.
2. $P_{diss}$, $\langle dT_{xx}\rangle$ и площадь петли **не зависят от амплитуды накачки**,
   хотя $P_{src}\propto \mathrm{amp}^2$ ($P^{code}\approx 3.3\,\mathrm{amp}^2$). Это
   самоподдерживающийся режим, не связанный с накачкой. Константа 190.8, которую вы заподозрили,
   действительно подозрительна, но формула $P_{diss}$ её не создаёт. Это реальное (численное)
   состояние полей.
3. Отсюда $\text{force\_per\_kW}=\langle dT_{xx}\rangle/P \propto 1/\mathrm{amp}^2$, а не
   $1/\mathrm{amp}$: $57238/605 = 94.6\approx 10^2$. «Тяга» не масштабируется с подводимой
   мощностью. Это постоянное смещение $\langle dT_{xx}\rangle$, порождённое тем же паразитным
   режимом, а не диссипативная тяга.
4. Для вакуумного промежутка $[-a,a]$ в 1D точно выполняется
   $\partial_t g_x + \partial_x T_{xx}=0$, откуда
   $\frac{d}{dt}\int_{-a}^{a}g_x\,dx = T_{xx}(-a)-T_{xx}(a)$. В периодическом установившемся
   режиме среднее левой части равно нулю, поэтому $\langle dT_{xx}\rangle\neq0$ означает одно из
   двух: либо режим не стационарен, либо это ошибка дискретизации (`Txx_at` смешивает $E$ и
   усреднённое по соседним ячейкам $H$; такая ошибка должна убывать при измельчении $N$).
   Значение $1.97$ при $P_{src}\sim10^{-3}$ не может быть линейным откликом на накачку.

### 7.4. Независимая проверка и вероятный механизм

Я переписал ту же схему (Yee + неявный шаг Прейзаха с демпфированным квазиньютоном, 4 итерации,
множитель 0.5) на numpy: N=160, $\omega_0=3$, amp=1, $\sigma_e=\sigma_m=3$; 1D-Прейзах с 201 и с
2001 гистероном и 2D-сетка 41×41. Выход на режим занимает около 6 периодов, дальше всё
стабильно до 5 знаков:

| | $P_{src}$ | $P_{diss}$ | $P_{num}$ | $P_{src}+P_{num}$ | $\langle dT_{xx}\rangle$ | $P_{src}/\Delta t$ |
|---|---|---|---|---|---|---|
| numpy-порт (1D, 2001 гист.) | 1.937·10⁻³ | 2.011·10⁻³ | 7.38·10⁻⁵ | **2.011·10⁻³** | −6·10⁻⁸ | 3.87 (ноутбук: 3.45) |

* Баланс замыкается **точно**, если учесть $P_{num}$ — мощность фиктивного источника из
  невязки неявного шага:
  $$F_{res}=H^{n+1}-H^n+4\pi(M^{n+1}-M^n)+4\pi\sigma_m\Delta t\,H^{n+1}+\Delta t(\nabla\times E)_z+4\pi\Delta t\,j_m,\qquad
    P_{num}=\int_{fer}\frac{\bar H\,F_{res}}{4\pi\Delta t}\,dx .$$
  Демпфированная итерация (множитель 0.5, 4 итерации) оставляет невязку $\sim 2^{-4}$ от
  начальной. Дискретный закон Фарадея нарушается на $F_{res}$, и схема **сама подкачивает
  энергию**: ≈4% при малых полях, больше при переключениях реле (у реле ступенька
  $4\pi\cdot2w_k$ в $B$, и точного решения неявного уравнения может не быть вовсе).
* В порте **нет** режима с $P_{diss}\approx190$ и $\langle dT_{xx}\rangle\approx2$:
  $\langle dT_{xx}\rangle$ на 7 порядков меньше, чем в ноутбуке. Возмущение $H\in[-3,3]$ в
  феррите без накачки (amp=0) затухает за 1–2 периода.

Значит, паразитный режим в ноутбуке возникает позже (у вас среднее считается по периодам
4500–8000), либо выводы получены другой сборкой `.so`: в папке `.so` нет, ноутбук запускали на
другой машине. Наиболее вероятный источник энергии — невязка $P_{num}$. Подтвердить это может
только прогон исходного кода. Новая программа пишет в `periods.csv` для каждого периода
$P_{src}$, $P_{diss}$, $P_{num}$, $\min/\max H_n$ и энергию поля, так что момент появления и
источник режима будут видны.

### 7.5. Что сделать, прежде чем делать выводы о тяге

1. `./mendrive -c examples/notebook_preisach_magnetic.cfg` (≈7–8 ч, можно `threads=6`).
   Посмотрите в `periods.csv`, когда $P_{diss}$ отрывается от $P_{src}$ и чему равно $P_{num}$.
2. Тот же прогон с `amp=0.1`: если $\langle dT_{xx}\rangle$ и $P_{diss}$ останутся ≈1.97 и ≈190,
   режим не связан с накачкой, и это артефакт.
3. `n_newton=20` (и/или меньше `dt_frac`): если $P_{num}\to0$ и режим исчезает, он численный.
4. Сходимость по сетке $N=80,160,320$: физическая величина $\langle dT_{xx}\rangle/P_{src}$
   должна сходиться, а артефакт дискретизации — убывать как $O(\Delta x^2)$.
5. Сравнивать тягу только при выполненном балансе: `rel_dif_src_plus_num` ≪ 1 и
   $\lvert P_{num}\rvert\ll P_{src}$.

Пока эти проверки не пройдены, `force_per_kW` из magnetic_right (от 57 238 до 605 «Н/кВт» в
код-единицах) **нельзя интерпретировать как диссипативную тягу**. Эта величина делит не
связанное с накачкой смещение $\langle dT_{xx}\rangle$ на мощность, завышенную в $1/\Delta t$ раз.


---

## 8. Моды и поляризации (`// MODE:`), соответствие ноутбукам `MenDrive_qnm*.ipynb`

### 8.1. Две независимые пары компонент

Модель одномерна: поля зависят только от $x$ ($k_y=k_z=0$). Уравнения Максвелла в Гауссовой системе
($c=1$, время в коде $\tau=ct$), записанные так же, как в ячейке 74 ноутбуков
($\mathrm{rot}\,\vec H=\partial_\tau\vec D+4\pi\vec j_e$, $\mathrm{rot}\,\vec E=-\partial_\tau\vec B-4\pi\vec j_m$),
распадаются на две пары:

$$\text{TE}=(E_y,H_z):\qquad \partial_\tau D_y=-\partial_x H_z-4\pi j_{e,y},\qquad \partial_\tau B_z=-\partial_x E_y-4\pi j_{m,z},$$

$$\text{TM}=(E_z,H_y):\qquad \partial_\tau D_z=+\partial_x H_y-4\pi j_{e,z},\qquad \partial_\tau B_y=+\partial_x E_z-4\pi j_{m,y}.$$

Названия TE/TM взяты из ноутбуков (`is_TE`: $H_z\ne0,\ E_z=0$; `is_TM`: $E_z\ne0,\ H_z=0$, ячейки 45, 78).
Тангенциальные к стенкам компоненты: $E_y,E_z,H_y,H_z$. Отсюда правило:

* $j_{e,y}$ и $j_{m,z}$ возбуждают **TE** $(E_y,H_z)$;
* $j_{e,z}$ и $j_{m,y}$ возбуждают **TM** $(E_z,H_y)$.

### 8.2. Таблица соответствия `driving_mode` ↔ `excitation_mode`

Сведена по всем трём ноутбукам `MenDrive_qnm.ipynb`, `MenDrive_qnm-mode2_je_z_TM.ipynb`,
`MenDrive_qnm-mode3_jm_z_TE.ipynb`. Набор `driving_mode_0..5` во всех трёх одинаков (ячейка 20).

| `driving_mode` ноутбука | ток | компонента | стенка | пара | $k_y$ (`use_phase_y`) | $K$ в профиле | где активен | `excitation_mode` C++ |
|---|---|---|---|---|---|---|---|---|
| `_0 'je_z, k_y=0, TM'` | $j_e$ | $z$ | левая (металл) | TM $(E_z,H_y)$ | 0 (`False`) | `sign_K_E_l*K_E_left_conductor` | профиль только в mode3 (новая строка), к уравнениям не подключён | `je_z_left` (новый) |
| `_1 'je_y, k_y=0, TE'` | $j_e$ | $y$ | левая | TE $(E_y,H_z)$ | 0 | `sign_K_H_l*K_H_left_conductor` | **`MenDrive_qnm.ipynb`** (активный) | `je_y_left` = `electric_left` |
| `_2 'je_z, TM'` | $j_e$ | $z$ | левая | TM | $\ne0$ (`True`) | `sign_K_E_l*K_E_left_conductor` | **`…mode2_je_z_TM.ipynb`** (активный) | `je_z_left` (в 1D $k_y=0$) |
| `_3 'jm_z, TE'` | $j_m$ | $z$ | правая (феррит) | TE | $\ne0$ | `sign_K_H_r*K_H_right_conductor` | **`…mode3_jm_z_TE.ipynb`** (активный) | `jm_z_right` = `magnetic_right` |
| `_4 'jm_z, k_y=0, TE'` | $j_m$ | $z$ | правая | TE | 0 | `sign_K_H_r*K_H_right_conductor` | только mode3, новая строка, не подключён | `jm_z_right` |
| `_5 'jm_y, k_y=0, TM'` | $j_m$ | $y$ | правая | TM | 0 | `sign_K_E_r*K_E_right_conductor` | только mode3, новая строка, не подключён | `jm_y_right` (новый) |
| — | $j_e$ | $y$ / $z$ | правая | TE / TM | — | — | нет | `je_y_right` / `je_z_right` |
| — | $j_m$ | $z$ / $y$ | левая | TE / TM | — | — | нет | `jm_z_left` / `jm_y_left` |

**Важно.** В исходном коде комментарий к `excitation_mode=1` («electric_left — MenDrive_qnm-mode2_je_z_TM»)
неточен: этот режим всегда возбуждал $j_{e,y}\to E_y$, то есть TE (`driving_mode_1` из `MenDrive_qnm.ipynb`),
а не $j_{e,z}$ (TM, mode2). Источника $j_{e,z}$ и $j_{m,y}$ в ядре не было. Теперь они есть: `je_z_left`,
`jm_y_right` и т. д. Старые имена и коды 0/1 сохранены и считают те же числа.

Строго сравнимы с 1D-FDTD моды с $k_y=0$ (`_0`, `_1`, `_4`, `_5`) при $k_z^2\ll\omega^2/c^2$. В ноутбуках
$k_z=\pi/L_z$ (`resonator_type="half"`, $L_z=100$ см), в 1D $k_z=0$. Для `_2`, `_3`
($k_y=(2m+1)\pi/(2b)\ne0$, `phase_y`) профиль по $x$ сравним при фиксированном $y=y_0$ только приближённо,
потому что $K_{vac}=\sqrt{\omega^2/c^2-k_y^2-k_z^2}$ отличается от $\omega/c$.

### 8.3. Параметр `polarization` (реализация TODO в `reset_state()`)

| значение | что обновляется в `step_once` |
|---|---|
| `auto` (по умолчанию) | только пара источника: `je_y`/`jm_z` → $(E_y,H_z)$, `je_z`/`jm_y` → $(E_z,H_y)$ |
| `TE`, `TM` | указанная пара (если она не совпадает с парой источника, печатается предупреждение: поля = 0) |
| `both` | обе пары (как в исходном ядре) |

* **LLG / Hybrid**: намагниченность $\vec M$ (Родригес + Гильберт) входит в $H_y$ и $H_z$, то есть
  гиротропия $\mu_{yz},\mu_{zy}$. Поэтому всегда обновляются **обе** пары, при явном `TE`/`TM` печатается
  `ВНИМАНИЕ`. В линейном приближении при `bias_orientation=x` обе компоненты $H_y,H_z\perp\vec H_0$ связаны.
  При `y` или `z` линейно активна одна пара, а вторая возбуждается только нелинейностью.
* **JA / Preisach / Preisach2D** — скалярные модели вдоль $z$: действуют только на $H_z$. $H_y$ в феррите
  линеен с `mu_yy_lin` (по умолчанию 1, как в исходном ядре). Для TM-режимов $H_z\equiv0$, поэтому
  неявный шаг феррита (самая дорогая часть) **пропускается целиком**.
* `reset_state()` по-прежнему выделяет и обнуляет все четыре массива (их читают $T_{xx}$, энергия и формат
  чекпойнта, который не менялся), а также буферы `*_new` и токи. Неактивная пара остаётся тождественно нулевой.
* Старый C API (`mendrive_create`, `mendrive_create_ex`) работает с `POL_BOTH` и прежним отображением
  кода (0 → magnetic_right, иначе electric_left). Поэтому `mendrive_ctypes.py` и ноутбук
  `MenDriveCpp_run.ipynb` получают те же числа. Новые функции: `mendrive_create_mode(...,
  polarization, src_k_phase, src_sigma_fix, eps_left, eps_ferrite, mu_yy_lin, mu_zz_lin)`,
  `mendrive_get_polarization`, `mendrive_get_mode_info`, `mendrive_get_N`, `mendrive_step`, `mendrive_get_fields`.

### 8.4. Новые параметры (все по умолчанию дают прежнее поведение)

| ключ | умолч. | смысл |
|---|---|---|
| `excitation_mode` | `magnetic_right` | + `jm_z_right, je_y_left, je_z_left, jm_y_right, je_y_right, je_z_right, jm_z_left, jm_y_left` |
| `ferrite_model` | `Preisach2D` | + `Linear`: $B=\mu H$ ($\mu_{yy}$=`mu_yy_lin`, $\mu_{zz}$=`mu_zz_lin`), магнитная проводимость `sigma_m_leak` |
| `polarization` | `auto` | `auto | TE | TM | both` (см. 8.3) |
| `src_k_phase` | 0 | $j=\text{amp}\cdot\text{ramp}\cdot e^{-\xi/\delta}\sin(\omega t-k\xi)$; комплексная амплитуда $\tilde j=i\,\text{amp}\,e^{iK\xi}$, $K=k+i/\delta$, т. е. $e^{+iK(x-a)}$ справа и $e^{-iK(x+a)}$ слева (анзацы ячеек 47/48). Чтобы воспроизвести скин-профиль ноутбука: `skin_depth`$=1/\mathrm{Im}K$, `src_k_phase`$=\mathrm{Re}K$ |
| `src_sigma_fix` | 0 | 1: член $4\pi\Delta t\,j_e$ делится на $(\varepsilon+2\pi\sigma\Delta t)$ — согласованная полунеявная схема (в исходном ядре этого деления нет, ошибка $O(\sigma\Delta t)$) |
| `eps_left`, `eps_ferrite` | 1 | $\varepsilon$ стенок (в ноутбуках $\varepsilon_r=16$) |
| `mu_yy_lin` | 1 | $\mu_{yy}$ феррита для $H_y$ (JA/Preisach/Preisach2D/Linear; в LLG/Hybrid игнорируется) |
| `mu_zz_lin` | 1 | $\mu_{zz}$ феррита для $H_z$ (только `Linear`) |

Новые параметры с kind=2 входят в подпись чекпойнта, **только если отличаются от умолчания**. Поэтому
`--resume` старых (уже идущих) прогонов работает. Формат `checkpoint.bin` не менялся.

## 9. Экспорт полей (`// FIELDS:`)

### 9.1. Параметры

| ключ | умолч. | смысл |
|---|---|---|
| `fields_snap_phases` | пусто | доли периода, напр. `0,0.25,0.5,0.75` → `fields_snap.csv` |
| `fields_snap_every_periods` | 0 | снимки каждые K периодов (0 — только последний полный период) |
| `fields_snap_from_period` | −1 | первый период снимков (−1 = `record_from_period`) |
| `fields_harm` | пусто | номера гармоник, напр. `1` или `0,1,2,3` → `fields_harm_n<k>.csv` |
| `fields_harm_periods` | 0 | накопление по последним K периодам (0 — окно средних, `last_frac`) |
| `time_convention` | `minus` | `minus`: $F=\mathrm{Re}[\tilde F e^{-i\omega t}]$, как в ноутбуках (`exp(-I*omega_t)`, ячейки 46/308/586); `plus`: $e^{+i\omega t}$ (комплексное сопряжение) |
| `fields_harm_normalize` | `none` | `source`: $\tilde F/j_0$, $j_0=i\cdot\text{amp}$ (`minus`) или $-i\cdot\text{amp}$ (`plus`) — комплексная амплитуда источника $\sin\omega t$; сравнимо с `je0_ext`/`jm0_ext` ноутбука |
| `fields_every_x` | 1 | прореживание узлов по $x$ |

Экспорт на динамику не влияет. Его состояние (накопители гармоник, размер `fields_snap.csv`) хранится
в отдельном `checkpoint_fields.bin`. При `--resume` с другими параметрами экспорта или со старым
чекпойнтом гармоники продолжаются, только если окно накопления ещё не началось.

### 9.2. Определения

Гармоники накапливаются на лету, ряд не хранится:

$$\tilde F_n(x)=\frac{c_n}{T_{rec}}\sum_k F(x,\tau_k)\,e^{\,s\,i n\omega\tau_k}\,\Delta t,\qquad c_0=1,\ c_{n\ge1}=2,\qquad s=\begin{cases}+1,&\texttt{minus}\\-1,&\texttt{plus}\end{cases}$$

что для `minus` есть $\tilde F=\frac{2}{T_{rec}}\int F e^{i\omega t}dt$ и $F(x,t)=\mathrm{Re}\sum_n\tilde F_n e^{-in\omega t}$.
Здесь $\tau_k$ — истинное время отсчёта по схеме Йи: $E$ — $t$, $H,M,B$ — $t+\Delta t/2$, токи — $t-\Delta t$.
Поправки учитываются точно.

Координаты: A-сетка $x_i=-(a+h_l)+i\,\Delta x$, начало в центре вакуумного зазора, стенки $x=-a$ и
$x=+a$, как в ноутбуках (вакуум $-a\le x\le a$). $x_{cm}=x\cdot$`length_unit_cm`. Величины B-сетки
($H$, $M$, $j_m$) интерполируются в узлы: $(F_{i-1}+F_i)/2$. `region` = −1 / 0 / 1 — левая стенка /
вакуум / феррит.

Величины (16): `Ey Ez Hy Hz Mx My Mz By Bz je_y je_z jm_y jm_z w Sx Txx`, где
$B=H+4\pi M$ (для JA/Preisach $M$ — скаляр вдоль $z$, `Mz`; в `timeseries.csv` он по историческим
причинам назван `My`; для `Linear` $M=(\mu-1)H/4\pi$),
$w=(\varepsilon E^2+\vec H\cdot\vec B)/8\pi$, $S_x=(E_yH_z-E_zH_y)/4\pi$,
$T_{xx}=(E_y^2+E_z^2+H_y^2+H_z^2)/8\pi$ (как `Txx_at()`: давление $-\sigma_{xx}$, в вакууме $=w$).
В квадратичных величинах $H$ центрирован по времени. Для них полезны $n=0$ (средние) и $n=2$.

### 9.3. Формат CSV

`fields_harm_n<k>.csv`:
```
# MenDrive fields_harm: n=1 omega_n=... omega0=... T_rec=... periods=... samples=... a=... h_l=... h_r=... amp=...
# time_convention=minus: F(x,t) = Re[F~(x) exp(-i n omega0 t)] ...
# normalize=...; excitation_mode=...; polarization=...; x: центр зазора, стенки x=-a,+a; x_cm=x*...
n,omega_n,i,x,x_cm,region,Ey_re,Ey_im,Ey_abs,Ey_arg,Ez_re,...,Txx_arg
```
`fields_snap.csv`: `period,phase_target,phase,t,i,x,x_cm,region,Ey,Ez,Hy,Hz,Mx,My,Mz,By,Bz,je_y,je_z,jm_y,jm_z,w,Sx,Txx`
(фактическая фаза `phase` отличается от `phase_target` не больше чем на $\Delta t/T$).
Линейные CSV (`linear1d.write_profiles_csv`, `export_sage_profiles`, `analytic_for_run`) имеют тот же
формат без `i, region`.

### 9.4. Команды

```sh
make && make test                                   # ctypes-сравнение + tests/test_fields_linear.py
./mendrive -c examples/linear_check_TM.cfg          # Linear-феррит, je_z_left, гармоника 1 и снимки
python3 linear1d.py out/linear_check_TM             # точный линейный профиль -> linear_analytic_n1.csv
python3 plot.py out/linear_check_TM --compare out/linear_check_TM/linear_analytic_n1.csv --fit none --norm none
./mendrive -c examples/mode_je_z_left_TM.cfg        # новая мода (JA, скан резонанса)
./mendrive -c examples/fields_export_preisach_magnetic.cfg && python3 plot.py out/fields_pr_mag
./mendrive -c examples/notebook_preisach_magnetic.cfg fields_harm=1 fields_snap_phases=0,0.25,0.5,0.75
```

`plot.py` строит в `out_dir/plots/` графики `fields_harm_n<k>_{E,H,B,M,j,dArg,energy}.png` в стиле ячеек
521–537 ноутбуков: Re, Im, Abs, Arg/π тангенциальных компонент по трём участкам (левая стенка,
вакуум, правая стенка), TM зелёным (`color="green"`), TE синим, плюс график «delta Arg H и E».
Кроме того строятся `fields_snap_<F>.png` (снимки по фазам) и `compare_n<k>_<F>.png` с `compare_n<k>.txt`.

## 10. Как сверять с линейной моделью

1. **Внутренняя проверка схемы.** Возьмите `ferrite_model=Linear` и сравните с точным решением
   `linear1d.solve_linear` (кусочно-однородные слои, источник $e^{iK\xi}$, $E=0$ на внешних границах, как в
   FDTD). Ожидаемая невязка без подбора фазы и нормы — 3–8 % при $N=80$ и 1.5–4 % при $N=160$; она убывает как
   $O(\Delta x)$ из-за ступенчатой границы. Это делает `make test-fields`; отдельно —
   `examples/linear_check_*.cfg`. Numpy-копия схемы `fdtd1d_ref.py` должна совпадать с `./mendrive`
   до ~1e-8 (это проверка C++ без компилятора).
2. **Профиль из ноутбука.** В конце `MenDrive_qnm*.ipynb` (после ячейки `subs_fields(...)`, ≈506) выполните
   ```python
   import sys; sys.path.insert(0, '/Users/oleksiy.drozdov/articles/mendrive_llg/MenDriveCpp/standalone')
   from linear1d import export_sage_profiles
   export_sage_profiles('linear_mode3.csv',
       {'Ey': (Ey_left_conductor_d, Ey_vacuum_d, Ey_right_conductor_d),
        'Hz': (Hz_left_conductor_d, Hz_vacuum_d, Hz_right_conductor_d)},
       A=A, h_l=h_conductor_l, h_r=h_conductor_r, y0=0)
   ```
   затем `python3 plot.py out_dir --compare linear_mode3.csv --xcol xa` (координата $x/a$, нормировка по
   max, общая фаза). Перевод единиц, если нужно совпадение и по абсолютной шкале: длина
   $L_0$=`length_unit_cm` см, $\omega_{code}=\omega L_0/c$, $\sigma_{code}=\sigma_{СГС}L_0/c$,
   $\sigma_{m,code}=\sigma_{m}L_0/c$, `eps_*`/`mu_*_lin` — те же числа, что в ноутбуке.
   Для управляемого прогона используйте `omega_drive_subs` (вещественная $\omega$), а не QNM-частоту.
3. **Нелинейный режим.** При малом `amp` гармоника $n=1$ должна стремиться к линейному профилю с
   эффективной $\mu$ феррита. Рост `amp` даёт отличие $n=1$ и ненулевые $n=2,3$. Это и есть мера
   нелинейности, которую стоит показывать рядом с линейной кривой.
4. Ноутбук-помощник **`compare_linear_vs_cpp.ipynb`** содержит эти шаги в виде пар «Markdown-ячейка с
   формулами → Code-ячейка с численной проверкой (OK/FAIL и невязка)».

## 11. Замечания к незакоммиченным изменениям `MenDrive_qnm-mode3_jm_z_TE.ipynb`

(репозиторий `/Users/oleksiy.drozdov/articles`, HEAD `5fbe669`; рабочая копия = индекс, blob `ff584d7`
против `879dcc9` в HEAD, т. е. изменения проиндексированы, но не закоммичены)

* Ячейка 21: добавлены `jm_ext_y_op/jm_ext_y_sym` и профили для `driving_mode_0`
  (`je_ext_z`, TM, левый, $K_E$), `driving_mode_4` (`jm_ext_z`, TE, правый, $K_H$) и `driving_mode_5`
  (`jm_ext_y`, TM, правый, $K_E$). **Выбор $K_E$/`sign_K_E_r` для `jm_y` в TM верен**: в TM все компоненты в
  стенке, включая $H_y\propto\partial_xE_z$, идут с $K_E$ (анзац $E_z$, ячейки 47/48). `sign_K_H_r`/`K_H_right`
  здесь были бы ошибкой.
* **Вероятная опечатка (унаследована от `driving_mode_3` и скопирована в `_4`, `_5`, а также в ячейку 280):**
  профили правой стенки записаны как `exp(-I*sign_K*K_right*(x + a))`, а анзац полей правой стенки
  (ячейки 48/276/279) — `exp(+I*sign_K*K_right*(x - a)) = exp(-I*K_r*(a - x))`. При $\mathrm{Im}K>0$
  профиль $|e^{-iK(x+a)}|=e^{\mathrm{Im}K(x+a)}$ растёт вглубь феррита и в точке $x=a$ в $e^{2a\,\mathrm{Im}K}$
  раз больше единицы. Тогда $|f_{jm\_ext\_D}|\propto e^{\mathrm{Im}K(x+a)}e^{-\mathrm{Im}K(x-a)}=\text{const}$, и
  $\int_a^\infty$ в ячейке 332 не сходится. Шаблоны `collect(e^(2*ImK_H_r*a - 2*ImK_H_r*x))` в ячейке 330
  рассчитаны на форму $(a-x)$. Исправление: `(x + a)` → `(a - x)` в правых профилях ячеек 21 и 280.
* `# if ...` в ячейке 21 — это комментарии: все присваивания выполняются безусловно (`je_ext_z_profile` —
  дважды для `_0`/`_2`, `jm_ext_z_profile` — дважды для `_3`/`_4`; сейчас тексты совпадают, ошибки нет).
* Новые режимы `_0`, `_4`, `_5` **не подключены к уравнениям**: источники добавляются только для `_1`, `_2`, `_3`
  (ячейки 54/75), силы — тоже (327–332). Для `_5` нужно
  `eq_rot_E_right_conductor_y = ... == (I*omega/c*B_r[1])[0] - 4*pi/c*jm_ext_y_sym` и сила
  $f_x=-\frac1c\,j_{m,y}D_z^*/2$ (по формуле $f_x=-\frac1c(j_{m,y}D_z-j_{m,z}D_y)$ из ячейки 327).
* Ячейка 527 (Arg E): внутри цикла по `y_d` стоит `pl = plot(...)` вместо `pl += ...`, поэтому строится только
  последний `y_d`. В ветке `is_TE` нет подстановки `(x=x, y=y_d)` у `Ey_*_d`.
* `calc_F_jmD_conj_right(..., logging=True)` — умолчание изменено на `True`. Похоже на отладочный остаток.
