# Пакет B: проницаемость / сверхпроницаемость / растворимость — что литература говорит о W

**Вопрос.** W = (S_A/S_B)² = K1^A (K2^A)² / (K1^B (K2^B)²) — отношение «констант Сивертса, измеренных через грань A и через грань B». Детальное равновесие требует W = 1; при W ≠ 1 асимметричная мембрана давала бы стационарный поток при тепловом H2 одинаковых p и T с обеих сторон.

**Пакет (19 файлов, все прочитаны полностью; nishimura1991 — скан, все 7 страниц отрендерены и прочитаны с изображений; формулы, искажённые текстовым слоем, в hatano2007 ур.(3), hatano2008 (ДР-формула) и alimov2018 ур.(2) проверены по рендеру страниц).**
Номера страниц: для журнальных статей — номер страницы PDF (p.N) или журнальный номер страницы, если он указан явно (например, p.19685); для Bacal&Wada — страница рукописи.

Условные обозначения: α — вероятность диссоциативного прилипания с растворением (sticking), k_r — константа рекомбинационной десорбции (J = k_r c²), S (K, K_S) — константа Сивертса (c = S√p), Z — газокинетический множитель (поток = Z p).

---

## 1. Hatano et al. 2006 — `hydrogen_metal_iteraction/hatano2006.pdf`

**Ссылка:** Y. Hatano, A. Livshits, Y. Nakamura, A. Busnyuk, V. Alimov, C. Hiromi, N. Ohyabu, K. Watanabe, "Influence of oxygen and carbon on performance of superpermeable membranes", *Fusion Engineering and Design* **81** (2006) 771–776. doi:10.1016/j.fusengdes.2005.06.368

**Суть.** Технологическая работа о сверхпроницаемых мембранах Nb и V (0,1 мм) для откачки трития. Оксидные плёнки на V/Nb разлагаются выше ~700 K (кислород уходит в объём), остаётся хемосорбированный монослой O. Коэффициент прилипания атомов α_H не зависит от объёмной концентрации кислорода C_O, а α_H2 убывает с C_O (рис. 2, без чисел в тексте). Карбиды Nb/V проницаемы для H. Обсуждается ресурс в выхлопном газе термоядерного реактора.

| Величина | Значение | Ед. | Система | T | Стр. |
|---|---|---|---|---|---|
| Начало разложения оксида (O уходит в объём) | ~700 (почти полностью к 900) | K | V поликр., XPS | — | p.772 |
| Растворимость O в V | 3 | ат.% | V | 473 K | p.772 |
| Предел растворимости O в Nb | <1,5 | ат.% | Nb | <973 K | p.773 |
| Диапазон опытов по проницаемости | 473–1413 K; p_up(H2)=0,1–10 Pa | — | Nb, V 0,1 мм | — | p.773 |
| C_O (легирование) | до 1,5 (AES-образец 2) | ат.% | Nb, V | — | p.773 |
| α_H vs C_O | «completely independent» | — | Nb | — | p.773, Fig.2 |
| Скорость откачки атомарного H карбидизированной Nb-панелью | 10^20 | H m^-2 s^-1 | Nb+NbC/Nb2C | >473 K | p.774 |
| Коэф. прилипания O2 на Nb | ~0,01 / 0,2 | — | Nb | 745 K / 1428 K | p.775 |
| Коэф. реакции CH4 / C2H6 на чистом Nb | 2·10^-6 / 3·10^-5 | — | Nb | 900 K | p.775 |
| Макс. перепад давлений на мембране | 10 | Pa | Nb, V | — | p.773 |

**Цитаты.**
- p.771: "Suprathermal hydrogen particles can penetrate into the bulk through such adlayers of non-metallic impurities, while the reemission of hydrogen molecules is strongly impeded. Therefore, almost all hydrogen particles sticking to the surface permeate through the membrane, irrespective of its thickness and temperature. It is possible to separate the permeating hydrogen particles from any impurity gases including He and compress them at the downstream side [2]."
- p.773: "The sticking coefficient of hydrogen atoms, αH, was completely independent of the bulk oxygen concentration, CO, while αH2 decreased with increasing CO."
- О детальном равновесии, асимметрии при тепловом H2 с обеих сторон, равновесной растворимости — **ничего нет**. Сжатие идёт только для надтепловых частиц (атомы с горячих Ta-нитей).

**Релевантность для W: низкая** — работа о надтепловом (атомарном) потоке и деградации поверхности; тепловой H2 с двух сторон и константы Сивертса через разные грани не обсуждаются.


## 2. Hatano et al. 2007 — `hydrogen_metal_iteraction/hatano2007.pdf`

**Ссылка:** Y. Hatano, K. Watanabe, A. Livshits, A. Busnyuk, V. Alimov, Y. Nakamura, K. Hashizume, "Effects of bulk impurity concentration on the reactivity of metal surface: Sticking of hydrogen molecules and atoms to polycrystalline Nb containing oxygen", *J. Chem. Phys.* **127**, 204707 (2007). doi:10.1063/1.2804874

**Суть.** Трубчатая Nb-мембрана (0,1 мм, ~90% граней (100)) с контролируемым содержанием кислорода C_O = 0,03–1,5 ат.%. α_H2 определён прямо из стационарной проницаемости (α_H2 = 2J_perm/J_inc при одинаковых поверхностях) и независимо — через k_r из времени запаздывания и принцип детального равновесия α = z k_r/S². **Ключевой для W момент:** измерена проницаемость **в двух противоположных направлениях** (газ снаружи при T_g=300 K vs газ внутри трубки при T_g≈T_S) — получены близкие α_H2; согласие α из прямого процесса и из k_r через детальное равновесие подтверждено. α_H2 падает на порядки с ростом C_O (∝C_O^-2 при высоких T), α_H ≈ 0,25 от C_O и T_S не зависит. Идентичность двух поверхностей проверялась ТДС: потоки с внешней и внутренней поверхности «всегда близки».

| Величина | Значение | Ед. | Система | T | Стр./ур. |
|---|---|---|---|---|---|
| α_H2 чистой поверхности | ≈0,25, не зависит от T_S | — | Nb лента | 300–1800 K | p.4, p.6 (Fig.3) |
| α_H2 после распыления (остаточный O) | 10^-3–10^-2 | — | Nb, C_O≈0,03 ат.% | — | p.4 |
| α_H2 при T_S=500 K (ссылка 11) | 4·10^-5 | — | Nb + сегрегир. O | 500 K | p.3 |
| α_H2 при max C_O | 2·10^-5 | — | Nb, C_O max (≈0,2–1,5 ат.%) | T_S=1300 K, T_g=300 K; аналогично T=1350 K термостат. | p.11, p.12 |
| Макс. α_H2 в опытах с O | 5·10^-3 | — | Nb | — | p.8 |
| α_H | ≈0,25 (оценка ≥0,1; литература 0,1–0,4; для O/Nb 0,2–0,3) | — | Nb+O | 300–1400 K | p.8, Fig.5 |
| Зависимость α_H2(C_O) | ∝C_O^-2 при T_S>T_S*; ∝C_O^-1,5 при T_S<T_S* | — | Nb+O | — | p.7, p.10–11, ур.(7),(8) |
| 2E_b (сорбция на «дырках» покрытия) | 10±2 | kJ/mol | Nb+O | T_S<T_S* | p.11 |
| 2E_b^reg (на регулярных O-узлах) | >100 | kJ/mol | Nb+O | ~1350 K | p.11 |
| Кажущаяся E_a при T_S>T_S* | >110 | kJ/mol | Nb+O | 800–1300 K | p.11 |
| |ΔH_seg| кислорода | >50 (лит.: 71; 46–71; 50–63; 63) | kJ/mol | O–Nb | — | p.3, p.11 |
| Предэкспонента α*_H2 | 0,15–50 (T>T_S*); 10^-2→10^-4 (T<T_S*) | — | Nb+O | — | p.11, ур.(9) |
| T_S* (излом Аррениуса) | 800–1300 (рост с C_O); ~1000 при C_O 0,1–0,2 ат.% | K | Nb+O | — | p.9–10 |
| Лаг-время τ / L²/6D | 86 s / 0,15 s | s | Nb 0,1 мм | 820 K | p.5 (Fig.2) |
| D_O в Nb | 1,5·10^-6 exp(−117 kJ/mol/RT) | m²/s | Nb | — | p.10 |
| Максимальное покрытие O | ≈NbO0,9 | — | Nb, C_O=2 ат.% | — | p.8 |
| Коэф. аккомодации энергии β (лит.) | 0,37 | — | V(111)+O+S | — | p.12 |
| α_H2(v=1) | <2·10^-3 | — | Nb+O, C_O=0,2 ат.% | 1300 K | p.12 |
| Доля H2(v=1) в равновесном газе | ~0,01 | — | — | 1300 K | p.12 |
| p_up(H2) | 0,1 (низкий C_O) – 10 (высокий C_O) | Pa | — | — | p.5 |
| Формулы | α=2J∞/J_inc (2); α=z k_r/S² (3); J(t)=J∞ tanh²(t/τ), τ=L/(4√(k_r J∞)) (4); α_H=4J_H/J_inc,H (5) | | | | p.5–6 |

**Цитаты (дословно).**
- p.5 (стр. журнала 204707-4): "Sticking probabilities can be found from such desorption experiments by using the detailed balance principle (e.g., Ref. 13). In the present work, the integral sticking coefficient was directly determined from permeation because the dissociative sticking is one of the sequential steps of the permeation process. If the sticking is the rate-limiting step in permeation, and if the two membrane surfaces are identical, the permeation flux in steady state, J∞perm, equals to a half of the “sticking flux” Jinc αH2, and hence αH2 = 2J∞perm/Jinc. (2)"
- p.5: "If the gas is introduced at the outside of the tubular membrane (Fig. 1), Tg ≈ ambient temperature. Thus, with varying the sample temperature, one can find the dependence of αH2 only on TS. If, however, the gas is introduced into the inner membrane hollow (Fig. 1), one will find αH2(T) under thermostatic conditions (Tg ≈ TS ≡ T). Thus, relative effects of Tg and TS on αH2 can be found by comparing the permeability in the two opposite directions. Another way to find αH2(T) corresponding to thermostatic conditions is to determine αH2 from the rate constant of recombinative desorption kr with using the detailed balance principle, αH2 = z kr / S², (3) where z is the gas kinetic theory coefficient and S is the solubility constant of hydrogen in Nb."
- p.5: "Identity of the states of two surfaces was checked by thermal desorption experiments. The membrane sample was exposed to H2 at a certain TS, and then thermal desorption experiments were carried out by raising TS. Desorption fluxes from the outside and inside surfaces were always close to the other."
- p.8 (204707-7): "The data were obtained from permeation through a tubular membrane in two opposite directions (Sec. II C 1). The temperature dependencies of two corresponding curves are basically similar, and they give the absolute values of αH2 close to the other; the rather small difference in the absolute value can be attributed to the error in pressure measurements by different gauges. Such observations, first, give the grounds to equally refer the data on αH2 obtained at varying TS (Figs. 3–5) to the thermostatic conditions. Secondly, this means that αH2 is determined by TS and not dependent on Tg, at least in the first approximation. These two conclusions were also derived by another way (Sec. II C 1): the values of αH2 evaluated from the steady-state permeation rate at Tg = const = 300 K were close to those estimated for Tg = TS, from the permeation lag time τ [see Eqs. (3) and (4) and Fig. 2]."
- p.12 (204707-11): "If the activated sticking is determined by TS, the desorption flux will mostly consist of thermal molecules having a temperature TS (in accordance with the detailed balance principle). The above has been confirmed for the oxygen-covered Nb (Sec. III B) by the agreement of αH2 found from the direct process (steady-state permeation) with varying TS at Tg = const = 300 K with that evaluated from the desorption rate constant kr, with applying the detailed balance principle. Similar agreement was observed also for other metals covered by a monolayer of nonmetallic impurity.1,22,26,27 It is to be noted, however, that D2 flux desorbing from V(111) covered by oxygen and sulfur at TS = 950 K (Ref. 4) was observed to consist of more than a half of hyperthermal D2 molecules having the average translational energy exceeding 2kTS by a factor of 3–4 [...]"
- p.13 (204707-12, Conclusions): "The fact that αH2 is determined by the bulk impurity concentration provides us with a means to control the barrier of dissociative sticking (and, correspondingly, the barrier of recombinative release) of H2 molecules." — т.е. барьеры адсорбции и десорбции меняются **вместе** (микрообратимость).
- p.12: "H atoms are mainly absorbed via the regular surface sites covered by oxygen and not through the holes in coverage." — атомы и молекулы идут разными путями (важно для того, что α_H/α_H2 может сильно отличаться без нарушения ДР).

**Замечание для W.** Работа явно опирается на α = (z/S²)·k_r (так напечатано в ур.(3), проверено по рендеру страницы; здесь S, видимо, определена как обратная величина к c/√p — в Hatano 2008 та же формула записана как α = k_r K_S² √(2πmkT)) с **одной и той же** S для обеих граней; поверхности внутри и снаружи трубки сделаны одинаковыми (контроль ТДС), поэтому прямого теста W≠1 (асимметричные грани) здесь нет. Сравнение «двух направлений» — это проверка роли T_g vs T_S при одинаковых гранях; расхождение «rather small», списано на манометры. Числового значения S (константы Сивертса) в тексте нет (дана ссылка на Fromm & Gebhardt).

**Релевантность для W: высокая** — единственный в пакете прямой эксперимент «проницаемость в двух противоположных направлениях» + проверка детального равновесия α = z k_r/S² (согласие), хотя и при симметричных гранях.


## 3. Hatano et al. 2008 — `hydrogen_metal_iteraction/hatano2008.pdf`

**Ссылка:** Y. Hatano, A. Busnyuk, V. Alimov, A. Livshits, Y. Nakamura, M. Matsuyama, "Influence of Oxygen on Permeation of Hydrogen Isotopes Through Group 5 Metals", *Fusion Science and Technology* **54**:2 (2008) 526–529. doi:10.13182/FST08-A1869

**Суть.** Влияние объёмного кислорода на поверхностные константы H на V (проницаемость, мембрана 0,1 мм, 180–1200 °C, p_up = 10^-3–5 Pa) и Ta (ТДС, фольга 0,02 мм), сравнение с Nb. Вероятность проникновения I_P/I_inc резко падает с C_O; чувствительность к O растёт с атомным номером (Ta > Nb > V). k_r Ta определён из 2I_des = k_r C_H², α_H2 пересчитан из k_r **по принципу детального равновесия** с табличной K_S. Оценка для бланкета Pb-17Li.

| Величина | Значение | Ед. | Система | T | Стр./ур. |
|---|---|---|---|---|---|
| Связь проницаемости и прилипания при одинаковых гранях | I_P/I_inc = ½ α_H2(C_O,T) | — | все | — | p.527, ур.(1) |
| Детальное равновесие | α_H2 = k_r K_S² √(2πmkT) | — | Ta (K_S из Fromm&Gebhardt) | — | p.527 |
| Десорбция | 2I_des = k_r C_H² (2-й порядок) | — | Ta | до 850 °C | p.527 |
| C_O, равновесная с Pb-17Li | 2,4 / 0,054 / 0,034 | ат.% | V / Nb / Ta | 700 °C | p.528 |
| I_P/I_inc при этих C_O | 1,3·10^-4 / 5,5·10^-4 / 3,7·10^-5 | — | V / Nb / Ta | 700 °C | p.528 |
| Диапазон I_P/I_inc на рис.1–2 | 10^-6–10^-2 | — | V, Nb, Ta | 200–1000 °C | p.527, Fig.1–2 |
| ΔS_sol, ΔH_sol кислорода | −101 J/(K·mol O), −423 kJ/mol O | | V | — | p.528 |
| то же | −96, −387 | | Nb | — | p.528 |
| то же | −96, −383 | | Ta | — | p.528 |
| Толщина трубок бланкета | 0,5 | mm | — | — | p.528 |

**Цитаты.**
- p.527: "If the properties of upstream and downstream surfaces of a tube are identical to each other, the permeation probability of H2 is described as φP/φinc = ½ αH2(CO, T), (1) where αH2(CO, T) is the sticking coefficient of H2 molecules [...]. According to the principle of detailed balance6, αH2(CO, T) = kr(CO, T) K_S² √(2πmkT) where KS is the solubility of hydrogen in Ta [...]. The value of KS reported in Ref. 9 was employed for the evaluation."
- p.527: "Under all conditions examined, the permeation rate (V) and desorption rate (Ta) were controlled by surface process and not by bulk diffusion process."
- p.528: "The pressures of T2 on both upstream and downstream sides of tubes could be changed with various parameters such as flow rate of Pb-17Li." (асимметрия здесь только технологическая — upstream к Pb-17Li, downstream к вакууму с H2O.)
- Об асимметрии констант Сивертса через разные поверхности — ничего; K_S принята табличной (объёмной) величиной, не зависящей от состояния поверхности, хотя k_r и α меняются на порядки.

**Релевантность для W: средняя** — явное использование ДР α = k_r K_S² √(2πmkT) с одной объёмной K_S при поверхностных константах, меняющихся на порядки (т.е. неявно W=1: K_S ≡ k_d/k_r не зависит от поверхности), но прямого теста нет.


## 4. Cherkez et al. 2020 — `hydrogen_metal_iteraction/cherkez2020.pdf`

**Ссылка:** D.I. Cherkez, A.V. Spitsyn, A.V. Golubeva, V.M. Chernov, "Deuterium permeation through the low-activated V–4Cr–4Ti alloy under plasma irradiation", *Nuclear Materials and Energy* **23** (2020) 100756. doi:10.1016/j.nme.2020.100756

**Суть.** Плазменно-стимулированная (PDP) и газовая (GDP) проницаемость D через мембрану V–4Cr–4Ti (100 мкм, установка ПИМ, ЭЦР-плазма, ~0,1 Pa D2). PDP-потоки на 1–2 порядка выше GDP; вероятность проникновения χ = 0,03–0,12 при 600–825 K. Гистерезис по температуре (локальный максимум при ~720 K на подъёме, поток ×2 выше, чем на спуске) объясняется состоянием оксидного слоя (TiOx) на поверхностях. Приведена таблица режимов RR/RD/DR/DD с асимметрией k_f ≪ k_b.

| Величина | Значение | Ед. | Система | T | Стр. |
|---|---|---|---|---|---|
| χ (вероятность проникновения PDP) | 0,03–0,12 | — | V–4Cr–4Ti, 100 мкм | 600–825 K | p.4–6 |
| Аппроксимации χ(T) | 2,9·exp(−23284/RT); 11,5·exp(−27268/RT) | — (E в J/mol) | V–4Cr–4Ti | — | p.5, Fig.5 |
| E_χ | 23 и 27 | kJ/mol | V–4Cr–4Ti | — | p.5–6 |
| Падающий поток атомов D | 9,24·10^19 (I=0,5 mA/cm², n=2,66, +30% нейтралов) | at.D m^-2 s^-1 | ПИМ | — | p.4 |
| Ионный состав H3+:H2+:H+ | 72,5:21:6,5 | % | ПИМ | p=6,1·10^-2 Pa | p.3 |
| T_e, T_i плазмы | ~6 eV, ~0,4 eV | | | | p.3 |
| Давление D2 на входе | ~0,1 (0,12–0,18) | Pa | | | p.3–4 |
| Гистерезис PDP | J(↑T) ≥ 2·J(↓T), max ~720 K | — | V–4Cr–4Ti | | p.4 |
| Расчётная χ для чистого V (из α лит.) | 0,38 | — | V, 2,9 ат.% O | 963 K, BIAS −100 V | p.5 |
| α ионов H2+ в O-покрытом V (лит. [43]) | 0,24 → ~1 при 0,7→300 eV/H | — | V+2,9 ат.% O | 963 K | p.5 |
| α атомов H в V (лит. [44]) | 0,15–0,22, не зависит от потока 10^17–10^20 H m^-2 s^-1 | — | V | 343–623 K | p.5 |
| α атомов H в Nb (лит.) | ≈0,21 (0,25 в [44]); не зависит от C_O и T | — | Nb | | p.5 |
| Оценка χ для симметричной O-легированной Nb-мембраны | 0,11–0,16 | — | Nb | | p.5 |
| Сжатие на сверхпроницаемых мембранах (лит.) | 10^3 [14] – 10^4–10^5 [13] | — | Nb, V | | p.2 |
| Скорость откачки (лит.) | 0,5–2,5 L cm^-2 s^-1 (10^-3–10^-5 Pa); ~0,45 в JFT-2M при 0,03 Pa | | Nb, V | | p.2 |
| E_p GDP (лит. [51]) | 14,9 → 46,6 после отжига (сегрегация Ti) | kJ/mol | V–4Cr–4Ti | | p.5 |
| Содержание O | 0,06 | ат.% | V–4Cr–4Ti | | p.5 |
| Порог охрупчивания | H/M ≈ 0,2; для V: ≥4·10^4 Pa при 673 K, ≥20 Pa при 423 K | | V | | p.1 |

**Цитаты.**
- p.2: "One can see, that only in two modes permeating flux is a linearly dependent on incident particle flux – RR and DD modes. [...] RR should be the most ordinary permeation mode for such materials, and one can see that if kf = kb (symmetric membrane case) Jp max could be up to 50% of the incident particles flux, while in case of asymmetric membrane kf ≪ kb permeating flux could be practically 100% of the incident flux (if α ~ 1). In order to obtain such a high ratio of the permeating hydrogen flux to the incident one in the experiments, it's necessary to provide the following conditions: inlet surface should have a barrier transparent for suprathermal particles and reducing the inverse gas flux, while the outlet surface should not have such a barrier."
- p.2: "In some cases it can adversely affect the symmetry of the superpermeable membrane and reduce system performance (e.g. W atomizer and Nb membrane)."
- p.5: "The authors assumed a symmetrical membrane in RR permeation regime, and equilibrium state between dissolved oxygen in the specimen bulk and it presence at the surfaces as a monolayer."
- p.5: "Thus, it seems, that rather pronounced temperature dependence of permeation probability χ (Eχ = 23 and 27 kJ/mol) could characterize outlet surface state (probably titanium and/or vanadium oxide layers and/or processes on it)."
- О детальном равновесии, константе Сивертса, равновесной растворимости через разные поверхности — **ничего нет**. Асимметрия k_f ≠ k_b обсуждается только для надтеплового (плазменного) входного потока.

**Релевантность для W: низкая** — асимметричные мембраны рассматриваются лишь для надтепловых частиц (там она легитимна); о тепловом H2 с обеих сторон и о W не сказано ничего.


## 5. Skovoroda et al. 2002 — `hydrogen_metal_iteraction/skovoroda2002.pdf`

**Ссылка:** A.A. Skovoroda, V.S. Svishchov, A.V. Spitsyn, V.L. Stolyarov, Yu.M. Pustovoit, V.D. Borman, V.S. Kulikauskas, A.M. Shipilin, "Plasma-driven superpermeation of hydrogen through Nb membranes: bulk effects", *J. Nucl. Mater.* **306** (2002) 232–240. PII S0022-3115(02)01324-7

**Суть.** Плазменно-стимулированная сверхпроницаемость H через Nb-мембраны 0,025 и 0,1 мм (DECR-плазма, 0,01–0,1 Pa) при T < 1000 K. Наблюдается немонотонная T-зависимость эффективности η, гистерезис, обратная пропорциональность толщине (т.е. диффузионное ограничение, в противоречии с классической «поверхностной» формулой η = γ/(γ+1), γ = K_L/K_0). Объяснение — фазовый переход β→α в Nb–H (~470 K) и транспорт по границам зёрен («каналы»). Отдельно (установка In2) — сорбция/десорбция молекулярного H2: при p > ~3·10^-4 Pa сорбция, при меньших — десорбция; обратимость при 400–800 K, необратимость при ~1000 K. Максимальная эффективность 25%.

| Величина | Значение | Ед. | Система | T | Стр./ур. |
|---|---|---|---|---|---|
| Эффективность сверхпроницаемости | η = γ/(γ+1), γ = K_L/K_0 | — | теория [1,2,7] | — | p.233, ур.(1) |
| D(H) в Nb | 5,9·10^-8 exp(−1300/T) | m² s^-1 | Nb | 180–850 K | p.233 |
| K (коэф. рекомбинации) | 4·10^-25 exp(−13000/T) | m⁴ s^-1 | Nb | 500–1000 K | p.233 |
| Закон Сивертса (форма) | j_H2 = K C² | — | Nb | — | p.232 |
| Порог появления PDP-сверхпроницаемости | ~570 | K | Nb 0,025 мм | — | p.237 |
| T_c фазового перехода β→α (акустика) | 470±20 | K | Nb–H | — | p.236 |
| Граница сорбция/десорбция H2 (In2) | ~3·10^-4 | Pa | Nb 0,025 мм | 400–800 K и >1000 K | p.237 |
| Поглощённое количество | ∝√P (при высоких p); ∝ толщине (×4 для 0,1 vs 0,025 мм) | — | Nb | 400–800 K | p.237–238 |
| Необратимость при высокой T | десорбируется только половина при охлаждении | — | Nb | ~1000 K | p.238 |
| Поток H0 / ионов / макс. проникающий / сорбция H2 | 0,8 / 2,5 / 0,8 / 0,05 | 10^20 atom m^-2 s^-1 | P=0,07 Pa, T_e=4,1 eV, j=17 A/m² | — | p.234, Table 2 |
| Оценка диффузионного потока | ~10^21 | atom m^-2 s^-1 | Nb 0,02 мм | 570 K | p.239 |
| Макс. эффективность PDP | 25 | % | Nb | — | p.239 |
| Деградация η | ~1 | % / час | Nb | — | p.237 |
| H в исходной фольге | ~9 | % | Nb | — | p.236 |
| Ta в Nb | 0,15 | % | — | — | p.235 |
| Давление на выходе | никогда >100 Pa (0,001–0,1 Pa в опытах) | Pa | — | — | p.235, p.237 |

**Цитаты.**
- p.232: "As is known [1,2], the permeation efficiency of equilibrium molecular hydrogen through vanadium-group (V, Nb, and Ta) metallic membranes is extremely low in comparison with that of the non-equilibrium hydrogen. [...] The physical reason of such a phenomenon, which is called 'superpermeation', consists in the combination of the large absorption and small desorption probabilities under conditions, such that the membrane bulk operates as a trap for gas particles."
- p.232: "Sivert's law, jH2 = KC², describes the release of hydrogen from the membrane volume. [...] The recombination coefficient K depends strongly on the surface cleanness and cannot be controlled with a high accuracy, especially for Nb membranes [8]."
- p.233: "Classical relation (1) shows that the permeation through a thin Nb membrane should possess the following features: (i) superpermeation occurs for γ ~ 1 or higher, i.e., when the SB energy on the inlet side is close to or higher than that on the outlet side; (ii) the superpermeation efficiency η does not depend on the membrane thickness, temperature, particle flux density, and diffusion coefficient; and (iii) the volumetric membrane properties do not directly affect the superpermeation efficiency."
- p.233: "Note that the theory [1,2] predicts the disappearance of superpermeation not only for thick impurity layers, but also for an ideally cleaned inlet surface (see (1))."
- p.237: "The reversibility of sorption and desorption processes is observed in experiments with a preliminary outgassed membrane. The amount of the absorbed gas increases as √P at high pressures."
- p.238: "There is no reversibility of the sorption and desorption processes in this case. Only one-half of the previously absorbed hydrogen is desorbed at cooling." (при ~1000 K)
- p.238 (Discussion): "1. The permeation efficiency of equilibrium hydrogen through thin Nb membranes is very low in comparison with the efficiency of plasma-driven permeation."
- О детальном равновесии и константе Сивертса через разные грани — **ничего нет**; асимметрия K_L/K_0 рассматривается только для неравновесного входного потока.

**Релевантность для W: низкая/средняя** — показывает, что «асимметричная» мембрана (γ=K_L/K_0) даёт сжатие лишь при надтепловом входе, а «равновесный» H2 проходит слабо; прямого теста W нет. Гистерезис и необратимость при ~1000 K — признак нестационарных/фазовых эффектов, а не стационарного потока.


## 6. Alimov et al. 2014 — `hydrogen_metal_iteraction/alimov2014.pdf`

**Ссылка:** V.N. Alimov, A.O. Busnyuk, M.E. Notkin, E.Yu. Peredistov, A.I. Livshits, "Substitutional V–Pd alloys for the membranes permeable to hydrogen: Hydrogen solubility at 150–400 °C", *Int. J. Hydrogen Energy* **39** (2014) 19682–19690. doi:10.1016/j.ijhydene.2014.09.122

**Суть.** Равновесные PCT-изотермы (метод Сивертса, 10–10^6 Pa, 150–400 °C, H/M от 4·10^-4 до 0,6) для сплавов V–xPd (x = 0; 5; 7,3; 9,7; 12,3; 18,8 ат.%), обе поверхности лент покрыты 2 мкм Pd «для ускорения диссоциации/ассоциации». Pd резко снижает растворимость (за счёт энтальпии при почти постоянном энтропийном множителе S*); закон Сивертса выполняется до H/M ≈ 0,3 при x ≥ 9,7 ат.%. Это **объёмная равновесная** растворимость; влияние поверхности на её значение не исследовалось (Pd-покрытие — лишь катализатор, поправка на H в Pd вычтена).

| Величина | Значение | Ед. | Система | T | Стр./ур./рис. |
|---|---|---|---|---|---|
| Закон Сивертса | C = S√P | (H/M)·Pa^-0,5 | V, V–xPd | 150–400 °C | p.19684–85, ур.(1),(2) |
| S(T) | S = S*_V exp(ΔH/RT), S*_V-xPd ≈ S*_V | | V–xPd | | p.19686, ур.(3) |
| Общая точка пересечения при 1/T→0 (S*) | ≈2·10^-6 (по рис.3; Pd ≈6·10^-6) | (H/M)·Pa^-0,5 | V и V–xPd | экстрапол. | Fig.3 (оценка по графику) |
| ΔH растворения H | ≈ −32,5 (V) → ≈ −16 (V–18,8Pd), почти линейно по x | kJ/mol | V–xPd | 150–400 °C | Fig.4 (по графику) |
| ΔH для V–Cr (лит.) | ≈ −33 → −26,5 (x=0→30 ат.%) | kJ/mol | V–xCr [18] | | Fig.4 (по графику) |
| E_H-alloy | 2,55–2,4 (меняется <15%) | eV | V–xPd | | p.19688, Fig.4 |
| E_H2 (энергия связи H в H2) | 2,25 | eV/atom | — | | p.19688 |
| S_V-10Pd/S_V | ≈0,13 | — | V–Pd | 300 °C | p.19687 |
| S_Nb-10Pd/S_Nb (лит. [11]) | ≈0,30 | — | Nb–Pd | 300 °C | p.19687 |
| Граница закона Сивертса | H/M ≈ 0,05 (V); до ≈0,3 (x ≥ 9,7) | — | | | p.19684–85 |
| (H/M)_max | ≈0,7–0,8 | — | V | | p.19687 |
| Давление при H/M=0,2 | ≈0,017 | MPa | V | 350 °C | p.19689 |
| Рост равновесного давления при замене V→V–9,7Pd | ×65 (H/M=0,2), ×30 (H/M<0,02) | — | | 350 °C | p.19689 |
| C(V) при 0,1 / 0,5 MPa | >0,3 / ≈0,5 | H/M | V | 400 °C | p.19683 |
| Переход хрупкость/пластичность (лит.) | H/M = 0,22 (P > 0,05 MPa) | | V | 400 °C | p.19683 |
| Растворимость V ≈ ½ Nb | ×2 | | | 300 °C | p.19688 |
| T разрушения лент в H2 | 150 (x=0); 200–250 (5); 150 (9,5); 250–300 (12,6); 150–200 (18,8) | °C | V–xPd | | Table 1 |
| Время дегазации | ≈15 | min | Pd-покрытые ленты | 200 °C | p.19684 |

**Цитаты.**
- p.19683: "Both surfaces of each specimen were electrolessly plated with 2 μm layer of Pd to promote dissociative/associative processes at the gas–metal interface [2,4]."
- p.19684: "Since the equilibrium content of hydrogen was measured for the whole specimen including the Pd coating, we subtracted the portion of H absorbed in Pd when the hydrogen concentration in the alloys was computed from experimental data [...]. This correction was not large: a few percent for the alloy with the maximum content of Pd (V-18.8Pd)".
- p.19684: "until the concentration of solute hydrogen, C, remains lower than H/M ≈ 0.05, the P–C curves for pure V well obey Siverts' law С = S_V √P (1) where S_V is Siverts' constant for pure V."
- p.19684: "The solubility of hydrogen in pure V observed in the present experiment was close to the reference data [14,15]" (данные с Pd-покрытием совпадают с литературными для V, полученными другими методами/поверхностями — косвенно в пользу независимости S от поверхности).
- p.19689: "[...] formation of V-xPd alloys in the range of 0 < x < 18.8 at% (bcc solid solutions) hardly can result in reduction of global permeance of the composite membrane because hydrogen solubility in the V–Pd alloy remains not lower than that in the coating of pure Pd."
- О детальном равновесии, направлении проницаемости, асимметрии — **ничего**.

**Релевантность для W: низкая/средняя** — даёт эталонные равновесные константы Сивертса объёма (≡ K1·K2² в терминах задачи); совпадение S(V) для Pd-покрытых образцов с литературой — слабый косвенный аргумент, что поверхность (катализатор) не меняет равновесную S. Прямого сравнения граней нет.


## 7. Alimov et al. 2018 — `hydrogen_metal_iteraction/alimov2018.pdf`

**Ссылка:** V.N. Alimov, I.V. Bobylev, A.O. Busnyuk, M.E. Notkin, E.Yu. Peredistov, A.I. Livshits, "Hydrogen transport through the tubular membranes of V-Pd alloys: Permeation, diffusion, surface processes and WGS mixture test of membrane assembly", *J. Membrane Science* **549** (2018) 428–437. doi:10.1016/j.memsci.2017.12.017

**Суть.** Трубчатые мембраны Pd–(V–κPd)–Pd (κ = 0–13,3 ат.%, стенка 145–244 мкм, Pd 1–2 мкм с обеих сторон), 300–500 °C, 0,001–1,2 MPa. Приведено общее уравнение стационарного потока через многослойную мембрану с **разными** α_in и α_out и константами Сивертса S_i каждого слоя (ур.1). Из опыта с чистым V найдено α(T) на Pd-покрытии (не-UHV): α = 0,6 exp(−50200/RT), α(400 °C) = 7·10^-5. Получены D(T) в сплавах. Изотерма j(√P_in − √P_out) расщепляется на две ветви — признак поверхностного ограничения. Сборка из 18 трубок для 1 кВт ТЭ.

| Величина | Значение | Ед. | Система | T | Стр./ур. |
|---|---|---|---|---|---|
| Общее ур. потока | √(P_in − j/(Z α_in)) − √(P_out + j/(Z α_out)) − 2j Σ L_i/(S_i D_i) = 0 | — | многослойная мембрана | — | p.429, ур.(1) (форма по ур.(2), проверено по рендеру) |
| Предел поверхностного ограничения | j ≈ ½ α Z (P_in − P_out) | — | α_in=α_out | — | p.430, ур.(6) |
| Предел диффузионного ограничения | j ≈ 0,5(√P_in − √P_out)/(2L_Pd/S_PdD_Pd + L_V/S_VD_V) | — | | | p.429, ур.(4) |
| Z_H2 | ≈3,2/√T | mol H2 m^-2 s^-1 Pa^-1 | | | p.436 (Notations) |
| α (H2 на Pd-покрытии, не-UHV) | 0,6 exp(−50200 J/mol / RT) | — | Pd-V-Pd | 300–500 °C | p.432, ур.(10) |
| α(400 °C) | 7·10^-5 | — | Pd-V-Pd | 673 K | p.431–432 |
| α на атомарно чистом Pd | ≈1, не зависит от T | — | Pd поликр. | | p.432 |
| D_V | 3,5·10^-8 exp(−4807/RT) | m²/s | V | | p.431 |
| S_V | 0,24 exp(32646/RT) | mol H m^-3 Pa^-0,5 | V | | p.431 |
| S_Pd·D_Pd | 3,82·10^-7 exp(−14000/RT) | mol H m^-1 s^-1 Pa^-0,5 | Pd | | p.431 |
| D_V-13,3Pd | 2·10^-7 exp(−17900/RT) | m²/s | | 300–500 °C | p.433, ур.(11) |
| D_V-8,9Pd | 7,9·10^-8 exp(−11900/RT) | m²/s | | | p.434, ур.(12) |
| D_V-4,3Pd | 4,1·10^-8 exp(−7200/RT) | m²/s | | | p.434, ур.(13) |
| S_V-13,3Pd | 4,86 / 5,97 / 7,57 / 9,98 / 13,8 | mol H m^-3 Pa^-0,5 | | 500/450/400/350/300 °C | Table 1 |
| j / j_s / C_e,in / D (V-13,3Pd, P_in=0,58, P_out=0,1 MPa) | 0,049/8,2/0,029/1,18e-8 (500); 0,051/4,3/0,036/1,0e-8 (450); 0,052/2,4/0,045/0,81e-8 (400); 0,053/1,2/0,059/0,63e-8 (350); 0,049/0,55/0,082/0,43e-8 (300) | mol H2 m^-2 s^-1; —; H/M; m²/s | | °C | Table 1 |
| S_V-8,9Pd | 9,13/11,6/15,3/20,6/29,4 | mol H m^-3 Pa^-0,5 | | 550/500/450/400/350 °C | Table 2 |
| j / j_s / C / D (V-8,9Pd, P_in=0,1, P_out=0,0014 MPa) | 0,072/2,3/≤0,015/1,34e-8 (550); 0,076/1,7/0,026/1,16e-8 (500); 0,078/0,90/0,033/1,01e-8 (450); 0,081/0,50/0,044/0,90e-8 (0,80e-8*) (400); 0,076/0,25/0,06/0,83e-8 (350) | | | | Table 2 |
| S_V-4,3Pd | 22,6/30,9/44,3/67,2 | mol H m^-3 Pa^-0,5 | | 500/450/400/350 °C | Table 3 |
| j / j_s / C / D (V-4,3Pd, P_in=0,1, P_out=0,0022 MPa) | 0,16/1,70/<0,1/1,34e-8; 0,18/0,90/<0,1/1,20e-8; 0,19/0,50/0,1/1,09e-8; 0,19/0,25/0,15/1,06e-8 | | | 500/450/400/350 °C | Table 3 |
| C_e при 1,15 MPa | H/M = 0,52 | | V | 400 °C | p.431 |
| D_V/D_V-13,3Pd | <3 (300 °C), <2 (400 °C) | | | | p.434 |
| Параметры WGS | 1,2 MPa общ., p_H2 0,56 MPa; извлечение 80% (теор. max 84%); ~12 slpm; 18 трубок, 772 см² | | V-5,2Pd | 350–370 °C | p.435–436 |
| Ресурсные испытания | ≈8000 | h | V-Pd | | p.435 |

**Цитаты.**
- p.429: "With taking into account the surface processes of dissociative absorption of H2 molecules and associative desorption of absorbed H atoms, the correlation between the density of steady state hydrogen flux through a multilayer metallic membrane, j, and the inlet and outlet pressures of hydrogen (Pin and Pout respectively) is [7,8,19,27] [Eq. (1)] where αin is the probability of dissociative sticking of H2 molecule to the inlet surface followed by hydrogen solution in the metal lattice (sticking coefficient), αout is the same for the outlet surface, Li is the thickness of i-layer, Si and Di are the hydrogen solubility constant and diffusivity for the i-layer respectively, ZH2 is the proportionality factor connecting hydrogen pressure and the density of incident flux".
- p.429: "It was taken αin = αout = α in Eq. (2) since both sides of the membrane are identical and interact with pure hydrogen in our experiment."
- p.432: "The sticking coefficients to the inlet and outlet surfaces are assumed to be equal (αin = αout = α) since both sides of the membrane are identical (plated with Pd) and interact with pure hydrogen or nitrogen in our experiment."
- p.432: "Of course, the obtained α(T) is no means a fundamental property of H2-Pd system. For comparison, the kinetics of H2 molecule dissociative absorption through the atomically clean surface of polycrystalline Pd is radically different: α does not depend on T and close to 1 [...]. On the other hand, the magnitude and temperature dependence of obtained α(T) is typical for any non-ultrahigh vacuum conditions when the surface of Pd is inevitably covered by a film of nonmetal impurities (usually about 1 monolayer in thickness)".
- p.431: "two values of j400°C correspond to one and the same difference √Pin − √Pout depending on the range of Pin. Correspondingly, two branches of j400°C are formed (Fig. 2) whereas, if the surface processes are insignificant (Eq. (4)), only one value of permeation flux can correspond to given value of difference √Pin − √Pout and only single isotherm can be formed."
- p.433: "Three factors complicate the determination of D by the indicated way: (1) some doubts in a complete identity of surface of Pd coating for different membrane samples as well as in a complete stability of the surface state in the course of experiments, [...]"

**Замечание для W.** Модель ур.(1) — ровно кинетическая схема задачи: поверхность входит **только** через α_in и α_out, а растворимость S_i — свойство объёма слоя. При P_in = P_out = P ур.(1) имеет единственное решение j = 0 при **любых** α_in ≠ α_out, т.е. авторы молча принимают W = 1 (одна S на объём, не зависящая от грани). Это модельное допущение, а не измерение: асимметричные (α_in ≠ α_out) мембраны в работе экспериментально не исследовались. Pd-покрытие с обеих сторон делает S_Pd общей для обеих граней.

**Релевантность для W: средняя** — даёт явную стандартную модель с α_in ≠ α_out и одной S (W=1 «по построению»), плюс реальные α, S, D; прямой проверки W нет.


## 8. Suzuki & Yukawa 2021 — `hydrogen_metal_iteraction/suzuki2021.pdf`

**Ссылка:** A. Suzuki, H. Yukawa, "Quantitative Evaluations of Hydrogen Diffusivity in V-X (X = Cr, Al, Pd) Alloy Membranes Based on Hydrogen Chemical Potential", *Membranes* **11** (2021) 67. doi:10.3390/membranes11010067

**Суть.** PCT (473–673 K) и газовая проницаемость (573–673 K, перепад давлений, permeate 10 kPa) сплавов V–Cr, V–Al, V–Pd (диски 0,37–0,56 мм, PdAg 200 нм с обеих сторон). Поток описан через градиент **химического потенциала** водорода: J = (RTB/2L)∫c·d ln(P/P0)/dc dc = (RTB/2L)·f_PCT, где концентрации на гранях c1, c2 берутся из **равновесных** PCT-кривых по давлениям на сторонах. Найдена подвижность B(T), Q, B0; DFT: H не может занимать T1-позиции рядом с Al и Pd (блокировка). Поверхностная кинетика не рассматривается (предполагается равновесие газ–объём на каждой грани).

| Величина | Значение | Ед. | Система | T | Стр./ур. |
|---|---|---|---|---|---|
| Поток через хим. потенциал | J = (RTB/2L)∫_{c2}^{c1} c·d ln(P/P0)/dc dc | mol H m^-2 s^-1 | все | | p.2, ур.(1) |
| Вант-Гофф | ln(P(c)/P0) = 2ΔH(c)/RT − 2ΔS(c)/R | | | | p.7, ур.(3) |
| ΔS_0,2 (при H/M = 0,2) | ≈ −47 J/(mol H·K) (−94 J/(mol H2·K)); не зависит от типа/количества легирования | | V–X | | p.9 |
| S° H2 газа (сравнение) | ~130 | J/(mol H2·K) | | | p.9 |
| ΔE вставки H (DFT, V16H1) | ≈ −0,36 (эксп. ΔH = −0,32) | eV/H атом | V | | p.15 |
| ΔE для V15Cr1H1 / V15Fe1H1 | ≈ −0,25 / −0,30 | eV | | | p.15 |
| ΔE для V15Al1H1(T1), V15Pd1H1(T1) | ≈ 0 (T1 блокирована) | eV | | | p.15 |
| B (подвижность), 10^13 mol m² s^-1 J^-1 | V-4Cr 11,0/7,8/4,4; V-14Cr 14,3/11,4/7,5; V-23Cr 13,8/11,7/8,5; V-5,5Al 11,5/7,7/4,8; V-16Al 7,6/5,8/3,8; V-20Al 4,7/3,4/2,4; V-25Al 2,0/1,6/1,0; V-10Pd 5,7/4,2/3,0; V-14,5Pd 5,6/4,1/3,3; Pd 10,0/–/6,0; Pd-23Ag 8,3/6,4/4,5 | | | 673/623/573 K | Table 5 |
| B0 (10^11 mol m² s^-1 J^-1) / Q (kJ/mol) | V-4Cr 23,3/29,8; V-14Cr 5,9/20,6; V-23Cr 2,3/15,6; V-5,5Al 18,9/28,5; V-16Al 4,0/22,0; V-20Al 2,4/22,0; V-25Al 1,1/22,2; V-10Pd 2,4/21,0; V-14,5Pd 1,1/16,6; Pd 1,4/14,9 (423–673 K); Pd-23Ag 7,3/24,4 (373–673 K) | | | 573–673 K | Table 6 |
| Толщины мембран | 0,366–0,562 | mm | | | p.3 |
| Объёмное расширение | ~2,6 Å³ на 1 (H/M) | | V | | p.4 |
| DBTC | ~0,2 H/M | | V, Nb | | p.1 |
| Давления опытов | feed 20–1000 kPa, permeate 10 kPa | kPa | | 573–673 K | Tables 2–4 |

**Цитаты.**
- p.2: "This equation is based on the diffusion equation based on the chemical potential of hydrogen [18]. Here, J is hydrogen flux through the membrane, [...] c1 and c2 are hydrogen concentrations at feed and permeation sides of the membrane, P is hydrogen pressure, P0 is the standard hydrogen pressure (101325 Pa)".
- p.7: "The pressure of gaseous hydrogen in equilibrium with hydrogen atoms in metals and alloys with the hydrogen concentration of 0.2 (H/M), i.e., DBTC (P (0.2)), were quantified from the PCT curves".
- p.9: "ΔS0.2 was almost constant at approximately −47 J mol H−1 K−1 (−94 J mol H2−1 K−1) and independent of types and additive amount of alloy elements. The absolute value of the partial molar entropy was close to the entropy of hydrogen gas in the standard state (~130 J mol H2−1 K−1) [26], indicating that the partial molar entropy change is mainly attributable to the loss of hydrogen gas with high entropy."
- О детальном равновесии, асимметрии поверхностей, направлении проницаемости — **ничего**.

**Релевантность для W: низкая** — чисто объёмная (термодинамическая) модель: движущая сила — разность химпотенциалов, концентрации на гранях — равновесные по PCT; при равных давлениях J ≡ 0. Это неявное W=1, но без проверки.


## 9. Hanke et al. 2020 — `hydrogen_metal_iteraction/hanke2020.pdf`

**Ссылка:** S. Hanke, C. Day, T. Giegerich, J. Igitkhanov, Y. Kathage, X. Luo, S. Varoutis, A. Vazquez Cortes, T. Härtl, A. Busniuk, A. Livshits, S. Merli, A. Schulz, M. Walker, K. Baumgärtner, J. Hofmann, "Progress of the R&D programme to develop a metal foil pump for DEMO", *Fusion Engineering and Design* **161** (2020) 111890. doi:10.1016/j.fusengdes.2020.111890

**Суть.** Обзор программы KIT по металлофольговому насосу (MFP) для DEMO на основе сверхпроницаемости (установка HERMESplus, трубчатая Nb-мембрана 0,1 мм, 1 см × 10 см, Duo-Plasmaline 2,45 ГГц). Описаны 5 стадий процесса (возбуждение → падение → поверхность → диффузия → рекомбинация). Моделирование ProVak3D: скорость откачки 23,9 vs 13,4 м³/с при вероятности проникновения 20% vs 5%. Подчёркнуто, что сверхпроницаемость требует **надтепловых** частиц; газовая проницаемость будет измеряться при выключенной плазме для получения кинетики поверхности.

| Величина | Значение | Ед. | Система | T | Стр. |
|---|---|---|---|---|---|
| Площадь фольги HERMESplus | 30 | cm² | Nb трубка | | p.2 |
| T нити (атомизатора) | ~2000 | K | | | p.2 |
| Нижний предел давления collisional plasma | «well below 10» | Pa | Duo-Plasmaline | | p.2 |
| T_e, n_e плазмы (лит.) | ~3 eV, ~10^18 m^-3 | | | | p.2 |
| Модель MFP | 2 м длина, Ø фольги 0,5 м, ядро плазмы 0,1 м | | | | p.3 |
| Вероятность проникновения / девозбуждения / отражения | 5%/15%/80% и 20%/0%/80% | | | | p.3 |
| Скорость откачки | 23,9 и 13,4 | m³/s | | | p.3 |
| Фактор разделения | 82 → 96 | % | | | p.3 |
| Минимальная T против гидридов | >171 | °C | V, Nb | | p.4 |
| Достигнутые потоки проникновения | до 1 | Pa·m³ s^-1 m^-2 | HERMESplus | | p.4 |
| Магнитное поле в месте MFP | ~1 T; ЭЦР-порог 87,5 mT; стенд до 250 mT (2 kA) | | | | p.5 |

**Цитаты.**
- p.1: "Characteristic for superpermeation is that the permeation can work against pressure gradients. Superpermeation can be achieved with the production of suprathermal hydrogen particles, which impinge on a thin metallic membrane".
- p.1–2: "Eventually, it is released into the gas phase after finding its way to a surface and recombining with another hydrogen atom in step ⑤. This can happen on both sides of the membrane. The probability for downstream release depends on the ratio of the energy barriers on both sides. This probability is always several orders of magnitude larger than the permeation probability in the diffusion limited regime for low implantation depths. Once the hydrogen molecule is released on the downstream side as a thermal molecule, it is by several orders of magnitude less likely to be absorbed into the lattice again as compared to the initial suprathermal hydrogen [6,7]. This yields an effective pumping, even possible against large pressure gradients [8]."
- p.3: "It should be noted that even for slightly unfavorable asymmetric energy barriers on up- and downstream side, the plasma driven permeation still works and still is, by orders of magnitude, stronger than any gas driven permeation at similar conditions."
- p.4: "In this case, the permeation flux can be nearly half of the absorption flux (for symmetric membranes), or almost equal to it (higher inlet than outlet surface barrier [3])."
- p.4: "Thus, conditions for classical gas driven permeation are provided, in which the permeated flux (significantly lower than in plasma driven permeation conditions) is mainly given by absorption and desorption rates for the different surface states and membrane temperatures."
- p.5: "The lack of reliability of the permeation fluxes that has been noticed by several different authors [6,23,33], despite the same incident flux and membrane conditions, is the second major risk. This is probably due to the instability of the surface state with time under influence of external agents."
- О детальном равновесии и W — **ничего**. Перепад давлений работает только с надтепловым источником.

**Релевантность для W: низкая** — инженерный обзор; утверждения об асимметрии относятся только к надтепловому входному потоку (внешняя энергия плазмы/нити).


## 10. Zheludkevich et al. 2008 — `hydrogen_metal_iteraction/zheludkevich2008.pdf`

**Ссылка:** M.L. Zheludkevich, A.G. Gusakov, A.G. Voropaev, E.N. Kozyrski, S.A. Raspopov, A.A. Vecher, "“Superpermeability” and “pumping” of atomic hydrogen through palladium membranes", *J. Membrane Science* **320** (2008) 528–532. doi:10.1016/j.memsci.2008.04.045

**Суть.** Pd-мембрана 0,1 мм (Ø 2 см), атомарный H из СВЧ-разряда (~10% диссоциации), 293–523 K, поток ~10^17 (атомы) / 10^18 (молекулы) см^-2 с^-1. Наблюдены «сверхпроницаемость» и «накачка» (рост p(H2) в закрытой камере за мембраной до 10100 Pa при 2 Pa атомарного H, 363 K). Термодинамическая модель: движущая сила — высокий химпотенциал атомарного H (ΔG°_дисс = 406,6 кДж/моль); предельное давление накачки — из равенства прямого (атомного) и обратного (молекулярного) потоков. Максимум проницаемости атомов при 373 K.

| Величина | Значение | Ед. | Система | T | Стр./ур. |
|---|---|---|---|---|---|
| ΔG°_298 диссоциации H2 | 406,6 | kJ/mol | газ | 298 K | p.529, p.532 |
| Степень диссоциации | ~10 | % | разряд 2450 MHz, 60 W | | p.529 |
| Проницаемость молекул (эмпир.) | J = 7,29·10^-2 exp(−15900/RT)(√P′ − √P″)·S/Δx | Pa cm³ s^-1 (S в cm², Δx в cm) | Pd 0,1 мм | 293–523 K | p.531, ур.(8) |
| Потоки (Pa cm³/s) атомы 2 Pa / H2 15 Pa / H2 10^4 Pa | 348 K: 25,02/0,34/5,23; 373: 30,82/0,49/6,30; 393: 20,83/0,63/9,20; 423: 7,07/0,89/17,7; 473: 2,28/1,43/n/a | | Pd | | Table 1 |
| Расчётное P*_H2 / коэф. накачки | 348 K: 64184 Pa / 32092; 373: 46622/23311; 393: 12638/6319; 423: 730/365; 473: 29/15 | | Pd, P_H = 2 Pa | | Table 2 |
| Измеренное накопленное p(H2) | 40 (473 K), 405 (423 K); 10100 (363 K, мембрана разрушилась) | Pa | Pd | | p.530, p.532 |
| Вероятность проникновения атома (max) | ~3·10^-3 | — | Pd | 373 K | p.531, p.532 |
| Коэф. накачки (лит. [3]) | 10^3–10^7 | — | Pd | | p.528 |
| Падающий поток атомов в ранних работах | 10^11–10^12 | atoms cm^-2 s^-1 | | | p.528 |

**Цитаты.**
- p.530: "As well known the driving force for hydrogen transfer through the membrane is the difference in chemical potential of chemisorbed hydrogen at two surfaces of membrane and in turn the difference in concentration of hydrogen atoms dissolved in the metal at two sides of membrane."
- p.531: "When the membrane is exposed to molecular hydrogen, dependence of the permeating flow on pressure can be described with the following equation: Jmol = D K_S (√P′H2 − √P″H2) S/Δx (2) where P′H2 and P″H2 are the pressure values of molecular hydrogen at the front and rear surfaces of membrane, respectively and KS is the Siverts constant for solvation of hydrogen in palladium." (одна K_S для обеих сторон)
- p.531: "This phenomenon is also due to the high chemical potential of atomic hydrogen, which creates a higher concentration of dissolved hydrogen near the front surface of membrane. Diffusion of hydrogen through the membrane will continue until the system will reach equilibrium state where the dissolved gas concentrations at the two surfaces will become equal."
- p.531: "Equilibrium will be achieved when the flow caused by atomic hydrogen acting upon the front surface is equal to the oppositely directed flow created by molecular hydrogen acting upon the rear surface."
- p.531: "The reason for this is that atomic hydrogen has a much chemical potential than its molecular counterpart at the same pressure and temperature. In turn the difference in chemical potential on the two sides of membrane is a driving force for hydrogen permeation."

**Релевантность для W: средняя** — явно термодинамическая трактовка: «накачка» возможна только благодаря химпотенциалу атомарного H (внешний источник неравновесности), предельное давление P* связано с ΔG°_дисс; при чисто молекулярном H2 с обеих сторон движущей силы нет (одна K_S). Это аргумент за W=1, хотя и не проверка разных граней.


## 11. Nishimura, Komaki, Amano 1991 — `hydrogen_metal_iteraction/nishimura1991.pdf` (скан без текстового слоя; все 7 страниц отрендерены и прочитаны)

**Ссылка:** C. Nishimura, M. Komaki, M. Amano, "Hydrogen Permeation Characteristics of Vanadium-Nickel Alloys", *Materials Transactions, JIM* **32** (1991) No. 5, 501–507.

**Суть.** Газовая проницаемость Pd-покрытых (0,1 мкм электролитического Pd с обеих сторон) мембран V, V–10Ni, V–15Ni, V–25Ni (0,565–1,814 мм) при 423–673 K, 10^4–3·10^5 Pa; PCT (Сивертс, 473–673 K, 10²–3·10^6 Pa). Диффузионное ограничение (L²∝ t_lag, J·L ∝ √P, n ≈ 0,5, нет зависимости от толщины). **Ключевое для W:** константа растворимости, полученная из проницаемости K = Φ/D, в **1/10–4/10 раза меньше** равновесной K из PCT. Авторы объясняют это неустановлением равновесия на входной границе / уменьшением эффективной площади из-за оксидной плёнки под Pd, а для одиночной Pd-мембраны (лит. [10]) K совпадают «quite good».

| Величина | Значение | Ед. | Система | T | Стр./ур./рис. |
|---|---|---|---|---|---|
| Поток | J = D(C_u − C_d)/L = Φ(P_u^½ − P_d^½)/L | mol H2 m^-2 s^-1 | | | p.502, ур.(1),(2) |
| Сивертс | C = K P^½; Φ = DK | K: mol H2 m^-3 Pa^-½ | | | p.502, ур.(3),(4) |
| Лаг-время | t_l = L²/6D; композит по Ash et al. ур.(8); Φ_T по Barrie et al. ур.(9) | | | | p.502 |
| K_Pd/K_bulk (отношение констант растворения) | ≈0,1 | — | Pd/V–Ni | 423–673 K | p.502 |
| Совпадение D по ур.(8) и (7) / Φ по (9) и (2) | 0,5% / 0,2% | отн. ошибка | L2=0,565 мм, L1=0,1 мкм | | p.502 |
| Φ чистого V (Pd-покрытый) | ~10^-7 – 10^-8 (на 1–2 порядка выше экстраполяции Heinrich) | mol H2 m^-1 s^-1 Pa^-½ | V | 473–673 K | p.503, Fig.1 |
| Растрескивание V | при p_up ≈ 1·10^4 Pa | Pa | V | 473 K | p.501, p.503 |
| Без растрескивания | до 6·10^4 (V-10Ni); до 2·10^5 (V-15Ni) | Pa | | 473 K | p.501, p.504 |
| Растрескивание aged V-15Ni | 1·10^5 Pa | | | 413 K | p.504 |
| Φ(V–10Ni)/Φ(Pd) | 6–15 | — | | | p.501, p.503 |
| Φ(V–15Ni)/Φ(Pd) | 2–5 (≈5 при 473 K) | — | | | p.501, p.507 |
| Показатель n в J = BP^n | 0,51 (573 K), 0,43 (448 K) | — | V–15Ni | | p.503 |
| D0 / E_D | V-10Ni: 2,3·10^-7 m²/s / 20; V-15Ni: 1,3·10^-6 / 30; V-15Ni aged: 8,2·10^-7 / 27 | kJ/mol H | | 423–673 K | p.504, Table 1 |
| Чувствительность расходомера | 1·10^-10 | mol H2/s | | | p.503 |
| K из проницаемости / K из PCT | 1/10 – 4/10 | — | V–10Ni, V–15Ni | 423–673 K | p.506, Fig.8 |
| Разница H/M (K_perm vs K_eq) | ≈0,2 | H/M | | 573 K, 1·10^5 Pa | p.507 |
| H/M при 673 K, 10^5 Pa | V 0,306; V-10Ni 0,126; V-10Co 0,020*; V-10Fe 0,032*; V-10Cr 0,069*; V-10Ti 0,445* (*экстрапол. от 1000 K) | | | | p.506, Table 2 |
| ∂(H/M)/∂N_M | Ni −1,8; Co −2,9; Fe −2,7; Cr −2,4; Ti +1,4 | | | | Table 2 |
| Блокировка позиций одним Ni | 1,8 | позиций | V–Ni | 673 K, 10^5 Pa | p.501, p.506 |
| Линейное расширение при H | 3,3% (V), 1,3% (V-15Ni) | | | 573 K, 10^5 Pa | p.506 |
| Активация Pd | отжиг на воздухе 573 K, 600 s | | | | p.502 |
| Активация для PCT | 3 цикла гидрирования при 673 K до 3·10^6 Pa; дегидрирование 723 K 10,8 ks; равновесие — без изменений 1,2 ks | | | | p.502 |

**Цитаты (по рендеру страниц).**
- p.502: "From the formula of Barrie et al., the appropriate equation for hydrogen permeation through that type of composite membrane is given by 1/Φ_T = 2L1/((2L1+L2)D1K1) + L2/D2K2 [...] The ratio of solution constants (K1/K2) is approximately 0.1 in the examined temperature range. [...] Thus, the palladium overlayer functions just to prevent the membrane from oxidation."
- p.506: "Some difference is observed between the K values determined by two kinds of methods: the K values determined by permeation measurements are in the range 1/10 to 4/10 of those determined by the gas-metal equilibration. The difference suggests that the true equilibrium of hydrogen solution was not established at the upstream interface even after the steady state permeation rate was reached."
- p.506: "As shown in Fig. 8, the solubilities determined from Φ and D values were much smaller than those determined by the gas-alloy equilibration. It is in sharp contrast with the quite good agreement of K values for a single palladium membrane reported by Amano et al.(10) Recently Tahara et al.(37) proposed a new boundary condition for the hydrogen permeation through the multi-layered membrane. Their boundary condition was introduced for a non-equilibrium condition, however, the detailed mechanism is not described in their paper."
- p.506–507: "[...] if the solubility determined by the gas-alloy equilibration includes a contribution from trap sites, there can be a difference between the two kinds of solubilities. The solubility difference shown in Fig. 8, however, is too large to be ascribed only to the contribution from the trap sites, i.e., the difference in H/M ratio reaches about 0.2 at 573 K under a hydrogen pressure of 1×10^5 Pa."
- p.507: "Thus, the difference between the K values determined by the two kinds of methods can be qualitatively explained well by reduction in the effective area of the membrane caused by the oxide film."

**Замечание для W.** Это расхождение «S через поверхность при потоке» vs «S равновесная» (×2,5–10) — **не** измерение равновесной S через разные грани: K_perm = Φ/D получено в стационарном потоке и определяется входным граничным условием (неравновесная граница/эффективная площадь). В равновесии (без потока) PCT на тех же Pd-покрытых образцах даёт объёмную K. Направление проницаемости не менялось (только up→down при P_d ≈ 0).

**Релевантность для W: средняя** — единственный в пакете случай, где «константа Сивертса, измеренная через поверхность» количественно отличается от равновесной; но это кинетический (потоковый) эффект на входной границе, и сами авторы так его и трактуют — прямого довода за W≠1 нет.


## 12. Livshits, El Balghiti, Bacal 1994 — `hydrogen_metal_iteraction/livshits1994.pdf`

**Ссылка:** A.I. Livshits, F. El Balghiti, M. Bacal, "Dissociation of hydrogen molecules on metal filaments in H⁻ ion sources", *Plasma Sources Sci. Technol.* **3** (1994) 465–472.

**Суть.** Вероятность диссоциации H2 на горячих нитях (W, Ta, Nb, W–Th) K_s выведена **из детального равновесия** в термостате (газ ⇄ поверхность): K_s = (γ_s/2)·√(K(T)/P) (лимитирует десорбция атомов) или K_s = α_H2 (лимитирует адсорбция). Т.е. K_s = √2·γ_s·η_e — кинетический коэффициент поверхности × термодинамическая равновесная степень диссоциации η_e, **не зависящая от поверхности**. Для чистых Ta, Nb, W при эмиссионных T K_s = 0,2–0,3; W–Th — существенно меньше. Предложен метод измерения атомарного H по селективному поглощению Ta-лентой, покрытой монослоем неметалла (α_H = 0,1–1, α_H2 подавлен на порядки).

| Величина | Значение | Ед. | Система | T | Стр./ур. |
|---|---|---|---|---|---|
| Константа равновесия H2⇄2H | K(T) = 760·T^1,72·exp(−…/T) (показатель экспоненты в скане нечитаем) | Torr | газ | | p.466, ур.(4) |
| Детальное равновесие на стенке | K_s Z_H2 P_H2 = ½ γ_s Z_H P_H (ур.7) ⇒ K_s = (γ_s/2)(K/P)^½ (ур.12); K_s = √2 γ_s η_e (ур.14); K_s = α_H2 (ур.15) | | термостат | | p.466 |
| Z_H = √2 Z_H2 | | | | | p.466 |
| γ_s (рекомбинация H на чистой поверхности) | ≈1, не зависит от T | — | переходные металлы | | p.467, сноска [10] |
| γ_s на «нечистых» поверхностях (лит.) | 0,05–0,3 (Beitel: 0,05–0,2) | — | | | p.472, [10] |
| γ_s для атомов на металле с монослоем неметалла | 0,1–0,3, не зависит от T | — | | | p.472, [10] |
| α_H2 на чистых переходных металлах | 0,1–0,3, не зависит от T; W 0,3; Pt 0,24; Nb, Ta 0,2 | — | | | p.467 |
| K_s при эмиссионной T | 0,2–0,3, не зависит от T и P | — | W, Ta, Nb | >2100 K при P < 3·10^-3 Torr | p.467, ур.(18) |
| γ_s Ta при 1500–1700 K | ≈0,5 (по согласию с расчётом при K_s,sat = 0,2) | — | Ta | | p.470 |
| Растворимость C | W 0,2 ат.% (2580 K); Ta 2,5 (2470 K); Nb 2 (2450 K) | ат.% | | | p.467 |
| C_O в Ta | C_O = 1,35·10^-2 P_H2O exp(…/T) | ат.% (P в Torr) | Ta | | p.468, ур.(19) |
| θ_O/(1−θ_O) | 0,4 exp(…/T)·C_O | | Ta | | p.468, ур.(20) |
| P_H2O для θ_O/(1−θ_O)=1 | 4·10^-4 | Torr | Ta | 2470 K | p.468, ур.(21) |
| α_H (Ta с монослоем) | 0,1–1, не зависит от T металла; принято 0,2 | — | Ta | | p.468, p.470 |
| Ta-лента | 50 мкм, 420 см²; T_r 1300 → 700 K | | | | p.468 |
| Макс. n_H/n_H2 | 0,22 (T_H=T_H2) → 0,09 (T_H2=500 K, T_H=3000 K) | — | гибридный H⁻ источник | | p.471 |
| Остаточные газы | H2O 5·10^-?, CO+N2 1·10^-?, CH4 1,5·10^-? Torr (порядки в скане нечитаемы) | Torr | | Ta нить 2400 K | p.472, [28] |
| Давление опыта | ~3·10^-3 | Torr | | | p.471 |

**Цитаты.**
- p.465: "In order to determine the dissociation probability Ks of the H2 molecule in one collision with the surface let us examine the equilibrium. The gas is in a vessel (a thermostat), the wall of which is made of the material of interest at temperature sufficient to bring about noticeable [7] dissociation."
- p.466: "At equilibrium the dissociation and recombination rates are equal: Ks Z_H2 P_H2 = ½ γs Z_H P_H. (7)"
- p.466: "That is, the probability of dissociation Ks differs from the thermodynamically determined one for the given T and P value ηe only by the kinetic coefficient √2 γs, which has a clear physical meaning."
- p.466: "[...] the dissociation probability is determined by the values of P and T and by the kinetic coefficients γs and αH2, which are independent of P (and often also of T) and are related to the particular surface."
- p.472, сноска [7]: "A system with hot filaments, cold walls and a plasma is obviously a non-equilibrium one. However, under conditions of a small coverage or a hot surface, which is always the case, and also for the usually non-activated processes of atomic recombination and molecule sticking, the lack of equilibrium does not influence γs. In particular, this follows from analysis of the relation obtained by Brennan [8], and also from comparison with experimental data."
- p.468: "A monoatomic film of an impurity can reduce the absorption probability of hydrogen molecules, αH2, by orders of magnitude with respect to the case of a clean surface; as a result, molecular hydrogen practically stops being absorbed. In the meantime this film does not prevent absorption of atomic hydrogen: the thermal hydrogen atoms continue to be absorbed with probability αH = 0.1–1, independently of the metal temperature."

**Релевантность для W: средняя** — методологически прямой аналог: в термостате отношение «прямой/обратный» поток через любую поверхность фиксировано термодинамикой (η_e), а поверхность входит только кинетическим множителем; это та же логика, что требует W=1. Мембран и растворимости здесь нет.


## 13. Musyaev et al. 2005 — `hydrogen_metal_iteraction/ref/musyaev2005.pdf`

**Ссылка:** R.K. Musyaev, B.S. Lebedev, S.K. Grishechkin, A.A. Yukhimchuk, A.A. Busnyuk, M.E. Notkin, A.A. Samartsev, A.I. Livshits, "Tritium Superpermeability: Experimental Investigation and Simulation of Tritium Recirculation in “Prometheus” Setup", *Fusion Science and Technology* **48**:1 (2005) 35–38. doi:10.13182/FST05-A874

**Суть.** Цилиндрическая Nb-мембрана (Ø 100 мм, h 180 мм, стенка 0,1 мм) с Ta-атомизатором внутри (3 ленты, до 2300 K, 120 см² ≈ 20% площади мембраны). Продемонстрированы откачка H/D/T (скорость ∝ m^-½), отделение от He, сжатие и замкнутая рециркуляция трития. **Важно для W:** после остановки откачки выходной камеры давление на входе растёт из-за обратной *молекулярной* проницаемости; из отношения приращений давлений оценена «степень асимметрии мембраны» ~10³. При выключении атомизатора давления **возвращаются к исходным уровням**.

| Величина | Значение | Ед. | Система | T | Стр. |
|---|---|---|---|---|---|
| Скорость откачки мембраной | 2,5 (H) / 1,8 (D) / 1,5 (T) | l cm^-2 s^-1 | Nb 0,1 мм + Ta атомизатор | | p.35, p.37 |
| Отношение скоростей H/T, H/D | ≈1,7 и ≈1,4 (∝ √m) | — | | | p.37 |
| T атомизатора (макс.) | 2300 | K | Ta | | p.36 |
| T лент «в режиме ожидания» | 200–500 | °C | Ta | | p.36 |
| Площадь атомизатора | 120 cm² (~20% мембраны) | | | | p.36 |
| Степень асимметрии мембраны (inlet permeating flow / outlet flow) | ~10³ | — | Nb, пассивированная поверхность | | p.36 |
| Падение давления на входе при включении атомизатора | «few order of value» | | | | p.36 |
| Поток на входе в опыте по инвентарю | ~10^16 D/T s^-1 | | | | p.37 |
| Надтепловая энергия | ≥1 eV | | | | p.35 |

**Цитаты.**
- p.35: "The term “superpermeability” taken in this paper implies that metallic membranes of a macroscopic thickness may be virtually “transparent” to suprathermal (≥1 eV) hydrogen particles, comparable an opening of the same area. Virtually the whole flux of gas being implanted in metal lattice, or even just sticking to metal surface, passes through the membrane, irrespective of its temperature and thickness."
- p.36: "On switching off of the atomizer, the inlet and outlet pressures are returned to their starting levels."
- p.36: "[...] pumping of the outlet chamber with the external pump was ceased (point of time t2) to result in the outlet pressure rising due to the accumulation of permeating tritium at the outlet side. That is accompanied with an inlet pressure rise because of an enhanced molecular driven permeation through the membrane back from the outlet to the inlet side. Comparison of the pressure increment at the input and output side was used to calculate the degree of the membrane asymmetry (relation between the inlet permeating flow and the outlet flow), which depend on degree surface passivation and in this case was ~10³."

**Замечание для W.** Наблюдается (a) обратный *молекулярный* поток в сторону **низкого** давления (out→in), т.е. молекулярный H2 течёт «по градиенту», (b) без атомизатора исходные уровни давлений восстанавливаются — стационарного перепада при одном лишь тепловом H2 не возникает (в пределах этой качественной фразы; установка с откачкой обеих камер, прямого теста «равные давления с обеих сторон» нет). «Асимметрия ~10³» — кинетическая (отношение констант рекомбинации/прилипания граней), не W.

**Релевантность для W: средняя** — прямо фиксирует сильную кинетическую асимметрию граней (~10³) и при этом отсутствие «самонакачки» без надтеплового источника; количественного теста W нет.


## 14. Musyaev et al. 2008 — `hydrogen_metal_iteraction/ref/musyaev2008.pdf`

**Ссылка:** R.K. Musyaev, A.A. Yukhimchuk, B.S. Lebedev, A.O. Busnyuk, M.E. Notkin, A.A. Samartsev, A.I. Livshits, "Study of Hydrogen Isotopes Superpermeation Through Vanadium Membrane on “Prometheus” Setup", *Fusion Science and Technology* **54**:2 (2008) 523–525. doi:10.13182/FST08-A1868

**Суть.** То же, что [13], но мембрана из чистого V (Ø 100 мм, h 180 мм, 0,1 мм), Ta-атомизатор 1300–2100 K, p_in 10^-7–10^-4 mbar. Сверхпроницаемость тепловых атомов H/D/T, сжатие > 10³, замкнутая рециркуляция трития. **Ключевое для W:** из обратной *молекулярной* проницаемости (out→in) в режиме сжатия определены **разные коэффициенты прилипания H2 на двух гранях одной мембраны**: α_inlet = 2·10^-5 – 5·10^-5, α_outlet ~ 2·10^-4; причина — осаждение Ta с атомизатора на входную грань. При выключении атомизатора давления возвращаются к исходным.

| Величина | Значение | Ед. | Система | T | Стр. |
|---|---|---|---|---|---|
| α_H2 входной грани | 2·10^-5 – 5·10^-5 | — | V 0,1 мм, грань с осаждённым Ta | T мембраны не указана (нагрев излучением) | p.524 |
| α_H2 выходной грани | ~2·10^-4 | — | V, та же мембрана | | p.524 |
| Отношение α_out/α_in | ~4–10 (расчёт по данным статьи) | — | | | p.524 |
| Макс. сжатие p_out/p_in | > 10³ | — | V | | p.524, Fig.3 |
| Скорость откачки (макс.) | ~2,3 (H) / 1,6 (D) / 1,2 (T) | l cm^-2 s^-1 | V | | p.525 |
| Проводимости откачки камер | 18 (вход) и 15 (выход) | l/s для H2 | | | p.523 |
| T атомизатора | 1300–2100 (режим ожидания 500–700) | K | Ta | | p.524 |
| Давление на входе | 10^-7 – 10^-4 | mbar | | | p.524 |

**Цитаты.**
- p.524: "Reducing the atomizer temperature to its initial value (“switching off” of the atomizer) causes the inlet and outlet pressures to return to their initial levels."
- p.524: "The pumping of the outlet chamber with the external pump was ceased after permeating flux reached the steady state level. As a result, the outlet pressure increased (see Fig. 3) due to the accumulation of permeating hydrogen in the outlet chamber (compression mode). That was accompanied with an inlet pressure rise because of enhanced molecular driven permeation through the membrane back from the outlet to the inlet side. Comparison of the pressure increment in the input and output sides was used to calculate the important membrane parameter of maximum compression (ratio of outlet to inlet pressures). As one can see in Fig.3, the compression value exceeds three orders of magnitude."
- p.524: "The obtained data also demonstrated that the membrane behaves strongly asymmetric – a major part of atomic hydrogen absorbed by membrane at the inlet side is released to outlet side. Comparison of the pressure increment in the input and output sides permitted the calculation of the degree of the membrane asymmetry (relation between the inlet permeating flow and the outlet flow), which in turn depends on the sticking coefficients (experimental value: α_inlet = 2·10^-5 to 5·10^-5 and α_outlet ~ 2·10^-4). The probable cause of the asymmetry is deposition of tantalum onto the inlet membrane surface from atomizer."

**Замечание для W.** Это единственное в пакете *численное* измерение различия кинетики тепловых H2 на двух гранях (α различаются в 4–10 раз). Но различие α само по себе W не задаёт: W = 1 при любом α_A ≠ α_B, если k_r меняется согласованно (α ∝ S²k_r). Метод пересчёта «приращения давлений → α» в статье не описан явно; судя по формулировке, он опирается на стандартную модель с единой S (W=1), т.е. W не проверялся. Выключение атомизатора → возврат к исходным давлениям (качественно, при непрерывной откачке).

**Релевантность для W: средняя/высокая** — реальная асимметричная (по α ~×4–10) мембрана V с тепловым H2 на выходе; сжатие возникает только при работе атомизатора. Прямой проверки равенства S_A = S_B (W) нет.


## 15. Nakamura et al. 2000 — `hydrogen_metal_iteraction/ref/nakamura2000.pdf`

**Ссылка:** Y. Nakamura, S. Sengoku, Y. Nakahara, N. Suzuki, H. Suzuki, N. Ohyabu, A. Busnyuk, M. Notkin, A. Livshits, "Deuterium pumping experiment with superpermeable Nb membrane in JFT-2M tokamak", *J. Nucl. Mater.* **278** (2000) 312–319. PII S0022-3115(99)00243-3

**Суть.** Две трубчатые Nb-мембраны (Ø 1,5 см, L 9,4 см, стенка 0,2 мм) в закрытом диверторе JFT-2M; внутренняя (выходная) поверхность обработана дробеструйно для увеличения площади рекомбинации. Откачка D пропорциональна давлению в диверторе, макс. 7,3·10^19 D m^-2 s^-1 при 0,03 Pa. В стенде PMTD: MDP экспоненциально падает с T, PDP от T не зависит. **Важно для W:** асимметрия мембраны («отношение констант скорости десорбции с двух сторон») измерена ТДС после насыщения: поток десорбции с выходной стороны вдвое больше, чем с входной; при равенстве концентраций у обеих граней → A_s ≈ 2.

| Величина | Значение | Ед. | Система | T | Стр./ур. |
|---|---|---|---|---|---|
| Асимметрия A_s (k_des,down/k_des,up) | ≈2 (ТДС: поток вниз = 2× поток вверх) | — | Nb, выходная грань дробеструйная | нагрев после насыщения | p.316 |
| f_s = A_s/(1+A_s) | 0,67 | — | | | p.317, ур.(2) |
| Макс. плотность потока проникновения | 7,3·10^19 | D m^-2 s^-1 | | p_div = 0,03 Pa, T_M ≈ 650 °C | p.312, p.317 |
| Макс. скорость откачки | 2,8·10^17 | D/s | | | p.316 |
| Эффективная площадь | 38 | cm² (на трубку, ½ поверхности) | | | p.317 |
| f_g (геом. фактор) | ~0,3 | | | | p.317 |
| f_R (отражение D 2–5 eV) | 0,4–0,6 | | Nb + O | | p.317 |
| n(D атомов) в диверторе | 0,9–2,1·10^17 | m^-3 | | | p.318 |
| Поток атомов с поверхности дивертора | 5,6–8,4·10^20 | D m^-2 s^-1 | | | p.318 |
| Достигнутый ранее поток / теор. предел | 3·10^21 / ~10^23 | m^-2 s^-1 | | | p.313, p.318 |
| Устойчивость PDP к ионам | до 70 eV | | Nb | | p.313 |
| T мембраны | нагрев до 870 °C, при разряде ~650 °C; гидриды при T_M < 200 °C | | | | p.314 |
| Проводимость диафрагмы | 2,6 l/s для D2 | | | | p.314 |
| Плазма PMTD | n_e 10^16 m^-3, T_e несколько eV, p ~10^-2 Torr | | | | p.315 |

**Цитаты.**
- p.312: "The superpermeability originates from an activation barrier for the dissociative absorption of molecular hydrogen (and consequently also a barrier for the desorption of absorbed atoms through their recombination into molecules)."
- p.313: "The probabilities of desorption through the up- and downstream sides depend in this case only on the relative height of up- and downstream barriers. If the downstream barrier is lower than the upstream one (favorable membrane asymmetry), virtually the whole incoming flux may pass through the membrane with no dependence on its temperature and thickness."
- p.315–316: "We also investigated a membrane asymmetry (i.e. the ratio of desorption rate constants at both sides) using one more characteristic feature of superpermeation regime: the uniform distribution of absorbed hydrogen over the membrane thickness [2–4]. [...] Then the membrane was heated up again and the thermodesorption fluxes from both membrane surfaces were measured. The desorption flux from the downstream side was twice that from the upstream one. Taking into account the equality of hydrogen concentration at the both membrane sides, we come to the membrane asymmetry of around 2."
- p.315: "The molecular driven permeation (MDP) rate exponentially decreases with decreasing membrane temperature as it is usually observed, while the plasma driven permeation (PDP) does not depend on the metal temperature".
- p.316: "When closing the shield to prevent the membrane from facing the plasma, we observed no increase in the downstream pressure, that is, no permeation flux through the membrane".

**Замечание для W.** Измерено отношение k_r двух граней (~2), т.е. кинетическая асимметрия; равновесной S через каждую грань не измеряли. Без плазмы (экран закрыт) — потока нет (но и газа на входе почти нет: < 10^-5 Pa).

**Релевантность для W: средняя** — прямое измерение асимметрии k_r^B/k_r^A ≈ 2 у реальной мембраны; ничего о W (S) и о тепловом H2 с двух сторон при равных давлениях.


## 16. Ohyabu et al. 1995 — `hydrogen_metal_iteraction/ref/ohyabu1995.pdf`

**Ссылка:** N. Ohyabu, A. Komori, K. Akaishi, N. Inoue, Y. Kubota, A.I. Livshits, N. Noda, A. Sagara, H. Suzuki, T. Watanabe, O. Motojima, M. Fujiwara, A. Iiyoshi, "Innovative divertor concepts for LHD", *J. Nucl. Mater.* **220–222** (1995) 298–301. SSDI 0022-3115(94)00432-3

**Суть.** Обзор диверторных концепций LHD (helical divertor, local island divertor, сценарии радиационного охлаждения, высокотемпературного дивертора, «SP-operation»). Схемы откачки: углеродные листы (ёмкость 4·10^17 см^-2 при 1 кэВ; 60 м² → 2,4·10^23; хватает на 240 с при 10^21 /с от NBI 20 МВт, 125 кэВ) и — одной фразой — разработка металлической мембранной откачки на основе сверхпроницаемости.

| Величина | Значение | Ед. | Система | T | Стр. |
|---|---|---|---|---|---|
| LHD | R = 3,9 м, B = 4 T | | | | p.298 |
| Ток катушек LID | 200 kA на катушку при B = 3 T | | | | p.299 |
| Ширина перехода closed→open | ~2 мм (LID) vs ~50 мм (HD) | mm | | | p.300 |
| T_edge в HT-режиме | ~5 keV | | | | p.298 |
| Ёмкость C-листов | 4·10^17 cm^-2 (H 1 keV); 2,4·10^23 на 60 m² | | C | <200 °C при разряде; отжиг 700–1000 °C | p.301 |
| Поток пучка NBI | 1·10^21 /s (20 MW, 125 keV) | | | | p.301 |
| Эффективность откачки C-листами | ~40% | | | | p.301 |

**Цитаты.**
- p.301: "For a steady state hydrogen pumping, we are developing a metal membrane pumping, which utilizes superpermeability of the particular metal [12] and are planning a test to demonstrate practicality of the membrane pump."
- О детальном равновесии, растворимости, асимметрии мембран, тепловом H2 — **ничего**.

**Релевантность для W: низкая** — физика плазмы/дивертора; мембраны упомянуты одной фразой.


## 17. Day et al. 2016 — `hydrogen_metal_iteraction/ref/day2016.pdf`

**Ссылка:** Chr. Day, B. Butler, T. Giegerich, P.T. Lang, R. Lawless, B. Meszaros, "Consequences of the technology survey and gap analysis on the EU DEMO R&D programme in tritium, matter injection and vacuum", *Fusion Engineering and Design* (2016), article in press, doi:10.1016/j.fusengdes.2016.03.008 (G Model FUSION-8629).

**Суть.** Системно-инженерный выбор технологий топливного цикла DEMO (Direct Internal Recycling, Pd-Ag пермеатор, TCAP, пеллет-инжекция, ртутные насосы, металлофольговый насос). Металлофольговый насос на сверхпроницаемости получил 70% и включён в программу R&D как средство DIR. Физики мембран — одна фраза.

| Величина | Значение | Ед. | Система | Стр. |
|---|---|---|---|---|
| Оценка: Pd-Ag пермеатор | 82 | % | очистка выхлопа | p.4, Table 2 |
| Оценка: мембранный реактор | 74 | % | | p.4 |
| Оценка: TCAP | 72 | % | удаление протия | p.5 |
| Оценка: криодистилляция | 64 | % | | p.5 |
| Оценка: MgH2+1 mol% Nb2O5 / Zr-Co | 63 / 55 | % | хранение T | p.6 |
| Оценка: сорбционные слои (He-теплоноситель) / CECE | ~70 / 65 | % | | p.6 |
| Оценка: пеллеты | 75 | % | | p.7 |
| Оценка: крионасосы / металлофольговый / паро-диффузионный / жидкостно-кольцевой | 50 / 70 / 70 / 91 | % | | p.7–8 |
| Длительность импульса DEMO | 2 h | | | p.2 |
| Концентрация PEG | возможно >20% | | | p.4 |
| Скорость пеллет | ≥1 km/s; blower guns 300–500 m/s; JET HFPI 30–50 Hz | | | p.7 |
| Криоплант ITER на крионасосы | 25–30% | | | p.7 |

**Цитаты.**
- p.7–8: "Metal foil pumps (70%): Metal foil pumping is a technology that is not (yet) used commercially (R&D topic 14). It exploits the principle of superpermeation [50]. This makes use of the fact that metal foils can show a high permeability for hydrogen atoms at a certain atomic energy, which recombine to form molecules at the back side of the foil. This pump works only for hydrogen isotopes as a pump; it is not applicable for pumping helium and other heavy gases."
- О детальном равновесии, растворимости, асимметрии, тепловом H2 — **ничего**.

**Релевантность для W: низкая** — программный/инженерный документ; сверхпроницаемость упомянута как принцип для энергичных атомов.


## 18. Courteille et al. 1993 — `hydrogen_metal_iteraction/ref/courteille1993.pdf`

**Ссылка:** C. Courteille, J. Bruneteau, F.P.G. Valckx, Z. Sledziewski, M. Bacal, "Primary electron drift in a volume hybrid multicusp H⁻ ion source", *Rev. Sci. Instrum.* **64** (1993) 3265–3269. doi:10.1063/1.1144288

**Суть.** Описание большого гибридного многоостроечного источника H⁻ (тот самый, где Livshits et al. 1994 мерили диссоциацию на Ta-нитях). Наблюдён дрейф первичных электронов (grad B + кривизна) вдоль линий остроёв, направление которого чередуется от нити к нити; «D-нити» в 3,5 раза эффективнее «U-нитей». Расчёт траекторий в приближении ведущего центра. К металлам/водороду в смысле растворения/проницаемости отношения не имеет.

| Величина | Значение | Ед. | Система | Стр. |
|---|---|---|---|---|
| Камера | Ø 44 см, h 35,5 см; 16 Sm-Co магнитов; 3500 G на поверхности магнитов | | | p.3265 |
| Нити | 16 W, Ø 0,5 мм, на окружности Ø 33 см, 5,5 см от стенки | | | p.3265–66 |
| Поле и градиент у нити | ~130 G, 40 G/cm | | | p.3266 |
| Мощность / давление | до 6 kW dc; 1–8 mTorr (оптимум 2–3) | | | p.3265 |
| Дрейф 100-eV электронов | v_D,max = 3,7·10^7 cm/s; расчёт 2,1·10^7 и 1,65·10^7 cm/s | | | p.3266, p.3269 |
| Эффективность D/U нитей | 3,5 | — | | p.3267 |
| Прирост тока H⁻ при сдвиге U-нити на 15 см | 25% | | | p.3268 |
| Плавающий потенциал | −15 … +2 V | | | p.3266 |
| Зеркальное отношение | 22–27 | | | p.3268–69 |

**Цитаты.** Мест о детальном равновесии, растворимости, проницаемости, асимметрии мембран — **нет**.

**Релевантность для W: низкая (нулевая)** — физика плазменного источника.


## 19. Bacal & Wada 2020 — `hydrogen_metal_iteraction/ref/10.1088@1361-6595@ab6881.pdf` (64 стр., прочитано полностью; стр. 54–64 — список литературы)

**Ссылка:** M. Bacal, M. Wada, "Negative ion source operation with deuterium", *Plasma Sources Sci. Technol.* (2020), accepted manuscript, doi:10.1088/1361-6595/ab6881 (PSST-103426.R1).

**Суть.** Исторический обзор источников отрицательных ионов H⁻/D⁻ (двойная перезарядка, объёмные, Cs-добавка, поверхностно-плазменные, конвертерные) и изотопного эффекта (DEA-сечения, ларморовские радиусы, работа выхода, сопутствующие электроны). Металлы как растворители/мембраны водорода не рассматриваются; Livshits et al. 1994 цитируется только как источник данных о плотности атомарного H.

| Величина | Значение | Ед. | Система | Стр. (рукописи) |
|---|---|---|---|---|
| Сродство к электрону H / D | 6082,99 / 6086,2 cm^-1 (0,7542 / 0,7546 eV) | | | p.10 |
| Работа выхода Mo / min с Cs / насыщение | ~4,6 / ~1,5 (0,17 ML Cs) / ~2,1 eV (>0,36 ML) | eV | Mo+Cs | p.17 |
| Макс. эффективность отрицательной ионизации (Rasser) | 4% W(110), 40% Cs/W, 15% Cs | | | p.17 |
| Кинетическая энергия ухода ионов (Hiskes) | 0,77 (Cs), 1,25 (Na) | eV | | p.18 |
| Отрицательная ионизация тепловых H | <1,5·10^-4 ион/атом | | W-печь 1800–2400 K | p.20 |
| Макс. сечение DEA | 28,3·10^-16 cm² (v′=1, J=29) | | H2 | p.15 |
| Порог экзоэргичности DEA | 3,994 eV внутренней энергии; скорость ~10^-8 cm³ s^-1 | | H2 | p.14 |
| Равновесная доля D⁻ в Cs / Rb | 35% (<500 eV) / <30% | | | p.25 |
| Конверсия (Dyachkov) | 4,2% Mg, 5,3% Li, 9,5% Na | | | p.23 |
| T_e(D2)/T_e(H2) | 1,2 | | 3 mTorr | p.33 |
| n(D⁻)max/n(H⁻)max | ≈0,53 | | Camembert III, 80 A | p.33 |
| n⁻/n_e | 0,06 (D2) vs 0,12 (H2); с Cs: 0,14 (D2, 1 mTorr) vs 0,25 (H2) | | | p.37, p.39 |
| Плотность H2/D2 в разряде (CARS) | 71% / 42% от идеального газа (540 K, 55 μbar) | | | p.36 |
| V_p−V_f | 3,532 kT_e (H) / 3,878 kT_e (D) | | | p.35–36 |
| Целевые плотности тока SPIDER | 355 A/m² H⁻, 285 A/m² D⁻ | | | p.7–8 |
| Cs для ITER | 20 mg/h на источник; 40 g на печь | | | p.45 |

**Цитаты.** Единственное место, близкое по смыслу к микрообратимости (в газовой фазе, не на металле), p.10: "Reactions (2) and (3) are the time reversed process of the DEA: e + H2 → H⁻ + H (4) and e + D2 → D⁻ + D (5)". О растворимости, проницаемости, мембранах, асимметрии — **ничего**.

**Релевантность для W: низкая (нулевая)** — газоразрядная физика отрицательных ионов.

---

# Сводка по пакету B относительно W

## Главный вывод
**Ни одна работа пакета не измеряет равновесную константу Сивертса отдельно через грань A и через грань B одной асимметричной мембраны, и ни одна не ставит опыт «тепловой H2 одинакового давления и температуры с обеих сторон асимметричной мембраны».** Прямого экспериментального значения W нет. Все количественные модели пакета закладывают W = 1 по построению: одна объёмная S на слой, а поверхность входит только через кинетику (α, k_r). Ни одного наблюдения, которое противоречило бы W = 1, в пакете нет.

## За W = 1 (прямо или косвенно)
1. **Hatano et al. 2007 (JCP)** — самое прямое из имеющегося.
   - (a) Проницаемость трубчатой Nb-мембраны измерена **в двух противоположных направлениях**. Получены близкие α_H2(T), «rather small difference» авторы относят к погрешности манометров.
   - (b) α_H2, найденная из стационарного прямого потока, совпала с α = (z/S²)·k_r, вычисленной из k_r (по времени запаздывания) **с одной табличной S**. Это проверка детального равновесия на одной поверхности.
   - (c) ТДС с внешней и внутренней поверхности дала близкие потоки.
   - Ограничение: грани делались одинаковыми. Это тест ДР и симметрии, а не тест W при разных гранях.
2. **Musyaev et al. 2008 (V) и 2005 (Nb)** — реальные **сильно асимметричные** мембраны:
   - по коэффициентам прилипания тепловых H2: α_in = 2–5·10^-5, α_out ~ 2·10^-4 (Ta, осаждённый на входную грань);
   - «степень асимметрии» ~10³ (Nb).
   При этом в обеих работах при выключении атомизатора давления возвращаются к исходным: «Reducing the atomizer temperature to its initial value (“switching off” of the atomizer) causes the inlet and outlet pressures to return to their initial levels» (2008, p.524); «On switching off of the atomizer, the inlet and outlet pressures are returned to their starting levels» (2005, p.36). Сжатие (>10³) и перекачка вверх по градиенту есть только при работе атомизатора (надтепловые атомы). В режиме сжатия обратный *молекулярный* поток идёт out→in, т.е. по градиенту. Это качественное свидетельство: асимметрия α не порождает «самонакачки». Тест не количественный: обе камеры откачиваются, равных давлений не создавали.
3. **Nakamura et al. 2000** — асимметрия k_r двух граней ≈ 2 измерена ТДС *при допущении равных концентраций у обеих граней* (равновесие в объёме). Без плазмы потока нет.
4. **Явные модели с одной объёмной S при разных α_in, α_out:**
   - Alimov 2018, ур.(1): √(P_in − j/Zα_in) − √(P_out + j/Zα_out) − 2jΣL_i/(S_iD_i) = 0. При P_in = P_out единственное решение j = 0 при любых α_in ≠ α_out;
   - Hatano 2008: α = k_r K_S²√(2πmkT) с табличной K_S, хотя α и k_r меняются на порядки;
   - Zheludkevich 2008: J = DK_S(√P′ − √P″), «накачка» только за счёт химпотенциала атомарного H (ΔG°_дисс = 406,6 кДж/моль);
   - Suzuki & Yukawa 2021: движущая сила — разность химпотенциалов, концентрации на гранях берутся из равновесных PCT;
   - Nishimura 1991: J = Φ(√P_u − √P_d).
5. **Livshits et al. 1994** — методологический аналог: вероятность диссоциации на любой поверхности в термостате K_s = √2·γ_s·η_e. Поверхность даёт только кинетический множитель, равновесная степень диссоциации η_e задана термодинамикой и от поверхности не зависит.
6. **Alimov 2014** — S(V), измеренная на образцах с Pd-покрытием, близка к литературной S(V), полученной другими авторами на других поверхностях. Это слабый аргумент, что катализатор на поверхности не меняет равновесную S.
7. **Skovoroda 2002, Cherkez 2020, Hanke 2020, Hatano 2006.** Везде, где обсуждается выгодная асимметрия (γ = K_L/K_0, k_f ≪ k_b, «higher inlet than outlet barrier»), она работает **только при надтепловом входном потоке**. Явные формулировки:
   - «permeation efficiency of equilibrium molecular hydrogen … extremely low» (Skovoroda 2002);
   - «Superpermeation can be achieved with the production of suprathermal hydrogen particles» (Hanke 2020).

## Против W = 1 (или требующее осторожности)
Прямых аргументов против **нет**. Есть три наблюдения, которые можно ошибочно принять за «S зависит от поверхности», но и в самих статьях они трактуются иначе:
1. **Nishimura 1991.** K, полученная из проницаемости (Φ/D), в **1/10–4/10** раза меньше равновесной K из PCT на тех же Pd-покрытых V–Ni (разница ~0,2 H/M при 573 K, 10^5 Pa). Авторы объясняют это неравновесием на входной границе под потоком и уменьшением эффективной площади из-за оксида под Pd. Это стационарный кинетический эффект, а не равновесная S через грань.
2. **Skovoroda 2002.** Гистерезис и немонотонность PDP; при ~1000 K сорбция/десорбция необратима («only one-half … desorbed»). Это нестационарные, фазовые и дефектные эффекты, а не стационарный поток в равновесии.
3. **Hanke 2020** («lack of reliability of the permeation fluxes») и **Hatano 2007** (зависимость α от термической истории, медленное, ~1 ч, установление равновесия O на поверхности). Состояние граней дрейфует. Любой опыт по проверке W должен это контролировать: дрейф даёт переходные потоки, не связанные с W.

## Нейтрально (к W отношения не имеют)
Ohyabu 1995, Day 2016, Courteille 1993, Bacal & Wada 2020 — физика дивертора, топливного цикла и источников ионов. Hatano 2006 и Cherkez 2020 содержат полезные α, но только для надтепловых частиц. Suzuki 2021 и Alimov 2014 дают объёмные S, D и ΔH, ΔS — справочные значения S = K1·K2² без информации о гранях.

## Полезные числа для оценки «каналов» в модели
| Параметр | Значение | Источник |
|---|---|---|
| α_H2 на чистом Nb | ≈0,25 (300–1800 K) | Hatano 2007 |
| α_H2 на Nb+O | 4·10^-5 (500 K) … 2·10^-5 (1300 K, max C_O); ∝C_O^-2 | Hatano 2007 |
| α_H2 на Pd-покрытии (не UHV) | 0,6 exp(−50,2 кДж/моль/RT); 7·10^-5 при 400 °C | Alimov 2018 |
| α_H2 на двух гранях одной V-мембраны | 2–5·10^-5 (вход, с Ta) / ~2·10^-4 (выход) | Musyaev 2008 |
| k_r^B/k_r^A одной Nb-мембраны | ≈2 | Nakamura 2000 |
| «Асимметрия» Nb-мембраны | ~10³ | Musyaev 2005 |
| 2E_b на «дырках» O-покрытия / на регулярных O-узлах | 10±2 / >100 кДж/моль | Hatano 2007 |
| K (рекомбинация) Nb | 4·10^-25 exp(−13000/T) m⁴/s | Skovoroda 2002 |
| S_V | 0,24 exp(32646/RT) mol H m^-3 Pa^-½ | Alimov 2018 |
| ΔH растворения H в V | ≈ −32,5 кДж/моль (по графику) | Alimov 2014 |
| ΔS_0,2 для V–X | ≈ −47 J/(mol H·K), не зависит от легирования | Suzuki 2021 |
| α_H (атомы) | 0,1–0,4 (≈0,25), не зависит от T и C_O | Hatano 2007, Cherkez 2020 |
| γ_s (рекомбинация атомов H) | ≈1 на чистой поверхности; 0,1–0,3 на поверхности с монослоем | Livshits 1994 |

## Чего в пакете нет (и что нужно искать в других пакетах)
- Опыта «p_A = p_B, T_A = T_B, асимметричные грани → измерение потока или перепада давлений» (с закрытыми камерами).
- Измерения равновесной S (PCT) на одном образце с разными покрытиями или гранями при прочих равных, с точностью лучше нескольких процентов.
- Ссылок на подповерхностные состояния и порядок кинетики, кроме второго порядка рекомбинации (J = k_r c²) у Hatano 2008 и Skovoroda 2002 и излома Аррениуса у Hatano 2007.
