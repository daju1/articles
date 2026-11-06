# Пакет A: работы группы А.И. Лившица (1990–2023) — что они говорят о W = (S_A/S_B)²

**Вопрос.** Может ли асимметричная металлическая мембрана в тепловом H2 одинакового давления и температуры по обе стороны давать стационарный поток? В кинетической модели газ ⇄ поверхность ⇄ объём ⇄ поверхность ⇄ газ всё сводится к W = (S_A/S_B)² = K1^A (K2^A)² / (K1^B (K2^B)²); детальное равновесие требует W = 1.

**Как читались файлы.** Все 10 PDF прочитаны полностью, до конца, по локальным копиям из tmp/lit/src/hydrogen_metal_iteraction/. Страницы с формулами и таблицами, у которых плохой текстовый слой, отрендерены в картинки и прочитаны глазами: 1990 с.81; 1991 с.13 и 15; 1992 с.162; 1995 с.260; 1997 с.1207; 2002 FST с.883–885; 2004 с.1664; 2023 с.14–15. Числа, снятые с графиков «на глаз», так и помечены. Мои собственные оценки помечены отдельно («моя оценка»).

**Обозначения у авторов.** α (у них также α_0/α_L, α_u/α_d, α_up/α_down) — вероятность диссоциативного поглощения теплового H2 за один удар (коэффициент прилипания). k_r — константа рекомбинационного выхода, j_r = k_r c². S или K — константа Сивертса. z = Z = (2πmkT)^(−1/2). ξ — вероятность поглощения надтепловой частицы. В терминах задачи: α ~ K1·(доля, ушедшая в объём); k_r ~ α Z/S². Во всех статьях пакета принято α = k_r S²/Z с ОДНОЙ и той же S для обеих сторон, а это и есть W = 1.

---

## 1. Livshits, Notkin, Samartsev 1990

**Ссылка:** A.I. Livshits, M.E. Notkin, A.A. Samartsev, "Physico-chemical origin of superpermeability — large-scale effects of surface chemistry on 'hot' hydrogen permeation and absorption in metals", *J. Nucl. Mater.* **170** (1990) 79–94.

**Суть.** Обзорно-экспериментальная статья: сверхпроницаемость (χ ~ 0.1–1 для «горячего» H — атомов/ионов) объясняется тем, что выход растворённого H требует ассоциативной стадии, которая замедляется неметаллическими монослоями (C, S, O) на поверхности. Феноменологическая модель с двумя поверхностями 0 (вход) и L (выход): α_0, α_L — вероятности диссоциативного растворения H2, k_r = αZ/S² с ОДНОЙ и той же объёмной константой Сивертса S для обеих сторон (т.е. W = 1 заложено в модель). Эксперименты: Pd 20 мкм с монослоями C и S, аустенитная сталь под D⁺ с H2S/O2, Nb с Pd-покрытием на выходной стороне (асимметрия α_L ≫ α_0). Явно сформулировано термодинамическое требование: при равных давлениях и изотермии — нет потока; проницаемость для молекулярного равновесного H2 одинакова в обоих направлениях (ур. 10a,10b симметричны по α_0↔α_L).

**Числовые параметры**

| Величина | Значение | Единицы | Система | T | Где |
|---|---|---|---|---|---|
| Z | 1.2×10^19 | H2/(см² с Па) | H2, 20 °C | 293 K | с.80, прил. |
| α_0*, α_L* (предэксп.) | 0.1–1 | — | общее | — | с.80, ур.(2) |
| σ (доля активных центров) | 10⁻³ … 1 | — | по литературе | — | с.81 |
| E_a1 (барьер растворения) | 0 … ~40 | кДж/моль H | по литературе | — | с.81 |
| D* | Fe 0.93; Ni 6.7; Pd 5.2; Nb 0.5; V 0.35; ASS 85 | 10⁻³ см²/с | табл.1 | — | с.81 |
| S* | Fe 2.4; Ni 2.6; Pd 3.3; Nb 0.61; V 1.0; ASS 2.4 | 10^17 ат.H см⁻³ Па⁻¹ᐟ² | табл.1 | — | с.81 |
| E_a3 (уровень растворённого H относ. барьера выхода) | Fe 36; Ni 54; Pd 15; Nb −28; V −28; ASS 67 | кДж/моль H | табл.1 | — | с.81 |
| ΔH растворения | Fe 24; Ni 17; Pd −10; Nb −39; V −33; ASS 8.8 | кДж/моль H | табл.1 | — | с.81 |
| E_a1 для Fe (по Waelbroeck) | 22 | кДж/моль H | Fe, «чистая» пов. | — | с.83 |
| σ для Fe (по Waelbroeck) | 3×10⁻³ | — | Fe | — | с.83 |
| α (чистый Pd) | 0.3–1 = const(T) | — | Pd чистый | — | с.84, 85 |
| α (очищенные пов.) | 0.1–0.6 = const(T), E_a1=0, σ=1 | — | Pd, Nb | — | с.84 |
| α (C-Pd) | 0.15 exp(−16.8 кДж/моль /RT); E_a1 = 8.4 кДж/моль H, σ=1 | — | Pd 20 мкм + монослой C | ≤500 K | с.85 |
| α (S-Pd) | 0.1 exp(−34.4 кДж/моль /RT); E_a1 = 17.2 кДж/моль H, σ=1 | — | Pd + монослой S | >300 K | с.85 |
| ξ (растворение тепловых атомов через C, S) | 0.2 = const(T) | — | C-Pd-C, S-Pd-S | — | с.86 |
| χ (атомы, сверхпроницаемость) | 0.1 = const(T) | — | C-Pd-C, S-Pd-S | — | с.86 |
| j_i для условия (9) | <10^15 (C-Pd-C), ≲10^17 (S-Pd-S) | см⁻² с⁻¹ | Pd 20 мкм | >300 K | с.85 |
| α_0 | 8×10⁻³ exp(−39 кДж/моль /RT) | — | аустенитная сталь (до ионной бомбардировки) | ~570–940 K | с.88, сноска |
| k_r0 | 1.5×10⁻¹⁸ exp(−71.5 кДж/моль /RT) | (ед. см⁴/с) | аустенитная сталь | — | с.88 |
| ξ(E) | 0.2, 0.4, 0.7, 1.0 при E = 10, 10², 10³, 5×10³ эВ | — | ионы | — | с.88 |
| Рост χ при O2 4×10⁻⁴ Па | ~10³ раз, χ до ~0.2 | — | ASS 100 мкм | 940 K | с.88 |
| χ при асимметрии Nb + O2/H2S | ≈1 | — | Nb | 570–1370 K | с.91 |
| Рост потока после очистки выходной стороны | ≈×2 (чуть меньше) | — | Pd 20 мкм | 440 K | с.86, рис.3 |

**Ключевые формулы (с.81–82):** k_r = αZ/S²; для равновесного H2: χ = min{α_0α_L/(α_0+α_L); (SD/L)/(Z√p)} (10a,b); симметричная мембрана χ = α/2 (10c).

**Цитаты**

- p.80: "From balancing the fluxes, it follows that k_r = αZ/S², where S is the equilibrium (Sieverts') constant for dilute solution of hydrogen in the metal contacting gas phase hydrogen"
- p.82: "Permeation in the regime χ ∝ 1/L is often spoken of as 'diffusion controlled', which seems to us somewhat unlucky. Such a term ought to point up that surface processes have no effect upon permeability, but this applies to an equilibrium gas only (eq. (10b)), whereas for a 'hot' gas, permeability indispensably depends on the characteristics of upstream surface (eq. (4a))"
- p.82: "The steady state concentration in a sample of metal plunged into a 'hot' gas essentially depends upon the surface properties too, in contradistinction to equilibrium concentration c_s = S√p = S*√p exp(−ΔH/RT)."
- p.82 (§2.4 Non-symmetric membranes): "Seen from eqs. (4) is that a non-symmetric membrane (α_0 ≠ α_L) may exhibit conspicuous differences in the permeability to 'hot' gas for the two opposite directions of gas flow (i.e. for α_0 and α_L interchanged). This trait is peculiar to non-equilibrium conditions, whereas the permeability to the molecular hydrogen gas is the same in either direction for both the extremes given by eqs. (10a) and (10b). We have to dwell on the latter point, since arguing about it resumes intermittently in the literature [34,35]."
- p.82: "Permeation differences in the different directions are not entirely prohibited for a non-symmetric membrane separating an equilibrium gas at some specific pressure from the vacuum, or from the gas at another pressure, if more complicated mechanisms of hydrogen-metal interaction are involved than those assumed to obtain eqs. (4) and (10) *. But any kinetic formalism being contemplated to account for an asymmetric permeation should meet one requirement stipulated by the principles of thermodynamics: no regular gas flow may exist under isothermic conditions, if hydrogen pressures are equal on the both sides of a membrane."
- p.82, сноска: "To be accurate, the difference is possible even in this case up to about 60% for the permeation regime transitional from eq. (10a) to eq. (10b) [6,36]."
- p.82–83: "This requirement appears to be disregarded by alleging that the deposition of palladium overlayers only on the upstream surface of the originally symmetric membranes of stainless steel [37] or niobium [38] should bring about a sharp rise of their permeability to molecular hydrogen owed simply to the easier hydrogen dissociation and dissolution characteristic of palladium. A similar argument [34] will infer a non-zero gas flow, if applied to the situation of equal pressures at either side of the asymmetric membrane. Again, the possibility of easier dissolution on the upstream side is included in eq. (10a). The quantities α_0 and α_L, pertaining to the two membrane sides, are equipollent there, and a two-fold permeation increase maximally may result from enhancing the rate of dissolution at only one side of a symmetric membrane."
- p.83: "In a real experimental set-up, the temperature of the gas acting on a metal sample may differ from the metal temperature. … A palladium membrane having a very thin copper coating on one surface may, therefore, display an asymmetric permeability to a gas that is 'colder' than the membrane. The permeation will be smaller when the gas acts on the copper-covered surface which repulses the slow molecules falling from the gas. The asymmetry is to disappear on equalizing the gas and the metal temperature."
- p.84: "A membrane, if originally symmetric (α_0 = α_L), remains so, even though only the upstream side is subject to the action of atomic hydrogen."
- p.86: "With the barrier (E_a1)_L vanishing during the cleaning treatment, the time, τ, for hydrogen atoms to stay inside the membrane falls by orders of magnitude, and the steady state permeation increases by approximately a factor of two. This complies with eqs. (4a) and (5) when α_L = α_0 changes to α_L ≫ α_0." (рис. 3: «the increment here is slightly less than two-fold, because an undesigned initial passivation hardly ever is exactly symmetric»)
- p.88, сноска: "…the x values were measured … for thermal molecules. The latter was found to be independent of the driving pressure, p … Consequently: (i) x = α/2 (eq. (10c)) for H2 molecules, as long as α_0 = α_L before the ion bombardment"
- p.91: "Coating the downstream surface by a palladium overlayer after cleaning results in asymmetry beneficial to superpermeability (α_L ≫ α_0) even when the upstream niobium surface is clean, since intrinsically α_Pd > α_Nb, as demonstrated by Pick et al. [69]."

**Чего нет:** прямых измерений проницаемости H2 в двух противоположных направлениях через асимметричную мембрану с количественным сравнением в этой статье нет; нет измерения S через разные грани; нет наблюдения потока при равных давлениях. Подповерхностные состояния явно не вводятся (есть только барьеры E_a1, E_a2, E_a3 на энергетической диаграмме рис.1). Порядок десорбции: 2-й (j_r = k_r c²), с оговоркой, что для ASS возможен j_r ∝ c (с.82).

**Релевантность для W: высокая** — прямо постулирует (через k_r = αZ/S² с единым S и термодинамический аргумент) W = 1 и указывает, что асимметрия проницаемости для H2 возникает только при неравновесии газа (горячий газ, газ холоднее мембраны); при этом W=1 — постулат модели, а не измерение.



---

## 2. Livshits, Notkin, Samartsev, Grigoriadi 1991

**Ссылка:** A.I. Livshits, M.E. Notkin, A.A. Samartsev, I.P. Grigoriadi, "Large-scale effects of H2O and O2 on the absorption and permeation in Nb of energetic hydrogen particles", *J. Nucl. Mater.* **178** (1991) 1–18.

**Суть.** Мембрана Nb 0.1 мм, вход — пучок 20 кэВ H2⁺ (10 кэВ/протон, 1.5×10^15 H⁺/см²с), напуск H2O или O2 на входную сторону; выходная поверхность «натуральная» (O/C). Оксидная пассивация входа резко увеличивает поглощение и проницаемость (до сверхпроницаемости j_L/j_i ≈ 0.75). Извлекаются константы рекомбинации k_r0 (вход, под пучком) и k_rL (выход); используется соотношение k_r = 2Zα/S² с единой объёмной S (т.е. W=1 по построению). Обсуждается максимальная степень сжатия g_max и обратный поток через мембрану 2Zp(1/α_0+1/α_L)⁻¹ — симметричная по α формула для молекулярного H2. Измерений с молекулярным H2 с обеих сторон и в двух направлениях нет.

**Числовые параметры**

| Величина | Значение | Единицы | Система | T | Где |
|---|---|---|---|---|---|
| D (H в Nb) | 5×10⁻⁴ exp(−10450 Дж/моль H /RT) | см²/с | Nb | — | с.9, 14 |
| S (Nb) | 7×10^16 exp(39700 Дж/моль H /RT) (на с.13 напечатано «7×10⁶» — опечатка оригинала; на с.15 — 7×10^16) | ат.H/(см³ Па¹ᐟ²) | Nb | — | с.13, 15 |
| Z | ≈1.1×10^19 | мол. H2/(см² с Па) | H2 | комн. | с.13 |
| α* (чистый Nb) | 0.1; E_0 = 0 | — | Nb чистый | — | с.13 |
| α (чистый Nb) | ~0.1, очень слабая T-зависимость | — | Nb | — | с.7 |
| k_r (чистый Nb, расчёт по k_r=2Zα/S²) | 3×10⁻²⁴; 1×10⁻²²; 6×10⁻²¹ | см⁴ (ат.H)⁻¹ с⁻¹ | Nb чистый | 500; 623; 838 K | с.13 |
| ΔH_c хемосорбции (Pick) | −56.4 | кДж/моль H | Nb, монослой H | — | с.11 |
| ΔH_c (собств. неопубл.) | −75 | кДж/моль H | Nb лента 10 мкм, очищ. 2400 K | — | с.11 |
| j_max десорбции через монослой | 10^15 exp(−2ΔH_c/RT); оценка 3×10^13 | ат.H/(см² с) | Nb | 400 K | с.11 |
| Уровень E_3 растворённого H ниже H в H2 | ~30 | кДж/моль H | Nb (V, Ta) | — | с.9 |
| Сверхпроницаемость j_0∞/j_i ; j_L∞/j_i | 0.25 ; ~0.75 | — | Nb + H2O 1.3×10⁻³ Па | 838 K | с.8 |
| Поток сегрегации примесей на поверхность | ~2×10^13 / ~1×10^12 | ат./(см² с) | Nb | 810 / 670 K | с.13, сноска 11 |
| Коэф. прилипания O2, H2O при оксидировании под пучком | ~1 | — | Nb | — | с.10 |
| Re-emission ratio Ψ_∞ | max{L²ξj_i k_r0/D²; k_r0/k_rL = α_0/α_L} | — | — | — | с.15, ур.(7) |
| g_max | βξ/α_0 = (βξ/α_0*) exp(2E_0/RT) | — | — | — | с.15 |

**Цитаты**

- p.1: "The quantity 2E_0 represents the activation energy for dissociative binding of a thermal H2 molecule impinging on the metal surface. The surface barrier also hinders the reverse process of the recombinative evolution of dissolved hydrogen, as it is clear considering the balance of fluxes at thermal equilibrium."
- p.3, сноска 1: "For two-layer membranes, χ > 1/2 is possible even with both the inlet and outlet surfaces being clean, on the condition that the outlet surface is characterized by a higher sticking coefficient, α, for hydrogen molecules. So, permeation with χ > 1/2 was observed for a Nb-Pd membrane having the both surfaces clean [8], in accordance with that α_Pd > α_Nb is known [36]." (речь о «горячем» водороде)
- p.13, сноска 10: "The formula [2,5,8,21] is obtained by combining Sieverts' law, C_s = S√p, and the balance equation, 2Zpα = k_r C_s², where C_s is the equilibrium concentration of hydrogen dissolved in metal corresponding to the hydrogen pressure, p … Note that eq. (1) is valid for the metals occluding hydrogen both exothermically and endothermically."
- p.13: "We would like to warn against a simplistic approach to interpretations in the terms of eq. (1) of the temperature dependence of k_r … the values of both the pre-exponent α* and the activation barrier E_0 may depend on metal temperature, because the state of the surface is determined by a dynamic balance between such processes as ion sputtering, oxidation …"
- p.12: "Unlike k_r0, the rate constant k_rL is not changing after the ion beam initiation, and it is not influenced by the presence of H2O or O2 on the upstream side."
- p.15, сноска 13: "More accurately, the backstream flux equals 2Zp(1/α_0 + 1/α_L)⁻¹ for the C_0∞ = C_L∞ regime (e.g. ref. [8]). Considering that α_0 < α_L (needed for superpermeability) one gets the estimate in the text."
- p.15: "Thus Ψ_∞ can never be smaller than k_r0/k_rL." (k_r0/k_rL = α_0/α_L — отношение констант рекомбинации двух сторон приравнивается отношению α, т.е. S одинаково для обеих сторон.)

**Чего нет:** нет измерений с молекулярным H2 при равных давлениях; нет измерений S через разные поверхности; подповерхностные состояния упоминаются лишь качественно («oxidation … affects … subsurface layers», с.9; ловушки радиационных дефектов, с.12).

**Релевантность для W: средняя** — использует и явно выводит k_r = 2Zα/S² из «балансa потоков при равновесии» с единой S, т.е. принимает W=1, и даёт формулу обратного молекулярного потока, симметричную по α_0, α_L; прямых проверок нет.



---

## 3. Livshits, Notkin, Samartsev, Busnyuk, Doroshin, Pistunovich 1992

**Ссылка:** A.I. Livshits, M.E. Notkin, A.A. Samartsev, A.O. Busnyuk, A.Yu. Doroshin, V.I. Pistunovich, "Superpermeability to fast and thermal hydrogen particles: applications to the pumping and recycling of hydrogen isotopes", *J. Nucl. Mater.* **196–198** (1992) 159–163. (Файл: livishits1992.pdf)

**Суть.** Короткий обзор/прикладная работа: сверхпроницаемость Nb, Pd, Fe для ионов и тепловых атомов; предложение использовать мембраны в откачных трактах ИТЭР. Модельный эксперимент: цилиндрическая мембрана Nb 100 мкм (∅10 см, 20 см) + атомизаторы (ленты Nb, Ta, W при ~1900 K); откачка ~1.5 л/(с·см²), сжатие ~10⁴. Механизм — неметаллический барьер на входе; поток вызван неравновесным (атомарным/ионным) водородом, получаемым от внешнего источника энергии (горячий атомизатор). Молекулярный H2 при равных условиях с обеих сторон не исследуется.

**Числовые параметры**

| Величина | Значение | Единицы | Система | T | Где |
|---|---|---|---|---|---|
| Энергия ионов | 19.2 кэВ H3⁺ (6.4 кэВ/протон); 200 эВ H3⁺ (67 эВ/протон) | — | Nb 100 мкм | 893 K; 838 K | с.159–160, рис.1–2 |
| Вероятность поглощения H-атома за удар | ≈0.2, не зависит от T и от химии поверхности | — | Pd + монослой S или C | от комн. T | с.160, рис.3 |
| Отношение проницаемостей атомы/молекулы | до ~10^10 (текст плохо распознан: «101~»; по контексту — 10^10) | — | Fe 20 мкм, неполностью очищ. | — | с.160, рис.4 |
| Скорость откачки | 800–1000 л/с при отношении площадей атомизатор/мембрана 0.15; ≈1.5 л/(с·см²) мембраны; ≈10 л/(с·см²) атомизатора | — | Nb 100 мкм, атомизаторы Nb/Ta/W | — | с.162 |
| Степень сжатия | ≈10⁴ (без заметного падения скорости откачки) | — | Nb | — | с.162, рис.7 |
| Вероятность диссоциации H2 на атомизаторе (с рис. 8, на глаз) | выходит на ≈0.3 (W), ≈0.2 (Nb, Ta) при T ≳ 1700 K; ~0 при 1200 K | — | W, Nb, Ta при p_H2 = 1.3×10⁻⁴ Па | 1200–2100 K | с.162, рис.8 |
| Проектные: удельная скорость | 12–14 | м³ с⁻¹ м⁻² | ИТЭР (оценка) | — | с.162 |
| Проектные: мощность на единицу скорости | 1.7×10⁻³ | МВт/(м³ с⁻¹) | — | — | с.162 |
| Проектные: площадь / длина / мощность / сжатие / запас T | 160 м²; 1.6 м; 5 МВт; 10⁴–10⁵; 5–10 г при сжатии 10³ | — | — | — | с.162 |

**Цитаты**

- p.159: "A potential barrier may be built up at a metal surface, originating from a nonmetallic impurity layer, that effectively impedes the dissociative uptake of molecular hydrogen by the metal. Such a barrier suppresses still further the metal permeability to hydrogen, small as it is by itself. Suprathermal particles (including thermal H atoms) are able to freely pass the barrier. The nonmetallic layer impedes, however, the back release of absorbed hydrogen particles through inhibiting their thermal recombination at the surface …"
- p.159: "With only thermal H2 molecules being released at the outlet side, an automatic compression becomes possible of the permeating hydrogen due to a very small permeability of the membrane to molecular hydrogen gas."
- p.160: "The probability of absorption of a H atom, per single impact on the membrane surface (≈0.2), does not depend on the metal temperature … and is insensitive to the chemical composition of the surface. … The membrane permeability to H atoms exceeds its permeability to molecular hydrogen gas by several orders of magnitude."
- p.162: "Thermalized H2 molecules will be atomized at the incandescent surfaces specially introduced." (рабочее тело — неравновесный атомарный H от источника ~1900 K)

**Чего нет:** нет обсуждения детального равновесия, нет S/k_r, нет опытов с H2 с двух сторон.

**Релевантность для W: низкая** — сжатие и поток здесь явно питаются внешним неравновесием (ионы, атомы от горячего атомизатора), а не асимметрией мембраны при тепловом H2; косвенно демонстрирует, что без «горячего» водорода сжатия не получают.



---

## 4. Livshits, Notkin, Pistunovich, Bacal, Busnyuk 1995

**Ссылка:** A.I. Livshits, M.E. Notkin, V.I. Pistunovich, M. Bacal, A.O. Busnyuk, "Superpermeability: critical points for applications in fusion", *J. Nucl. Mater.* **220–222** (1995) 259–263.

**Суть.** Оценка верхнего предела плотности потока при сверхпроницаемости (теоретически ~10^19 см⁻² с⁻¹ для металлов V группы), экспериментально достигнуто 3×10^16 (мембрана Nb 600 см², атомизаторы Nb/Ta) и 4×10^17 см⁻² с⁻¹ (имитационные опыты по поглощению на ленте Ta). 3000-часовой ресурсный тест в присутствии CH4, CO, O2, H2O — деградации нет. Энергетическая диаграмма (рис.1) с барьером рекомбинации E_r (вход, неметаллическая плёнка) и E_r^out (выход), уровнем растворённого H E_b, ямой адсорбции E_ad. Источник потока — атомарный H с атомизатора (неравновесие). Тепловой H2 с двух сторон не рассматривается.

**Числовые параметры**

| Величина | Значение | Единицы | Система | T | Где |
|---|---|---|---|---|---|
| Порог «надтепловой» энергии частиц | ≳1 | эВ | общее | — | с.259 |
| j_max (2-й порядок рекомбинации) | ∝ (1/L²) exp[2(E_r − E_b)/RT] | — | — | — | с.260, ур.(1) |
| E_b (уровень растворённого H относ. ½H2) | ≈ −0.3 (Nb, Ta); +0.15 (Pd) | эВ | — | — | с.260 |
| j_max (монослой, Эли–Ридил) | ≈10^15×10^13 exp[−2(E_r−E_ad)/RT]; ≥10^20 при T≥800 K, если E_r−E_ad ≤ 0.6 эВ | см⁻² с⁻¹ | — | — | с.260, ур.(2) |
| Энергия H+Me относ. ½H2+Me | 2.2 | эВ | диаграмма | — | с.260, рис.1 |
| D (V группа) | ≥10⁻⁴ | см²/с | V, Nb, Ta | ≥800 K | с.260 |
| k_r^out (выходная поверхность) | «≥10^24» напечатано — очевидно ≥10⁻²⁴ | см⁴/с | Nb | — | с.260 |
| C_max (допустимая) | 5 | ат.% | Nb | — | с.260 |
| j_max теоретический | ~10^19 | см⁻² с⁻¹ | Nb, L=0.3 мм | — | с.260 |
| Ранее исследованные потоки | ≤5×10^14 | см⁻² с⁻¹ | — | — | с.260 |
| Достигнутый поток проникновения | ≈3×10^16 | см⁻² с⁻¹ | Nb 0.1 мм, 600 см² | 850–900 K | с.260–261, рис.2 |
| Поток поглощения (имитация) | до 4×10^17; вероятность поглощения 0.1–0.3, не зависит от T | см⁻² с⁻¹ | Ta лента 0.05 мм, 400 см² | 850–900 K | с.261 |
| Удельная скорость откачки | 1.5–2 л/с на см²; 16 м³/с на м² | — | Nb мембрана + атомизатор | — | с.259, 262 |
| Подавление разложения CH4 монослоем O | 3–5 порядков | — | Ta, Nb | — | с.261, рис.3 |
| Давления в ресурсном тесте | H2O ~5×10⁻⁶; CO 3×10⁻⁵ (фон); CO ≈6×10⁻⁴; CH4 ≈2×10⁻³ (напуск) | Па | Nb мембрана; атомизаторы 1800–2100 K | 800–950 K | с.262, рис.4 |
| ИТЭР-оценки | разделение 90–99.9 %; мощность 0.05–5 МВт; запас T ≤1.5 г; 2000 м³/с; поток (0.2–3)×10^16 ат/см²с; сжатие 10³–10⁴; давление у мембраны ~5×10⁻² Па | — | — | — | с.262–263 |

**Цитаты**

- p.260: "On having entered the metal, the particles 'cool down', and atoms equilibrated in the state of solid solution are now effectively prevented by the surface barrier from associative reemission. Due to such a 'surface blocking', the dissolved atoms may reach the opposite side of the membrane due to diffusion with a higher probability than be reemitted back, and thus their major part may pass through the membrane, if the outlet barrier is lower than the inlet one (as in Fig. 1)."
- p.260: "As the dissolved atoms may overcome the recombination barrier, E_r, only through their association in pairs (a second-order process), the blocking efficiency falls with the flux density rise until superpermeability disappears."
- p.261: "It should be remembered however that the surface covered with monoatomic nonmetal layers is, typically, a natural state of thermodynamic equilibrium for a metal surface in vacuum: viz., there are mechanisms which 'automatically' maintain a more or less monoatomic coverage through both removing the excessive atoms via their dissolution in bulk metal and replenishing the missing ones via the surface segregation and reactions with ambient gases."
- p.261: "Many monoatomic nonmetal films are absolutely stable under the action of thermal atomic hydrogen … That has proved being true for the systems O-Ni, O-Fe, S-Pd, C-Pd, as well as for O-Nb and O-Ta."

**Чего нет:** детальное равновесие, S через разные грани, опыты с H2 с обеих сторон — отсутствуют.

**Релевантность для W: низкая** — диаграмма рис.1 показывает разные барьеры на входе и выходе при едином уровне растворённого H (E_b, ΔH) в объёме (т.е. асимметрия кинетическая, а не термодинамическая), но вопрос W напрямую не обсуждается.



---

## 5. Livshits, Ohyabu, Notkin, Alimov, Suzuki, Samartsev, Solovyov, Grigoriadi, Glebovsky, Busnyuk, Doroshin, Komatsu 1997

**Ссылка:** A. Livshits, N. Ohyabu, M. Notkin, V. Alimov, H. Suzuki, A. Samartsev, M. Solovyov, I. Grigoriadi, A. Glebovsky, A. Busnyuk, A. Doroshin, K. Komatsu, "Applications of superpermeable membranes in fusion: the flux density problem and experimental progress", *J. Nucl. Mater.* **241–243** (1997) 1203–1209.

**Суть.** Цилиндрическая мембрана Nb 100 мкм (~1000 см²) + раскалённый атомизатор Ta (проволока 1.5 мм, 125 см² из 200 см²) в немаслянотщательном (не-UHV) вакууме: скорость откачки ~10³ л/с (~1 л/с·см²), рекордная плотность потока 3×10^17 H/см²с при 3×10⁻² Торр, сжатие до 10³. Ключевое экспериментальное наблюдение: при холодном атомизаторе (только молекулярный H2 с входной стороны) проникновения не замечено; поток и сжатие появляются только при включении атомизатора (неравновесный H). Атомизация описана формулой с использованием константы равновесия H+H⇄H2 при T атомизатора (детальный баланс на горячей поверхности). 2000-часовая кампания без деградации.

**Числовые параметры**

| Величина | Значение | Единицы | Система | T | Где |
|---|---|---|---|---|---|
| Порог «надтепловой» энергии | ~1 | эВ | — | — | с.1203 |
| Мембрана | ∅12.5 см × 25 см, ~1000 см², 100 мкм | — | Nb | — | с.1204 |
| Атомизатор | Ta проволока 1.5 мм, 200 см² (использовано 125 см²) | — | Ta | T_a = 1910–2170 K | с.1204–1206 |
| Остаточное давление | 5×10⁻⁷ | Торр | H2O, органика (масляная система) | — | с.1205 |
| S_m | ~10³ л/с; ~1 л/с на см² мембраны; ~8 л/с на см² атомизатора | — | — | — | с.1203, 1205 |
| J_m max | 3×10^17 | ат.H/(см² с) | Nb | T_m = 1020 K | с.1205–1206 |
| p max | 3×10⁻² | Торр | — | — | с.1206 |
| Сжатие (p_down/p_up) | до 10³ | — | — | — | с.1205 |
| Пары (T_a, T_m) | (2170 K, 1020 K); (2080 K, 960 K); (?, 910 K) | K | — | — | рис.2 |
| α_H2 = γ (принято) | 0.3 | — | атомизатор | — | с.1208 |
| Атомная эффективность δ | 0.6–0.7 | — | система | — | с.1208 |
| β (вероятность атомизации) | α_H2 при p≪p*; 19.5·γT_a^0.86/√p·exp(−217 кДж/моль/RT_a)·(293/T_a)^¼ при p≫p* | — | — | — | с.1207, ур.(6),(7) |
| Энергия активации β (высокие p) | 217 (=½ энергии связи H2); измерено 200–250 | кДж/моль | — | — | с.1208 |
| τ_p, τ_D | 2 с; 0.15 с | с | Nb, J_m = 2×10^17 | 1020 K | с.1208 |
| D | 1.5×10⁻⁴ | см²/с | H в Nb | 1020 K | с.1208 |
| S_m (л/с) | J_m (мол/с)/p (Торр)/3.3×10^19 | — | — | — | ур.(1) |
| Кампания | ~2000 ч; T_m 500–700 °C; p 10⁻⁵–10⁻² Торр; пик. поток атомов ~10^18 см⁻²с⁻¹; флюенс ~10^22 см⁻²; H2O ~10⁶ Л; масло 10⁵–10⁶ Л | — | — | — | с.1209 |

**Цитаты**

- p.1203–1204: "The barrier for the dissociative absorption of hydrogen molecules at the membrane inlet surface is responsible for the phenomenon [1,2]. … The thin non-metallic overlayers exhibit a highly selective permeability to suprathermal hydrogen, thus inhibiting both the back release of the incoming particles that 'cool down' in the metal lattice and the back permeation of molecular hydrogen evolved at the downstream side. A membrane with symmetric surfaces may let through a half of the flux absorbed, and virtually the whole incoming flux may pass, if the evolution is easier at the downstream surface [1,2,6]."
- p.1204: "It is important that a metal surface covered with nearly monatomic non-metal layers (a 'real' surface) is thermodynamically stable and it reveals a remarkable capability for self-maintenance under the action of a wide variety of agents [10]."
- p.1205: "No permeation could be noticed when hydrogen was admitted at the upstream side with a cold atomizer. Switching on of the atomizer results in an upstream pressure drop and in a steep downstream pressure rise … it may be orders of magnitude higher (e.g. 10³ times in our experiments) than the upstream pressure (i.e. the membrane compresses the permeating gas). Thus a steady-state membrane pumping due to the atomic-driven permeation is observed."
- p.1207, сноска 2: "An equilibrium constant for the H + H ↔ H2 reaction referring to a pressure at the reaction temperature was employed [18] to derive Eq. (6)."
- p.1208: "The condition τ_p ≫ τ_D means that the absorbed atoms cross the membrane many (τ_p/τ_D) times before desorption and it is actually this fact that gives rise to superpermeation (if an additional requirement is also fulfilled that the gas evolution through the downstream surface is at least no more difficult than the evolution back to the upstream side)."

**Чего нет:** численных α, S, k_r для мембраны нет; детальное равновесие для мембраны не обсуждается; H2 с обеих сторон при равных p — нет.

**Релевантность для W: средняя** — содержит прямое экспериментальное наблюдение: асимметричная (сверхпроницаемая, вход пассивирован) мембрана при чисто молекулярном H2 (холодный атомизатор; газ, правда, комнатной температуры, а мембрана — горячая) не даёт заметного потока; сжатие возникает только при подводе энергии через атомизатор. Количественной верхней оценки W это не даёт (порог обнаружения не указан).



---

## 6. Livshits, Sube, Solovyev, Notkin, Bacal 1998

**Ссылка:** A.I. Livshits, F. Sube, M.N. Solovyev, M.E. Notkin, M. Bacal, "Plasma driven superpermeation of hydrogen through group Va metals", *J. Appl. Phys.* **84** (1998) 2558–2564. doi:10.1063/1.368418

**Суть.** Трубчатые мембраны Nb и V (∅1 см, 18 см, стенка 0.1 мм) с оксикарбидным монослоем в водородной плазме (n ≈ 5×10^10 см⁻³, T_e ≈ 0.6 эВ); энергия ионов меняется смещением мембраны (от единиц эВ до 250 эВ). Последовательно наблюдаются проницаемость при тепловых молекулах (очень мала), при атомах с раскалённых Ta нитей (рост > на порядок) и при плазме (ещё >10×). Измерены максимальное давление сжатия P_cm (баланс прямого «плазменного» потока и обратного молекулярного потока) и отсюда α_0 (вход под плазмой) и асимметрия α_L/α_0 ≈ 10; α_L = 1.2×10⁻⁴. Явно записаны формулы α = k_r S²/Z с ОДНОЙ S для обеих сторон и обратный молекулярный поток Z P_c α_0α_L/(α_0+α_L), симметричный по α_0, α_L. Стабильная сверхпроницаемость до E_k ≈ 50 эВ (порог физического распыления O/C-монослоя).

**Числовые параметры**

| Величина | Значение | Единицы | Система | T | Где |
|---|---|---|---|---|---|
| Энергия свободного атома H относ. H в H2 | 2.2 | эВ | — | — | с.2558, рис.1 |
| Порог надтепловой энергии | ~1 | эВ | — | — | с.2558 |
| α* (предэксп. в α = α*σ exp(−2E_r/RT)) | 0.1–1 | — | общее | — | с.2559, ур.(2) |
| Связь α ↔ k_r | α_0 = k_r0 S²/Z; α_L = k_rL S²/Z (одна S) | — | — | — | с.2559, ур.(3) |
| Z | (2πmkT)^(−1/2) | — | — | комн. T газа | с.2559 |
| ΔH, E_b | ΔH < 0, E_b < 0 для V-группы | — | V, Nb, Ta | — | с.2559, рис.1 |
| ΔH растворения | 0.4 (по модулю; «dissolution enthalpy») | эВ | H–Nb | — | с.2562 |
| α (чистая металлическая пов., нижняя типичная граница) | 0.1 | — | — | — | с.2560 |
| Плазма | n ≈ 5×10^10 см⁻³; T_e ≈ 0.6 эВ; ток 25–30 А, 60 В; p = 1–3 мТорр; атомов — несколько %; T атомов — десятые эВ | — | — | — | с.2560 |
| Примеси | ~10⁻⁶ Торр (CO), H2O 2×10⁻⁷ Торр | Торр | — | — | с.2560 |
| Мембраны | Nb, V; ∅1 см, 18 см, 0.01 см | — | — | 550–1000 °C | с.2561 |
| p_up | 1.6 | мТорр | V | 800 °C | с.2561 |
| ΔP_c от молекул | ≈0.07 | мТорр | V | 800 °C | с.2561 |
| Рост потока: нити Ta 2000–2100 °C | > на порядок; затем плазма (25 А, 62 В): ещё >×10 | — | V | 800 °C | с.2561 |
| j_p | ≈3×10^18 × P_c [Торр] ат/см²с (при откачке ~2 л/с) | — | — | — | рис.3 |
| j_p / j_i(атомы) | 0.1–0.3 | — | V/Nb | — | с.2561 |
| Сжатие q | ≈10³ | — | V | 800 °C | с.2561 |
| j_p max | >10^17; 2×10^17 | см⁻² с⁻¹ | Nb, V | — | с.2562, 2564 |
| Теоретический предел j_p | ≈10^19 | см⁻² с⁻¹ | — | — | с.2562 |
| Плавающий потенциал | ≈ −3 В отн. плазмы | В | — | — | с.2562 |
| Энергия рекомбинации H2⁺ | [15.4 + (0.4×2)] эВ/ион | эВ | Nb | — | с.2562 |
| Доля ионов в падающем потоке | ≤10 % от потока атомов | — | — | — | с.2562 |
| Стабильность сверхпроницаемости | до E_k = 60–70 эВ (j_p); α_0 неизменна до ~50 эВ | эВ | Nb | — | с.2563–2564 |
| α_L (выход) | 1.2×10⁻⁴ (не зависит от смещения) | — | Nb, натуральная пов. | — | с.2563, рис.7 |
| α_L/α_0 | до ≈10 | — | Nb под плазмой | — | с.2564, рис.7 |
| Плотность ионного потока | ≈10^16 | протон/(см² с) | — | — | с.2564 |
| Длительность кампании | ~300 | ч | Nb, V | 550–1000 °C | с.2562 |
| Ионные пучки (ссылки): деградация при | E > 70 | эВ/протон | Nb | — | с.2560 |

**Ключевые формулы:** j_p = Σ ξ_k j_ik · α_L/(α_0+α_L) (1); обратный молекулярный поток = Z P_c α_0α_L/(α_0+α_L); P_cm = (1/Z)·Σξ_k j_ik · (1/α_0) (6).

**Цитаты**

- p.2558: "The downstream barrier E_rL − E_b also reflects the absorbed atoms, and thus the atoms 'oscillate' by means of random-walk diffusion between the two boundaries. As a result, the absorbed atoms 'forget' their starting position and are being desorbed on either of the two sides of the membrane depending only on the relative heights of the corresponding surface barriers."
- p.2558: "(b) The permeating hydrogen can be automatically compressed by orders of magnitude, i.e., it can pass through the membrane from the lower pressure side to the side where the hydrogen pressure is higher by several orders. This is related to the fact that the absorbed atoms are thermally desorbed from the metal at its moderate temperature (<1000 °C) only when associated into molecules and that the superpermeable membrane is virtually impermeable to thermal molecules."
- p.2559: "The sticking coefficients are related to the often employed (e.g., Refs. 10 and 11) rate constants of recombinative release, k_r, via: α_0 = k_r0 S²/Z; α_L = k_rL S²/Z, (3) where S is the hydrogen solubility constant"
- p.2559: "It should be noted that, despite these arguments in favor of the possible stability of the monolayer film, the asymmetry condition (5) appears to be very severe, since energetic hydrogen particles are incident on the upstream side only, whereas the downstream side is not subject to anything similar."
- p.2559: "…a metal surface covered with a nonmetal monolayer may be a thermodynamically stable state of the surface under vacuum conditions"
- p.2560: "In particular, an initially symmetric membrane (α_0 = α_L) always remains symmetric under the action of thermal atoms even at a temperature as low as 20 °C"
- p.2561: "In spite of the high membrane temperature, T, (800 °C), the permeation due to hydrogen molecules is very small (the corresponding pressure increment, P_c, is about 0.07 mTorr and is not visible on the scale of Fig. 3), as it usually occurs when a metal surface is passivated with a nonmetallic film."
- p.2561: "Then the downstream pumping is cut off (by closing a valve—see Fig. 2). The downstream pressure rises to reach a maximum level, P_cm, corresponding to the balance of the back permeation flux (driven by molecular hydrogen) and the direct plasma driven flux."
- p.2563: "The back permeation flux caused by the thermal molecules accumulated in the downstream chamber is Z P_c α_0α_L/(α_0+α_L), where Z is the kinetic theory coefficient at room temperature, and P_c is the downstream pressure measured with a gauge at room temperature."
- p.2564: "One can see on Fig. 7 that the asymmetry factor, α_L/α_0, at the plasma driven superpermeation reaches ≈10. The existence of this favorable asymmetry is revealed quite clearly by comparing the permeation and compression graphs (Fig. 6) and taking into account Eqs. (1) and (6); this was verified with a number of other independent methods."

**Замечание к вопросу W.** Модель предполагает, что молекулярный поток при p_up (1.6 мТорр, газ комнатной T, мембрана 800 °C) — ненулевой (ΔP_c ≈ 0.07 мТорр); это прямой поток при перепаде давления, а не поток при равных давлениях. При перекрытии откачки (статическое сжатие) баланс достигается при P_cm ≫ p_up только при включённой плазме. Важно: газ в опытах комнатной температуры, а мембрана горячая — условия неизотермичны (авторы сами указывают Z «at room temperature»).

**Чего нет:** S отдельно по граням не измерялась; W не проверялся; опыта при равных давлениях без плазмы нет.

**Релевантность для W: средняя** — дана количественная асимметрия α_L/α_0 ≈ 10 для реальной мембраны и явная формула детального баланса α = k_r S²/Z с единой S; обратный молекулярный поток симметричен по α, так что при этой модели асимметрия по α не даёт потока при равных p. Но это модельное предположение, а не независимая проверка W.



---

## 7. Livshits, Hatano, Watanabe 2002 (Fusion Sci. Technol.)

**Ссылка:** A.I. Livshits, Y. Hatano, K. Watanabe, "Superpermeability in fusion technology: tritium accumulation and compression", *Fusion Sci. Technol.* **41** (3P2) (2002) 882–886. doi:10.13182/FST02-A22711 (файл livshits2002.pdf; формулы прочитаны по отрендеренным страницам, т.к. текстовый слой формул испорчен).

**Суть.** Теоретическая статья о сжатии и накоплении трития в сверхпроницаемой мембране Nb. Ключевой вывод: максимальное давление сжатия p_d,max = j_iξ/(zα_u) определяется ТОЛЬКО потоком надтепловых частиц и коэффициентом прилипания H2 на входной стороне и не зависит от выходной поверхности, материала, толщины. Вывод опирается на (1) условие, что «газ на выходной стороне в равновесии с растворённым водородом» p_d,max = c_b²/S², и (2) «принцип детального равновесия» k_ru = zα_u/S² — с одной и той же S для обеих сторон. Стационарная концентрация c = [(j_iξ + p_d z α_d)/(α_u+α_d) · S²/z]^½. Оптимальная асимметрия α_d0/α_u0 ≈ 10.

**Числовые параметры**

| Величина | Значение | Единицы | Система | T | Где |
|---|---|---|---|---|---|
| Порог надтепловой энергии | ~1 | эВ | — | — | с.882 |
| Макс. плотность потока SPM | ~10^19 | см⁻² с⁻¹ | V-группа | — | с.882–883 |
| Диапазон strong reflection | до ~10^19 см⁻²с⁻¹; до ~1000 °C; толщины до 1 см | — | V, Nb, Ta | — | с.883 |
| Уровни на диаграмме рис.1 | 29 кДж/моль (подписано у E_b относ. ½H2, под барьером входа) и 39.7 кДж/моль (|ΔH|) | кДж/моль | H–Nb | — | с.883, рис.1 |
| E_b | < 0 | — | V, Nb, Ta | — | с.883 |
| α_u (Nb, легированный 1.5 ат.% O) | 5.5×10⁻⁴ exp(−41.8 кДж/моль /RT) (= 5.5×10⁻⁴ exp(−5000/T)) | — | Nb+O, вход | — | с.884, рис.2; с.885 |
| 2E_ru | 41.8 (эксп.); 62.7 (×1.5, гипотет.) | кДж/моль | Nb+O | — | с.885, рис.3 |
| j_iξ | 1×10^16 (рис.2); 2×10^16 (рис.3) | (D/T или H)/(см² с) | ИТЭР | — | с.884–885 |
| p_d,max | от ~100 Торр (≈400 °C) до ~1 Торр (≈850 °C) (с рис.2, на глаз) | Торр | Nb+O | 400–850 °C | рис.2 |
| ε (рабочее p_d/p_d,max) | 0.1 | — | — | — | с.884 |
| α_d0/α_u0 | 10 (принято «разумной» асимметрией) | — | — | — | с.885 |
| S (H в Nb) | S_0 = 0.6×10^18 H/(см³ Торр^½); S = 0.6×10^18 exp(39.7 кДж/моль /RT) | H/(см³ Торр^½) | Nb | — | с.885 |
| ΔH | −39.7 | кДж/моль | H–Nb | — | с.885 |
| f_m (удельная скорость откачки) | ≈1.5 | л/(см² с) | Nb + атомизатор | — | с.885 |
| ИТЭР-оценка | J ≈ 10^22 (D/T)/с; F_He ≈ 10^5 л/с; A = 10; L = 0.01 см; p_d = 0.1 Торр; F_m = 10^6 л/с; G — несколько г при T ≥ 600 °C | — | — | — | с.885–886, рис.4 |
| T мембраны в опытах с атомизатором | 700–750 °C; с экраном ~1000 °C | °C | Nb | — | с.886 |
| Изотопный эффект растворимости H/D в V | «practically negligible» | — | V | — | с.885 |

**Ключевые формулы (с.883–885):** p_d,max = c_b²/S² (1); j_iξ = k_ru c_b² (2); k_ru = zα_u/S² (3, «detailed balance principle»); p_d,max = j_iξ/(zα_u) (4); α_u = α_u0 exp(−2E_ru/RT) (5); j_iξ + p_d zα_d = c²(k_ru + k_rd) (баланс); c = [(j_iξ + p_d zα_d)/(α_u+α_d) · S²/z]^½ (7); S = S_0 exp(|ΔH|/RT) (10); T_max(p_d) = (2E_ru/R)/ln[p_dα_u0 z/(j_iξε)] (12).

**Цитаты**

- p.883: "If there is no pumping at the downstream side, downstream pressure, p_d, increases up to a value of p_d,max corresponding to a flux balance between the permeation of molecules from the downstream side and that of suprathermal hydrogen from the upstream side. At the balance, gas at the downstream side is in equilibrium with absorbed hydrogen: p_d,max = c_b²/S². (1)"
- p.883: "With using the detailed balance principle, k_ru may be expressed through the probability of molecule absorption (sticking coefficient) α_u as [5] k_ru = zα_u/S², (3)"
- p.884: "One should note that the maximum compression pressure p_d,max depends only on the implanted flux density and properties of upstream surface, i.e. on α_u (without any dependence on membrane material and thickness, state of downstream surface, regime of permeation, etc)."
- p.884: "However, the net permeation flux is zero at p_d = p_d,max."
- p.884: "The condition of preferable release of absorbed hydrogen through the downstream side is α_d ≫ α_u. This favorable membrane asymmetry can be realized practically, e.g. by formation of a certain roughness at the downstream surface."
- p.884: "If the downstream surface differs from the upstream one only by its greater roughness (made specially for the purpose of a favorable asymmetry), then E_rd = E_ru, the asymmetry being due to a difference in the pre-exponential factors: α_d0 ≫ α_u0."
- p.883: "As a result, hydrogen is desorbed to the left or to the right depending only on relative surface barrier heights"
- p.885: "The solubility constant for tritium in Nb is unknown. But usually the isotope effect in hydrogen solubility is rather small (e.g., it is practically negligible for the solubility of H and D in vanadium [12])."

**Чего нет:** нет экспериментальных данных по S через разные поверхности и проверок W; всё — модельные следствия при W=1.

**Релевантность для W: высокая** — наиболее явная формулировка: детальное равновесие (k_r = zα/S²) с общей S для обеих сторон и утверждение, что выходной газ при балансе находится «в равновесии с растворённым водородом» p = c²/S²; отсюда p_d,max не зависит от выходной поверхности. Это прямо эквивалентно W=1 (принято как принцип, не измерено).



---

## 8. Livshits, Alimov, Notkin, Bacal 2002 (Appl. Phys. Lett.)

**Ссылка:** A.I. Livshits, V.N. Alimov, M.E. Notkin, M. Bacal, "Hydrogen superpermeation resistant to ion sputtering", *Appl. Phys. Lett.* **81** (2002) 2656–2658. doi:10.1063/1.1510961 (файл livshits_2002.pdf)

**Суть.** Трубчатая мембрана Nb (∅1 см, 18 см, 0.1 мм) в водородной плазме; энергия ионов задаётся смещением. Обычная деградация сверхпроницаемости при E_inc > ~50 эВ связана не с распылением O-монослоя как такового, а с полиатомным слоем NbC на «натуральной» поверхности: после декарбонизации (нагрев до 1550 °C, выход CO) сверхпроницаемость устойчива к распылению до ≥1 кэВ. Объяснение: O-покрытие обеих поверхностей находится в термодинамическом равновесии с растворённым O (Θ_OS^eq ≈ 1), равновесные вакансии O служат центрами десорбции H; симметрия выхода H сохраняется, если концентрации растворённого O у двух поверхностей равны и вакансии на входной поверхности остаются равновесными. NbC тормозит обмен O между поверхностью и решёткой.

**Числовые параметры**

| Величина | Значение | Единицы | Система | T | Где |
|---|---|---|---|---|---|
| Порог надтепловой энергии | ~1 | эВ | — | — | с.2656 |
| Доля потока через симметричную мембрану | ~1/2 | — | одинаковые поверхности | — | с.2656 |
| T нагрева до/после | ≤1280 °C (до); 1550 °C (декарбонизация); 1270 °C — максимум без выхода CO | °C | Nb | — | с.2656–2657 |
| T разложения NbC (реакция с растворённым O) | ~1400 | °C | Nb | — | с.2657 |
| Плазма | n ~10^11 см⁻³; T_e ~1 эВ; p_H2 = 3×10⁻³ Торр; примеси ~1×10⁻⁶ Торр (CO) | — | — | — | с.2657 |
| Поток CO при 1550 °C | ~3×10^15 CO/(см² с), спад ~5× за ~10 мин | — | Nb | 1550 °C | с.2657 |
| Удалено C (и O) | ~0.2 ат.% (~800 монослоёв) | — | Nb | — | с.2657 |
| Нанесённый слой C | ~10^18 C/см² (~10³ монослоёв); карбидизация при 1100 °C | — | Nb | 700 → 1100 °C | с.2657 |
| Устойчивость к распылению | до ≥1 кэВ (декарбонизированная) | эВ | Nb | — | с.2656 |
| Падение при наличии NbC | на порядок при E_inc > 50 эВ; за ~5 с | — | Nb | — | с.2657–2658 |
| Растворимость O в Nb | ~1 ат.% | — | Nb | рабочие T | с.2657 |
| Θ_OB (отн. объёмная конц. O) | 0.001–0.01 | — | Nb | — | с.2658 |
| Θ_OS^eq | ≈1 | — | Nb | все рабочие T | с.2657 |
| j_inc·γ | ≈1×10^14 | см⁻² с⁻¹ | Nb | >550 °C | с.2658 |
| ν | ≈10^13 | с⁻¹ | — | — | с.2658 |
| E_OD (диффузия O в Nb) | 115 | кДж/моль | Nb | — | с.2658 |
| ΔH_segr (O) | 72 | кДж/моль | Nb | — | с.2658 |
| D(O в Nb) | ≈10⁻⁹ | см²/с | Nb | 550 °C | с.2658 |
| Время удаления 10 % объёмного O при Θ_OB = 0.002 | ~1000 | с | Nb | — | с.2658 |

**Цитаты**

- p.2656: "The monolayer gives rise to a potential barrier impeding both the dissociative absorption and associative desorption of thermal H2 molecules."
- p.2656: "As a result, if the two membrane surfaces are identical (most typical case), ~1/2 of the implanted flux passes through the membrane, irrespective of its temperature and thickness."
- p.2656: "With E_inc exceeding the oxygen sputtering threshold, oxygen monolayer at the upstream surface gets damaged, which brings about an asymmetry of hydrogen release in favor of re-emission through the upstream surface."
- p.2657: "When there is no sputtering, the O coverage, Θ_OS, of the two membrane surfaces is in equilibrium with the O dissolved in the lattice: Θ_OS = Θ_OS^eq."
- p.2658: "Still, there is a certain equilibrium concentration of O vacancies in the O monolayer, Θ_VS^eq = 1 − Θ_OS^eq ≪ 1, at any given T and Θ_OB. These vacancies apparently serve as centers of thermal hydrogen desorption. The symmetry of hydrogen release is not violated if the approximate equality of O vacancy concentrations at the two surfaces is not perturbed by the upstream side sputtering."

**Чего нет:** S, k_r, α численно не приводятся; детальное равновесие для H не обсуждается; опытов с H2 с обеих сторон нет.

**Релевантность для W: низкая** — асимметрия здесь обсуждается только как кинетическая (концентрация центров десорбции на двух сторонах); косвенно полезно, что в равновесии (без распыления) обе поверхности Nb с растворённым O стремятся к одинаковому равновесному состоянию (симметрия), т.е. «естественная» асимметрия при термодинамическом равновесии примесей не возникает.



---

## 9. Livshits, Alimov, Notkin, Bacal 2004/2005 (Appl. Phys. A)

**Ссылка:** A.I. Livshits, V.N. Alimov, M.E. Notkin, M. Bacal, "Hydrogen superpermeation resistant to ion sputtering", *Appl. Phys. A* **80** (2005) 1661–1669 (online 23.03.2004). doi:10.1007/s00339-004-2595-0 (файл livshits2004.pdf)

**Суть.** Развёрнутая версия письма APL 2002. Трубчатая мембрана Nb (99.95 %, ∅1 см, 18 см, 0.1 мм) в плазме (~2 мА/см² ионов). Показано: (i) подповерхностный слой NbC/Nb2C ответственен за чувствительность сверхпроницаемости к распылению; (ii) после декарбидизации — устойчивость до 1 кэВ при T ≥ 520 °C; (iii) проницаемость для молекулярного H2 (плазма выключена) сильно зависит от концентрации растворённого O (после деоксидации при 2150 °C — рост на порядки, после легирования 1.8 ат.% O — падение на порядки), тогда как плазменная сверхпроницаемость от O не зависит. Модель: центры диссоциативной адсорбции/десорбции H2 — равновесные O-вакансии в O-монослое, α_eq = A(T)(Θ_VS^eq)^n = B(T)/Θ_OB^n. Явно сказано, что принцип детального равновесия связывает константы выхода с α на соответствующих поверхностях; асимметрия A определяется отношением α двух сторон; при отсутствии распыления термодинамическое равновесие поверхность–объём для O делает обе поверхности одинаковыми (α_up ≈ α_down).

**Числовые параметры**

| Величина | Значение | Единицы | Система | T | Где |
|---|---|---|---|---|---|
| Плазма | n ~10^11 см⁻³; T_e ~1 эВ; p_H2 ~3×10⁻³ Торр; объём 0.1 м³; примеси ~1×10⁻⁶ Торр (CO) | — | — | — | с.1662 |
| Ионный поток | ~2 мА/см² (H2⁺, H3⁺) | — | Nb | — | с.1662 |
| Порог распыления O-монослоя | ~50 эВ; деградация при 50–100 эВ; падение на порядки при 100–200 эВ | эВ | Nb | 550–750 °C | с.1661–1663 |
| Минимум вблизи | ≈10 эВ (хим. распыление O) | эВ | Nb | — | с.1662 |
| Нагрев, давший устойчивость | 1450 °C (несколько минут); ранее ≤1280 °C | °C | Nb | — | с.1661–1663 |
| Выход CO при 1450–1550 °C | ~10^18 | см⁻² | Nb | — | с.1663 |
| Предел растворимости C в Nb | 0.1 ат.% (~1400 °C); <0.1 ат.% при T<1400 °C | — | Nb | — | с.1663, 1668 |
| Толщина подповерхностного карбида | порядка 0.1 мкм | — | Nb | — | с.1663 |
| Нанесённый C | ~10^18 C/см²; карбидизация при 1200 °C | — | Nb | 700 → 1200 °C | с.1663 |
| Доп. O (легирование) | 1 и 2.5 ат.% (рис.2); 1.8 ат.% (рис.5–6) | ат.% | Nb | — | с.1662, 1664 |
| Деоксидация | 2150 °C | °C | Nb | — | с.1664 |
| p_in в опыте рис.5 | 4 мТорр H2 | — | Nb | 10³/T = 0.65–1.05 K⁻¹ | рис.5 |
| Проникновение H2 (молекулы), рис.5 (на глаз) | после деоксидации: ~3×10^15 (10³/T≈0.85) … ~5×10^14 (≈1.35) H2/(см² с); после легирования 1.8 % O: ~1.5×10^14 (≈0.65) … ~2×10^13 (≈1.0) | H2/(см² с) | Nb | — | рис.5 |
| То же, плазма (плавающий потенциал) | ~3×10^16, не зависит от T | H2/(см² с) | Nb | — | рис.5 |
| Оценка α по ур.(3) из рис.5 (моя, не авторов: j_inc(4 мТорр, 293 K) ≈ 6×10^18 H2/см²с) | ~5×10⁻⁴ (деоксидир., ~1000 K); ~6×10⁻⁶ (1.8 % O, ~1000 K) | — | Nb | ~1000 K | рис.5 (расчёт мой) |
| α чистого Nb | ≈0.1 | — | Nb | — | с.1666 |
| Θ_OB | 0.001–0.01 | — | Nb | — | с.1666 |
| N_me | 5×10^22 | см⁻³ | Nb | — | с.1667 |
| N_S | 10^15 | см⁻² | Nb | — | с.1667 |
| D*_O | 0.015 | см²/с | O в Nb | — | с.1667 |
| E_OD | 115 | кДж/моль | O в Nb | — | с.1667 |
| ν | ≈10^13 | с⁻¹ | — | — | с.1667 |
| ΔH_segr (O) | 42 (Hörz) и 71 (Farrell) | кДж/моль | O на Nb | — | с.1667 |
| j_inc γ | ≈1×10^14 | см⁻² с⁻¹ | — | — | с.1667 |
| E_th(Nb из NbC) | ≈300 эВ/протон (≈600 и 900 эВ для H2⁺, H3⁺) | эВ | NbC | — | с.1668 |
| Время распыления монослоя O | ~10^15/(j_inc γ) — несколько секунд | с | — | — | с.1668 |

**Ключевые формулы:** j_perm = j_inc ξ α_down/(α_up+α_down) (2); для молекул j_perm^H2 = 0.5 j_inc^H2 α (3, симметричная мембрана); Θ_OS^eq/(1−Θ_OS^eq) = Θ_OB/(1−Θ_OB)·exp(ΔH_segr/kT) (4); Θ_VS^eq = (1/Θ_OB) exp(−ΔH_segr/kT) (5); α_eq = A(T)(Θ_VS^eq)^n = B(T)/Θ_OB^n (6). (В ур.(1) напечатано A = α_up/(α_up+α_down), что не согласуется с (2) — по-видимому, опечатка; правильно α_down/(α_up+α_down).)

**Цитаты**

- p.1661: "Such a monolayer gives rise to a potential barrier impeding both the dissociative absorption and associative desorption of hydrogen in the form of thermal H2 molecules [1]."
- p.1661: "The membranes superpermeable to energetic hydrogen are at the same time virtually impermeable to thermal hydrogen molecules and to other gases. Because of this, the superpermeable membranes are able to automatically compress by orders of magnitude the permeating hydrogen"
- p.1665: "If there is no sputtering, a thermodynamic equilibrium between O-vacancies at the surface and O atoms at the surface and in the bulk is maintained. Therefore the inlet and outlet surfaces have equal concentration of O-vacancies, the membrane exhibiting a symmetry with respect to hydrogen evolution at either side."
- p.1665: "Under the experimental conditions of this work, hydrogen reemission by the metal virtually occurs only by way of the usual thermal desorption, even though the membrane is acted upon by energetic hydrogen gas [1, 19], and therefore A is determined by ratio of rate constants of the thermal release of absorbed hydrogen through the up- and downstream surfaces. The detailed balance principle links these rate constants with the probabilities of dissociative absorption of thermal H2 molecules (sticking coefficients) at the upstream (α_up) and downstream (α_down) surfaces."
- p.1665: "In the case of natural surfaces, a membrane typically is nearly symmetric, i.e., α_up ≈ α_down ≡ α."
- p.1665: "In the case of an interaction with thermal molecules, one should take ξ = α, and thus the flux of molecule driven permeation … through a symmetric membrane is given by: j_perm^H2 = 0.5 j_inc^H2 α … Thus the value of α can be easily found from measuring the molecule driven permeation."
- p.1666: "Then, the sticking coefficient, α_eq, for the dissociative absorption of H2 molecules through an equilibrium O coverage must be proportional to the equilibrium concentration of active centers … Thus, if the above idea about the role of surface O vacancies in the kinetics of hydrogen absorption-desorption is true, the molecule driven permeation, j_perm^H2, is expected to depend on concentration of solute O. The data presented in Fig. 5 demonstrate such a dependence"
- p.1664: "It should be noted that O content has no effect on plasma driven permeation (superpermeation) in both the sputter-resistant (Fig. 6) and sputter-sensitive states (Fig. 2), when there is no sputtering (e.g. at a zero bias)."

**Чего нет:** S по граням не измеряется; детальное равновесие используется, но не проверяется; опыта с H2 при равных давлениях с двух сторон нет. Подповерхностный слой здесь — NbC/Nb2C (для O, не как подповерхностное состояние H); через NbC толщиной ~1 мкм H проходит быстро (ссылка), сверхпроницаемость не нарушается.

**Релевантность для W: средняя** — явно опирается на детальное равновесие (константы выхода ↔ α той же стороны, при общей S), а также даёт важный косвенный довод: асимметрию α двух сторон авторы связывают с неравновесием (распылением); в термодинамическом равновесии O-покрытие и центры десорбции на обеих сторонах одинаковы. Это довод о кинетической симметрии, но о W прямо ничего не говорит.



---

## 10. Алимов, Буснюк, Кузенов, Передистов, Лившиц 2023

**Ссылка:** В.Н. Алимов, А.О. Буснюк, С.Р. Кузенов, Е.Ю. Передистов, А.И. Лившиц, «Рециркуляция изотопов водорода сверхпроницаемыми мембранами в термоядерных установках: мембраны из ванадиевого сплава с пониженной растворимостью водорода», в сб. (раздел «Кинетика и термодинамика взаимодействия…», Секция 1), 2023, с. 10–37. DOI: 10.53403/9785951505378_2023_10. УДК 538.931. (Полное название сборника и издательство на страницах файла не указаны — по DOI это издание ISBN 978-5-9515-0537-8; уточнить по DOI.)

**Суть.** Описание новой UHV плазменно-мембранной установки СУПМЕМ-1 и опыты с трубчатой мембраной из сплава V–6.5 ат.% Pd (150 мкм) при взаимодействии с молекулами H2 (MDP), тепловыми атомами H с Ta-нитей (ADP) и холодной плазмой (PDP). Проницаемость для H2 экспоненциально зависит от T; из неё найдена α_u = 0.012 exp(−4420/T) в ПРЕДПОЛОЖЕНИИ идентичности входной и выходной поверхностей (α_u = α_d). ADP и PDP не зависят от T (сверхпроницаемость). Декарбонизация (~1350–1500 °C) и растворение ~1 ат.% O делают мембрану устойчивой к распылению. Теоретическая часть повторяет соотношения Livshits–Hatano–Watanabe 2002: p_d,max = j/(zα_u); c = K√p_d·√[(p_d,max/p_d + α_d/α_u)/(1+α_d/α_u)] — с единой константой растворимости K для обеих сторон (W=1 по построению).

**Числовые параметры**

| Величина | Значение | Единицы | Система | T | Где |
|---|---|---|---|---|---|
| Порог надтепловой энергии | 1 | эВ | — | — | с.10 |
| Химическая энергия атома H относ. H2 | 2.25 | эВ | — | — | с.11, 21 |
| Доля прохождения тепловых атомов | ~20 % | — | — | — | с.11 |
| Предел j для SPM | ~10^19 | см⁻² с⁻¹ | — | — | с.11 |
| Условие сверхпроницаемости | j ≪ K²D²/(L² z α_u) (1); j ≪ (K*D*)²/(L² z α_u*) · exp[2(E_ru−E_b)/RT] (5) | — | — | — | с.14–15 |
| K_V | 2×10⁻⁶ exp(33 кДж/моль /RT) | (H/M)/Па^0.5 | H в V | — | с.17, ур.(18) |
| D_V | 3.5×10⁻⁴ exp(−5.0 кДж/моль /RT) | см²/с | H в V | — | с.17, ур.(19) |
| ΔH_V | −33 | кДж/моль | H в V | — | с.18 |
| α_u (V с растворённым O) | ≈10⁻⁴ (типично) | — | V + O | 400–700 °C | с.17 |
| j-ограничение (чистый V, L = 0.01 см) | ≪5×10^19 (600 °C); ≪2×10^19 (700 °C) | H2/(см² с) | V | — | с.17, (20)-(21) |
| j-ограничение (сплав, K/3, D/2) | ≪1.4×10^18 (600 °C); ≪1.0×10^18 (700 °C) | H2/(см² с) | V–Pd | — | с.18, (22)-(23) |
| ИТЭР-оценка | J = 170 м³Па/с = 4.5×10^22 (D2/T2)/с; S_SPM = 45 м² при j = 10^17; p_d = 1 Торр; G = 3.0 г (500 °C), 1.7 г (600 °C) | — | V 0.01 см | — | с.17 |
| q_Pd (эмпир. коэф.) | 1087 (в тексте единицы «кДж/(моль·ат.%)», но из ур.(25) 1087/R ≈ 130 следует Дж/(моль·ат.%)) | Дж/(моль·ат.%) | V–Pd | — | с.18–19 |
| K_{V–κPd}/K_V | exp(−1087 κ/(R T)) ≈ exp(−130κ/T) | — | V–Pd | — | с.19, ур.(25) |
| Снижение растворимости в ~3 раза | при ≈7 ат.% Pd | — | V–Pd | — | с.19, рис.3 |
| Нагрев установки | до 400 °C | — | — | — | с.20 |
| Ta-катоды | 5 шт., ∅0.8 мм, 150 мм | — | — | — | с.23 |
| Разряд | 50–70 В (82 В в опыте), до 10 А (2.5 А в опыте), p ≥ 0.1 Па; n ~10^10–10^11 см⁻³; T_e ~0.6 эВ | — | — | — | с.23–24, 28 |
| Мембрана | V–6.5 ат.% Pd; L = 121 мм; ∅6.1 мм; стенка 150 мкм | — | — | — | с.24 |
| Растворение оксидной плёнки V | при T ≥ 500 °C | — | V | — | с.24 |
| St. 1 | прогрев в UHV до 1000 °C | — | — | — | с.24 |
| Растворение O | p_O2 ~10⁻³ Па при ≈900–910 °C; введено ≈0.15 ат.%, затем +0.4 и +0.9 ат.%; итог St. 2: 1–1.3 ат.% O | — | V–Pd | — | с.25–27 |
| Декарбонизация | резкий рост выхода CO при ≈1100 °C; нагрев до 1350 °C; удалено ≈1.0×10^19 CO (≈0.5×10^18 C/см²); на рис.9: Q_des = 0.49×10^19 CO, «0.02 % C» | — | V–Pd | — | с.25–26 |
| p_in (MDP, ADP, PDP) | 3.5×10⁻³ Торр (3.47–3.54×10⁻³) | Торр | — | — | с.27 |
| Откачка | S_out = 4 л/с; S_in = 80 л/с (рис.10), 120 л/с (рис.11) | л/с | — | — | рис.10–11 |
| T Ta-нитей при MDP | ≈500 °C (ниже порога диссоциации) | °C | — | — | с.27 |
| α_u (из MDP, при α_u = α_d) | 0.012 exp(−4420/T); по рис.10 ~2×10⁻⁵ … 2×10⁻⁴ при 400–1000 °C; ~10⁻⁴ при 500–900 °C | — | V–6.5Pd, St.2 | 400–1000 °C | с.27–28, ур.(26), рис.10 |
| Формула извлечения α_u | j_p = 0.5 α_u z P_u | — | — | — | с.27 |
| T нитей при ADP | 1950 | °C | Ta | — | с.28 |
| Отношение площадей нитей / стенки | ~1/200 | — | — | — | с.28 |
| j_p (PDP) достигнутое | ≈1.5×10^16 | H2/(см² с) | V–6.5Pd | 650 °C, α_u = 1×10⁻⁴ | с.29, рис.11 |
| α_u для Nb и Ta с 1 ат.% O | 1×10⁻⁵ (Nb), 2×10⁻⁶ (Ta) | — | Nb, Ta | 700 °C | с.29 |
| Ионный ток при смещении | ~1 мА/см² | — | — | 650–750 °C | рис.12–13 |
| Минимум PDP / максимум / падение | −6 В / −30…−50 В / резкое падение при больших смещениях; порог физ. распыления O ≈ 50 эВ | В | St.1 | 650–750 °C | с.31–32 |
| Время установления PDP после смены смещения | ≈10 | с | St.1 | — | с.30 |

**Ключевые формулы (с.15–16):** p_d,max = j/(zα_u) (6); c = K√p_d · √[(p_d,max/p_d + α_d/α_u)/(1 + α_d/α_u)] (9); при α_d/α_u ≫ 1: c = K√p_d·√(p_d,max α_u/(p_d α_d) + 1) (11); при α_d/α_u ≫ p_d,max/p_d (12): c = K√p_d (13) — т.е. при сильной асимметрии растворённый H находится в равновесии с ВЫХОДНЫМ газом через ту же константу K.

**Цитаты (рус., оригинал)**

- с.13: «Барьер Eru является труднопреодолимым препятствием для абсорбции тепловых молекул Н2, и мембрана практически непроницаема для них.»
- с.13: «В результате абсорбированные атомы H «забывают» свое исходное положение и десорбируются налево или направо в зависимости только от относительной высоты поверхностного барьера. Если Eru > Erd, то большая часть имплантированных частиц водорода пройдет через мембрану независимо от температуры и толщины мембраны.»
- с.14: «…барьер отсутствует в случае атомарно-чистой поверхности переходных металлов (Eru = 0) [3, 26–28], но появляется на «реальной» металлической поверхности, т. е. на поверхности, покрытой монослоем неметаллической примеси … Если температура не слишком низкая, этот неметаллический монослой находится в динамическом равновесии с примесью, растворенной в объеме металла (явление поверхностной сегрегации)»
- с.15: «Если при стационарном проникающем потоке полностью прекратить отвод водорода с выходной стороны СПМ, его давление на выходной стороне будет расти до значения [25] p_d,max = j/(zα_u), (6) при котором обратный поток термолизованного молекулярного водорода, проникающего сквозь СПМ с выходной стороны, сравняется с прямым потоком.»
- с.27: «Левая ось на рис. 10 – это температурная зависимость вероятности диссоциативной абсорбции молекул Н2, αu, найденная из данных по плотности проникающего потока как j_p = 0,5 αu z P_u [1–3] в предположении, что входная и выходная поверхности идентичны, т. е. αu = αd.»
- с.32: «Предполагается, что причиной падения PDP под действием ионного распыления является нарушение монослойного покрытия О и, как следствие, нарушение симметрии мембраны, т. е. условия αu ≈ αd, в сторону преимущественного выделения поглощенного водорода обратно через входную поверхность, т. е. αu > αd.»
- с.33: «…распыление не полностью удаляет атомы О с поверхности, а генерирует лишь определенное количество кислородных вакансий в дополнение к определенному количеству равновесных вакансий, всегда присутствующих в кислородном монослое при данной температуре и полном количестве О в образце. Поскольку кислородные вакансии служат активными центрами диссоциативно-ассоциативных процессов при поглощении/выделении водорода…»

**Чего нет:** измерений растворимости/проницаемости через разные грани или в двух направлениях; проверки детального равновесия; опыта с H2 при равных давлениях. Газ — комнатной температуры (атомы термализованы на водоохлаждаемых стенках), мембрана 400–1000 °C — неизотермические условия.

**Релевантность для W: средняя** — использует W=1 (единая K) в теоретических соотношениях и измеряет α_u только в предположении симметрии; новых независимых свидетельств о W нет, но дан полезный набор свежих чисел (α_u(T), K, D для V и V–Pd).


---

# Сводка по пакету A (Лившиц и соавт.)

## Общее наблюдение
Ни в одной из 10 работ нет **прямого измерения** ни растворимости/константы Сивертса через разные грани, ни проницаемости для теплового H2 в двух противоположных направлениях через асимметричную мембрану при изотермических условиях, ни проверки детального равновесия в форме α = S²k_r/Z, где S, k_r и α измерены независимо. Также нет наблюдения потока или перепада давлений при молекулярном H2 одинакового давления с обеих сторон. Всюду W = 1 **заложено в модель**: одна S для обеих поверхностей, k_r = αZ/S². Иногда оно прямо названо следствием «detailed balance principle» (2002 FST, 2004) или требованием термодинамики (1990).

## За W = 1
1. **1990, с.82–83.** Прямо формулируется термодинамическое требование: «no regular gas flow may exist under isothermic conditions, if hydrogen pressures are equal on the both sides of a membrane». Утверждается, что проницаемость для равновесного молекулярного H2 «is the same in either direction» в обоих предельных режимах (10a, 10b). Авторы критикуют работы [37,38] (Swansiger & Bastasz 1979; патент Behr et al.) и Ионова [34]: из таких рассуждений следовал бы ненулевой поток при равных давлениях. Асимметрия по α даёт рост проницаемости для H2 максимум в 2 раза, без направленного эффекта.
2. **1990, с.83.** Асимметрия проницаемости для H2 (Pd с тонким Cu-покрытием) ожидается только тогда, когда газ «холоднее» мембраны. «The asymmetry is to disappear on equalizing the gas and the metal temperature».
3. **1991 (сноска 13), 1998 (с.2563).** Обратный молекулярный поток записан как 2Zp(1/α_0 + 1/α_L)⁻¹ = Zpα_0α_L/(α_0+α_L), то есть симметричен по перестановке α_0 ↔ α_L. При W=1 асимметрия по α одна потока при равных давлениях не даёт.
4. **2002 FST (с.883–884), 2023 (с.15).** При максимальном сжатии «gas at the downstream side is in equilibrium with absorbed hydrogen: p_d,max = c_b²/S²» и «with using the detailed balance principle, k_ru = zα_u/S²». Отсюда p_d,max = j_iξ/(zα_u) не зависит от выходной поверхности. По сути это утверждение W = 1.
5. **Экспериментальные косвенные данные.** Асимметричные (сверхпроницаемые) мембраны дают поток против перепада давлений и сжатие в 10³–10⁴ раз ТОЛЬКО при подводе неравновесного водорода: ионов, атомов с раскалённого атомизатора (~1900–2170 K), плазмы. 1997, с.1205: «No permeation could be noticed when hydrogen was admitted at the upstream side with a cold atomizer». 1998, с.2561: проникновение от одних молекул «very small» (ΔP_c ≈ 0.07 мТорр при p_up = 1.6 мТорр, 800 °C), а сжатие в ~10³ раз появляется только с плазмой. Ни в одном опыте не отмечено сжатия или потока от одного молекулярного H2. Оговорка: порог обнаружения не указан, и задача так не ставилась.
6. **2002 APL, 2004 (с.1665), 2023 (с.33).** Физико-химический довод: без распыления О-монослой на обеих поверхностях Nb (V) находится в термодинамическом равновесии с растворённым O. Поэтому концентрации О-вакансий (центров диссоциации/рекомбинации H2) на входе и выходе равны и мембрана симметрична. Асимметрию α_up ≠ α_down создаёт только неравновесное воздействие (распыление входной стороны). Это довод о симметрии кинетики в равновесии, а не о W напрямую. Но он согласуется с тем, что «естественная» асимметрия в равновесии не возникает.

## Против W = 1
- **В пакете ничего нет.** Ни одна работа не сообщает и не предсказывает поток при равных p и T. Единственная оговорка в 1990, с.82: различия проницаемости в разных направлениях «not entirely prohibited … if more complicated mechanisms of hydrogen-metal interaction are involved». Но речь там идёт о мембране между газом и вакуумом или газами разных давлений, и сразу добавлено требование отсутствия потока при равных давлениях. Сноска о различии «up to about 60%» в переходном режиме между (10a) и (10b) [6,36] тоже относится к ненулевому перепаду давлений (направленная разница при p ≠ 0 с одной стороны и вакууме с другой), а не к равным давлениям. Это не довод против W = 1.

## Нейтрально или методические предостережения
- **Асимметрия по α велика и хорошо документирована, но к W не относится напрямую.** α_L/α_0 ≈ 10 под плазмой (1998); α_d0/α_u0 ≈ 10 принимается как «разумная» (2002 FST); Nb–Pd с α_Pd > α_Nb (1990, 1991); пассивация C, S, O меняет α на порядки: α(C-Pd) = 0.15 exp(−16.8 кДж/моль/RT), α(S-Pd) = 0.1 exp(−34.4 кДж/моль/RT), α(Nb+1.5 % O) = 5.5×10⁻⁴ exp(−41.8 кДж/моль/RT), α(V–6.5Pd) = 0.012 exp(−4420/T). Всё это асимметрия K1 (и k_r) при общей S, то есть при W = 1 по построению.
- **α из опытов с H2 всегда извлекается в предположении симметрии** (α_u = α_d; χ = α/2; 1990, 2004, 2023). Значит, W эти измерения не проверяют.
- **Условия опытов неизотермичны.** Газ комнатной температуры (Z «at room temperature», 1998), мембрана 500–1000 °C. Авторы сами указывают (1990, с.83), что в таких условиях асимметрия проницаемости для H2 возможна и законна.
- **Подповерхностные состояния.** Явно упоминаются только O/C в приповерхностном слое и NbC/Nb2C (2002 APL, 2004), а также ловушки радиационных дефектов (1991). Подповерхностные состояния H как отдельная ступень кинетики в пакете не вводятся.
- **Порядок десорбции** везде второй (j_r = k_r c²). Для аустенитной стали оговорён возможный первый порядок (1990, с.82).
- **Степень сжатия.** g_max = βξ/α_0 (1991); P_cm = Σξj/(Zα_0) (1998); p_d,max = j_iξ/(zα_u) (2002, 2023). Достигнуто 10³–10⁴ (1992, 1997, 1998). Везде определяется только α входной стороны: это прямое следствие W = 1 в модели.
- **Теплоты и энергии.** ΔH (кДж/моль H): Fe 24, Ni 17, Pd −10, Nb −39 (−39.7), V −33, ASS 8.8 (1990; 2002; 2023). E_b: Nb, Ta ≈ −0.3 эВ, Pd ≈ +0.15 эВ (1995). Хемосорбция H на Nb: −56.4 (Pick) и −75 кДж/моль H (1991). ΔH_segr(O/Nb) = 42…72 кДж/моль (2002, 2004).

## Итог для вопроса о W
Пакет Лившица единодушно и явно (как постулат термодинамики и детального равновесия) стоит на W = 1. Экспериментально он показывает, что все наблюдавшиеся потоки против градиента и всё сжатие требуют внешнего неравновесного источника. Независимой количественной проверки W (измерения S через разные грани или потока при равных p и T) в пакете **нет**.
