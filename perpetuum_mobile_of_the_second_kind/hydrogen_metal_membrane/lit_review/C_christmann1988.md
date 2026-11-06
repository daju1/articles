# Пакет C: Christmann 1988 (обзор) — что он говорит о W = (S_A/S_B)²

Файл: `tmp/lit/src/hydrogen_metal_iteraction/christmann1988.pdf` (163 стр., текстовый слой есть; прочитан целиком, включая список литературы). Таблицы 3, 4, 9 и рис. 54 в текстовом слое повёрнуты или искажены, поэтому я отрендерил стр. 20, 21, 26, 27, 86 и 116 и прочитал их с изображений. Номера страниц ниже — журнальные; в этом PDF они совпадают с номерами страниц файла.

---

## 1. Библиографическая ссылка

K. Christmann, **"Interaction of hydrogen with solid surfaces"**, *Surface Science Reports* **9** (1988) 1–163. North-Holland, Amsterdam. DOI: 10.1016/0167-5729(88)90009-X. (Institut für Physikalische Chemie der Freien Universität Berlin; manuscript received in final form 13 April 1988.)

## 2. Суть (3–6 строк)

Большой экспериментальный обзор по хемо- и физисорбции H₂/H на монокристаллах металлов (и немного полупроводников). В нём разобраны: прилипание s₀ и его чувствительность к дефектам и ступеням; активированная и неактивированная диссоциация; теплоты адсорбции и энергии связи Me–H; адсорбционные места, колебания, диффузия, латеральные взаимодействия, 2D-фазовые переходы, реконструкции; кинетика ТПД (порядок, предэкспоненты), H₂–D₂-обмен, подповерхностный водород. **Объём металла, растворимость и проницаемость автор сознательно исключает** (стр. 4). Поэтому прямых данных о K₂ (поверхность⇄объём), о S или о сравнении растворимости через разные грани в обзоре **нет**. Для вопроса о W полезны три вещи: (а) явное использование детального равновесия на каждой грани для связи s₀ и ν_d (ур. 30–34); (б) локальный «закон Сивертса» для поверхности в равновесии адсорбции и десорбции (стр. 123); (в) качественная картина подповерхностных состояний, в которой поверхность меняет только барьеры, а объёмные уровни (теплота растворения) берутся как свойство объёма.

## 3. Ключевые дословные цитаты

### 3.1 Объём / растворимость исключены из обзора
- p. 4: *"…and, of course, heats of solution of hydrogen in bulk materials [13]. The problem of hydrogen interacting with bulk materials will, however, be ignored in this article because there exists a great many exhaustive review papers on this subject [4]."*

### 3.2 Детальное равновесие (detailed balancing) для адсорбции и десорбции
- p. 15: *"Particles leaving the surface again first recombine to H₂ molecules before they enter the gas phase. Actually there is no single experiment known which gave evidence for an atomic desorption. According to the detailed balancing concept the two dissociation reaction paths must show up also in the desorption path: The reversal of the indirect process (in which the molecules desorb from the precursor) will lead to a cosine distribution whilst directly desorbing molecules come off translationally "hot", i.e., in a sharper than cosine distribution. This was, among others, confirmed for the H/Ni(111) system, where the angular distribution of the desorbed molecules obeyed a cos⁵ power law [39,48,49]."*
- p. 19: *"Also, there are no anomalies with regard to the angular distribution of desorbing H₂ molecules … in all these cases a simple cosine distribution of the emitted H₂ molecules is observed. This means that the desorbing molecules are in thermal equilibrium with the surface, unlike molecules desorbing, e.g., from Ni(111) whose rate components are peaked normal to the surface up to cos⁵θ … As pointed out before, the concept of "detailed balancing" (for more details see, e.g., refs. [73,74]) would just predict this behavior. Only those particles come off translationally "hot", i.e., with an excess kinetic energy, which have passed an activation barrier for adsorption in the entrance channel, on their way from the gas phase into the adsorbed state, and vice versa."*
- p. 19 (Pd(100), десорбция из объёма — наблюдение, существенное для асимметрии путей): *"…Comsa et al. [75,76] who studied the velocity distribution of D₂ molecules desorbing from a Cu(100) or a Pd(100) surface by means of a time-of-flight technique and found also that part of the desorbing molecules showed fluxes that were strongly peaked normal to the surface. However, for the particular case of Pd there is no activation energy for adsorption of hydrogen as a detailed study of the H/Pd(100) system revealed [77] and can therefore not explain the observed effects. Instead, the fast molecules are absorbed in the Pd bulk, travel to the surface upon heating and recombine there directly, without equilibrating in the chemisorption well. Moreover, sulfur impurities played an important role in the penetration and recombination mechanism. This is just an example to demonstrate that care has to be taken when interpreting angular distributions of desorbing hydrogen."*
- p. 88: *"Other complications, as mentioned before, may arise from the existence of precursor states in the desorption path which, according to the "detailed balancing", should modify the thermal desorption traces considerably, as long as the temperatures are low enough to ensure a noticeable population of these states."*
- p. 89 (Kubiak et al., десорбция после **проникновения через кристалл Cu**): *"…Kubiak et al. [274] who studied the rotational and vibrational distribution of H₂ and D₂ recombinatively desorbing from clean Cu(110) and (111) surfaces following atomic permeation through the Cu crystals at T_s = 850 K and 960 K. … The rotational distributions of the desorbing H₂ (D₂) molecules were found to be non-Boltzmann and showed mean rotational energies which where 80%–90% of the surface temperature, T_s, regardless of the isotope and of the surface orientation. In contrast, the vibrational population ratio as given by the fraction of molecules with quantum numbers v₁ : v₀ was found approximately a factor of 100 greater than the ratio predicted by Boltzmann statistics at T_s. The conclusion was that detailed balance arguments applied to the vibrational distributions measured in recombinative desorption are not suited to predict correctly the dissociative adsorption probability as a function of vibration: It is concluded that these two processes are dynamically different at least for the Cu(110) and Cu(111)/H systems."*

### 3.3 Равновесие адсорбции и десорбции → связь s₀ и ν_d (аналог α = S²k_r/Z для поверхности)
- p. 84: *"It can be shown that, for the zero coverage limit, the pre-exponential factor is obtained by equating the rates of adsorption and desorption, r_ad and r_des, (which corresponds to the adsorption-desorption equilibrium)"* — ур. (30); далее ур. (31) H₂(g) + 2* → 2 H_ad, ур. (32) n_g/n_ad² = (f_g f_s,0 / f_ad² f_s,ad)·exp(−E_ad/kT) «for the case of non-activated adsorption … and the limit of small coverages and non-interacting particles».
- p. 85, ур. (34): *ν₂ = s₀ c̄ f_g / f′_ad²* («for the zero coverage limit ν₂ can be expressed by the initial sticking probability, s₀ and the partition function of the gas molecule, f_g, and the adsorbed phase, f_ad»). Здесь же оговорка, что при H-индуцированной реконструкции f_s,0 ≠ f_s,ad: *"This assumption, however, may be considerably wrong, if, e.g., the substrate reconstructs under hydrogen. Then, soft phonon modes may become excited and f_s,0 ≠ f_s,ad."*
- p. 85: *"It has been established that with hydrogen adsorption/desorption κ_d is of the order of 0.1, whereas with CO adsorption κ_d ≈ 1.0 [264]."*
- p. 86, ур. (38): *"As pointed out by Madix et al. [265] this yields for hydrogen ν₂ = 2.2 × 10²⁸ (f′_ad)⁻² [cm²/s·H-atom]."*

### 3.4 Поверхностный «закон Сивертса» в равновесии (раздел о H₂–D₂-обмене)
- p. 123: *"With the adsorption-desorption equilibrium being established we have d[H_ad]/dt = d[D_ad]/dt = 0, and we can express [H_ad] and [D_ad], respectively, by [H_ad] = ((k₁/k₋₁)[H₂(g)])^{1/2} and [D_ad] = ((k₂/k₋₂)[D₂(g)])^{1/2}."*
- p. 121 (равновесный состав не зависит от поверхности-катализатора): *"At equilibrium conditions and at temperatures above ~300 K the relative amounts of hydrogen, hydrogen deuteride and deuterium formed on active transition metal surfaces are related by 1 : 2 : 1 [394]…"*
- p. 125 (грань и ступени меняют скорость, но не энергию активации): *"Christmann and Ertl report on the same ΔE* = 10.5 kJ/mol for a stepped (997) and a flat Pt(111) surface, whereby the actual rate was a factor of 10 higher on the stepped surface [131]."* Подпись к рис. 55: *"Note that in both cases the same activation energy results whilst the pre-exponentials are different."*

### 3.5 Подповерхностные состояния и связь с объёмом
- p. 116: *"This subsurface state may be regarded as an intermediate between a chemisorbed and a dissolved hydrogen species and was proposed among others by Behm et al. [275] for hydrogen interacting with a Pd(110) surface. … This diagram is reproduced in fig. 54 and basically consists of a superposition of the potential energy curves of the dissociative chemisorption (cf. fig. 6) and of the dissolved hydrogen in the Pd bulk, whereby the respective bulk energies were taken from earlier work by Auer and Grabke [373] and Lynch and Flanagan [374]."*
- p. 116, подпись к рис. 54 (Pd(110)): *"…E_sol, heat of solution; E_diff, activation energy for bulk diffusion; E_ss, activation energy for desorption from the H subsurface sites; E_hydr, activation energy for hydrogenation; E_H, activation energy for desorption from the chemisorbed state; E_H₂, adsorption energy of the physisorbed H₂; E_diss, heat of dissociation of the H₂ molecule. The dashed line indicates the potential energy situation after the surface has been reconstructed by 1.5 monolayers of hydrogen. Clearly, the activation barrier for population of subsurface sites is removed by reconstruction."* (На самой диаграмме числовых значений нет — проверено по растру.)
- p. 117: *"Behm et al. pointed to the important influence of the surface (1×2) reconstruction of the Pd(110) surface at Θ_H = 1.5 ML which may lead to structural changes of the surface region of the Pd crystal so that further uptake of H may easily become possible which in turn causes the reconstruction to proceed."*
- p. 117 (Eberhardt et al., Ni/Pd/Pt(111)): *"…evidence for a (temperature dependent) occupation of subsurface sites with the (111) surfaces of Ni, Pd and Pt. … They specified their statements somewhat since it could be shown that at least for H/Ni(111) [50] and for Pt(111) [378] no such subsurface H species was detected by other methods (LEED, ion channeling)."*
- p. 118 (Ru(0001), Yates et al.): *"…in going from the vacuum side to the bulk metal a first deep potential populated by chemisorbed H is passed followed by a shallower endothermic subsurface potential which is then transforming to the endothermic bulk H sites. All potential minima are separated by appreciable activation barriers."* Там же: *"…a final conclusion is not yet reached."*
- p. 105 (Pd(111), EAM): *"Calculations of the classical diffusion barriers for hydrogen motion on a surface and between surface and subsurface sites reveal that the equilibrium is reached very rapidly, that is, within the time scale of the experiment."*
- p. 119 (Lagos): *"…the self-trapping energy of an interstitial impurity (e.g., a H atom) increases substantially when it approaches the surface, taking at the first layer of interstitial sites a value which is considerably larger than that of the bulk."*
- p. 119 (итоговая позиция автора): *"…whereever and whenever the small H atom adsorbs on a surface a strong perturbation of the local electronic structure occurs which also causes a local strain on the crystal lattice of the substrate. … in the most complicated case accompanied by a steady uptake (dissolution, permeation) of hydrogen leading to a progressive reconstruction and final hydride formation (as observed with Pd(110))."*
- p. 144 (Muscat, Ni(111)): *"Owing to its very low binding energy value (−9.5 eV below E_F) which is not found experimentally, a H subsurface site can be ruled out."*

### 3.6 Зависимость кинетики от грани и дефектов (кинетика, не равновесие)
- p. 5: *"…recent studies of the dynamical process of the H adsorption have revealed that there is one obstacle left for getting reproducible kinetic data, particularly for the zero coverage limit, namely, the characterization of kind and number of crystallographical defects…"*
- p. 18: *"…extremely carefully prepared Pt(111) surfaces revealed s₀ values as small as 10⁻⁴. … These data can be quantitatively understood in terms of a model in which the H₂ dissociation occurs exclusively at step sites."*
- p. 41 (заселённость ступеней по Больцману, ур. 10): *"…otherwise this population is governed by Boltzmann statistics according to N_steps/N_terrace = exp(−ΔE/kT) … which may wipe out the differences to a certain extent at higher temperatures."*

### 3.7 Что прямо искали, но в обзоре нет
Измерений растворимости или проницаемости через разные грани одного кристалла нет. Проницаемости в двух противоположных направлениях нет. Нет величин S, K₂, k_r в «мембранном» смысле, нет степени сжатия β. Нет наблюдений потока или перепада давлений при молекулярном H₂ одинакового давления с обеих сторон. Нет и явной проверки соотношения α = S²k_r/Z. Единственные «проникающие» эксперименты, которые упоминаются, — Kubiak et al. (десорбция после проникновения через Cu) и Comsa et al. (быстрые D₂ из объёма Pd(100)), оба в пересказе.

## 4. Числовые параметры

### 4.1 Прилипание, барьеры, теплоты (по тексту)

| Величина | Значение | Ед. | Система | T | Стр./ур. |
|---|---|---|---|---|---|
| Средняя тепловая скорость H₂ | ≈1775 | м/с | газ | 300 K | p. 11 |
| Типичный s₀ для H₂ | 0.01–0.1 | — | многие плотные грани | — | p. 17 |
| s₀ H₂/Pt(111), бездефектная | до 10⁻⁴ | — | Pt(111) | — | p. 18 |
| s₀ на «открытых» гранях | ≈1 | — | Ni(110), Pd(110), Co(10-10), Rh(110), Ru(10-10), W(111) | — | p. 19 |
| E_a адсорбции | ≈5 | кДж/моль | Ni(100) (Hamza & Madix, изотопный эффект) | — | p. 17 |
| Барьер (теория Jackson & Metiu) | ≈7 (ширина ≈0.6 Å) | кДж/моль | Ni(100), on-top | — | p. 17 |
| Барьер адсорбции (порог образования HD) | 10–20 (энергия пучка 6.7–45) | кДж/моль | Cu(110), Cu(310), Cu(100) | — | p. 29 |
| Теплота адсорбции (изотермы) | 39 → 41.5 | кДж/моль | Cu(311) | — | p. 28 |
| Энергия физисорбции H₂ | 3.5–15 | кДж/моль | разные | ~ЖГе | p. 23 |
| Глубина ямы хемосорбции атома H | ~1 | эВ | общее | — | p. 15 |
| Глубина атомной ямы | 500–600 | кДж/моль | «одна связь H–поверхность» | — | p. 23–24 |
| E_Me–H (типично) | ≈2.6 | эВ | переходные металлы | — | p. 25, 129 |
| Теплоты адсорбции H₂ (диапазон) | 50–150 | кДж/моль | хемосорбция | — | p. 148 |
| Энергии связи H–атом подложки | 240–300 | кДж/моль | — | — | p. 148 |
| ΔE_ad ступень − терраса | 5–10 | кДж/моль | Pd(111), Pt(111)/(997) | — | p. 41 |
| q_st | ≈85, рост на ≈12 | кДж/моль | Ni(110) (притяжение) | — | p. 72 |
| q_st − q_diff = RT | ≈2.5 | кДж/моль | — | 300 K | p. 71, ур. 19 |
| ΔE_diss(H₂) − ΔE_diss(D₂) | −8.1 (431.37 − 439.5) | кДж/моль | газ | — | p. 31, ур. 4 |
| Разность нулевых энергий ΔE₀ (H/D) | 1–3 | кДж/моль | H на переходных металлах | — | p. 32 |
| κ_d (трансмиссионный коэфф.) | ~0.1 (H₂); ≈1 (CO) | — | — | — | p. 85 |
| ν₂ (по Madix et al.) | 2.2×10²⁸ (f′_ad)⁻² | см²/(с·атом H) | — | — | p. 86, ур. 38 |
| Типичные ν_d 2-го порядка | 10⁻¹–10⁻³ | см²/(с·атом H) | разные | — | p. 87 |
| Время жизни: упругое / прямое неупругое / прекурсор / хемосорбция | ≤10⁻¹³ / ≤10⁻¹² / до 10⁻⁹ / 10⁻³–10³ | с | общее | — | p. 15 |
| E* обмена H₂–D₂ | 10.5 (одинаково); скорость ×10 на ступенчатой | кДж/моль | Pt(997) vs Pt(111) | — | p. 125 |
| E* обмена H₂–D₂ | ≈−40 (низкие T) / 15 (высокие T) | кДж/моль | Ni(111) (Schober) | — | p. 125 |
| Барьер газофазной реакции H₂ + D₂ → 2HD | 176 | кДж/моль | газ | — | p. 124 |
| Равновесное отношение H₂ : HD : D₂ | 1 : 2 : 1 | — | активные переходные металлы | >300 K | p. 121 |
| E_diff поверхностная | 14.6 (H) / 18.6 (D); D₀ = 2.5×10⁻³ / 8.5×10⁻³ см²/с | кДж/моль | Ni(100) | — | p. 66 |
| E_diff | ≈25; D₀ = 1.8×10⁻⁵ см²/с | кДж/моль | W | — | p. 66 |
| E_diff | 30 | кДж/моль | Ni (FEM-острие) | — | p. 66 |
| E_diff (физисорбированный T₂) | 0.84; D₀ = 3×10⁻⁴ см²/с | кДж/моль | Ni(111) | 15 K | p. 66 |
| Барьер миграции (кластеры, теория Siegbahn) | ≈22 | кДж/моль | Ni | — | p. 135 |
| Энергия связи H (EAM): поверхность / октаэдр. подповерхностный узел | 2.67 / 2.65 | эВ | Pd(111) | — | p. 142 |
| Тройной узел на поверхности стабильнее октаэдрического объёмного | 0.25 | эВ | Pd (кластер, LSD) | — | p. 146 |
| Барьер объёмной диффузии над окт. узлом | 0.34 (теор.) / 0.23 (эксп.) | эВ | Pd | — | p. 146 |
| d(Pd–H) теория | 1.72 | Å | Pd(111) тройной узел | — | p. 146 |
| Подповерхностный уровень | −9.5 эВ ниже E_F → подповерхн. место исключено | эВ | Ni(111) (Muscat) | — | p. 144 |
| Вклад d-электронов в энергию связи | ≈−0.4 (≈20 %); только jellium −1.7…−2.6 | эВ | Ni(111) | — | p. 144 |
| ΔE (EMT, Nørskov) | −2.7 (теор.) / −2.74 (эксп.) | эВ отн. своб. атома | Ni(100), Ni(111) | — | p. 141, табл. 12 |
| ħω₀ | 79 / 74 (Ni(100)); 134 / 139 (Ni(111)) теор./эксп. | мэВ | — | — | p. 141, табл. 12 |
| Энергия Si–H | 2.2 (Θ = 0.05) – 3.1 (Θ = 1) | эВ | Si(111)-7×7, атомарный H | — | p. 129 |
| Энергия обмена (Ising) | 3.81; либо w₁ = 2.13, w₂ = 1.07 | кДж/моль | Pd(100), T_c = 260 K | — | p. 101 |
| w₁, w₂, w₃, трио | −7.9, −1.3, −2.1, +3.35 | кДж/моль | Fe(110) | — | p. 102 |
| T_c упорядоченных фаз | 270 (Ni(111) c(2×2)); 260 (Pd(100)); 245/265 (Fe(110)); 85/105 (Pd(111) √3-фазы) | K | — | — | p. 99–103 |
| H–Ni и межслойное Ni–Ni | 1.84 / 1.15 ± 0.05 | Å | Ni(111) | — | p. 98 |
| Дипольные моменты H | ~0.05 Д (типично) | Д | табл. 6 | — | p. 48–49 |

### 4.2 Таблица 3 (s₀, только монокристаллы; прочитано по растру стр. 20–21)

| Грань | T, K | s₀ | Прим. |
|---|---|---|---|
| Fe(100) | 300 | 0.03 | |
| Fe(110) | 150 | 0.16 | |
| Co(0001) | 300 | 0.045 | |
| Co(10-10) | 120 | 1.0 | TDS |
| Ni(100) | 120 | 0.06–0.1 | TDS; рядом строка «activated adsorption» [57] |
| Ni(111) | 120 / 140 | 0.1; 0.05 | TDS; «activated adsorption» [39]; там же строка «~1» — по растру не удалось однозначно сопоставить с гранью |
| Ni(110) | 120 / 300 | 0.96 / 0.4 | TDS / MB |
| Cu(100), (111), (110) | — | < 0.1 | activated adsorption, MB |
| Mo(100) / (110) / (111) | 300 | 0.7 / 0.35 / 0.7 | |
| Ru(0001) | 100 | 0.4; 0.25; 0.3 | разные работы |
| Ru(10-10) | 100 | 1.0 | TDS |
| Rh(110) | 85 | 1.0 | TDS |
| Rh(111) | 100 | 0.65 | TDS |
| Pd(100) | 120 | 0.5 | TDS |
| Pd(110) | 120 | 0.7 | TDS |
| Pd(111) | 80 | 0.5 | LEED, TDS |
| Ta(100) | 300 | 0.04 | |
| W(100) | 300 | 0.65; 0.51; 0.6; 0.56 | одна — calibrated capillary |
| W(110) / (111) / (211) | 300 | 0.22 / 1 / 0.57 | |
| Ir(111) | 100 | 7×10⁻³ | |
| Ir(110) | 130 / 300 / 130 | 1 (β₂) / 1 / 0.007 (β₁) | |
| Pt(100)-hex | 78 / 300 / 300 | 0.15 / 0.17 / 0.07 | |
| Pt(110)-1×2 | 78 / 300 / 300 | 0.2 / 0.31 / 0.33 | |
| Pt(111) | 150 / 300 / 78 / 300 / 300 | 0.1 (activated) / 0.07 (MB) / 0.1 / 0.016 / 0.2 | |
| Pt(210) | 78 | 0.4 | |
| Pt(997) | 120 | 0.34 | TDS |

Nb в таблице 3 нет.

### 4.3 Таблица 4 (начальные теплоты адсорбции E_ad,0 и E_Me–H, кДж/моль; растр стр. 26–27)

| Грань | E_ad,0 | E_Me–H |
|---|---|---|
| Fe(110) / (100) / (111) | 109 / 100 / 88 | 271 / 265 / 260 |
| Co(0001) / (10-10) | 67 / 70 | 250 / 251 |
| Ni(100) / (110) / (111) | 95 / 90 / 95 | 264 / 261 / 264 |
| Cu(311) | 39 | 236 |
| **Nb(100)** | **111** | **272** |
| Mo(110) | 142 | 287 |
| Ru(0001) | 80; 120 (две работы) | 256; 276 |
| Rh(111) / (110) | 78 / 77 | 255 / 255 |
| **Pd(100) / (111) / (110)** | **99 / 88 / 102** | **266 / 260 / 267** |
| Ag(111) | ~15 | 241 |
| W(100) / (110) / (111) | 134 / 138 / 155 | 283 / 285 / 294 |
| Ir(110) / (111) | 96.2 / 53 | 264 / 243 |
| Pt(100)-hex / (110)-hex | 62.8 / 54.8 | 247.5 / 243 |
| **Pt(111)** | **42; 77; 75** (три работы) | 237; 255; 254 |

Для одного и того же металла теплоты адсорбции на разных гранях различаются на 10–40 кДж/моль: например, Pd 88–102, Fe 88–109, W 134–155.

### 4.4 Таблица 9 (ν_d второго порядка, см² с⁻¹ (атом H)⁻¹; растр стр. 86)

Fe(100) 1.5×10⁻²; Ni(100) 3×10⁰, 2.5×10⁻¹, 8×10⁻², 5×10⁻²; Ni(111) 2×10⁻¹, 2.3×10⁻²; Mo(110) 10⁻²; Mo(100) 5×10⁻² (β₂); Ru(0001) 3×10⁰ (state 1), 5×10⁻² (state 2), 4×10⁻³; Rh(111) 1.2×10⁻³; Pd(100) 1×10⁻²; Pd(110) 8×10⁻²; W(100) 4.2×10⁻² (β₂); W(110) 10⁻² (β₁), 1.4×10⁻² (β₂); Ir(111) 2×10⁻⁶; Ir(110) 1.5×10⁻² (β₂), 10⁻⁷ (β₁); Pt(111) 10⁻³, 3×10⁻⁹; Pt(110)-1×2 ~10⁻⁴ (β₂), 3×10⁻⁴ (β₁).

Pt(111) — известный спор (p. 87–88): разные группы получили 2-й порядок с низким ν и низким E_d* или 1-й порядок с высоким ν и E_d*. Пучковые данные Comsa поддержали более высокую E_d*. Концентрация дефектов «crucially influencing the adsorption/desorption kinetics».

### 4.5 Порядок кинетики (p. 78–79, 87–88)
- **Адсорбция 2-го порядка (Ленгмюр):** Rh(111), Ir(111).
- **Адсорбция 1-го порядка:** Ni(111) при 300 K, Ir(110)-1×2 (β₂), Pt(110)-1×2, Ru(0001) («N»-состояние). Объясняется островками или тем, что лимитирует захват в прекурсор.
- **Прекурсорная кинетика:** Pd(100) (рис. 28, 170 K), Ni(110) ниже 200 K, Pd(110), Ru(10-10), Co(10-10), Rh(110), W(100).
- **Десорбция 2-го порядка** (высокоэнергетические состояния): Ni(111), Ni(100), Pd(110), Pd(111), Rh(111), Pt(100), Fe(110), Ru(0001), W(110).
- **Псевдо-1-й порядок:** β₁ на Ni(110), α-состояния на Rh(110) и Pd(110). Для Pd(110)/H и Ni(110)/H кинетику десорбции модифицируют реконструкция и/или заселение подповерхностных узлов (p. 88).

## 5. Оценка релевантности для вопроса о W

**Средняя.** Обзор подробно показывает, что константа «газ ⇄ поверхность» K₁ (через s₀, ν_d, E_ad) сильно зависит от грани, и явно выводит её из детального равновесия на каждой грани (ур. 30–34, p. 123). Но о K₂ (поверхность ⇄ объём), растворимости, проницаемости и сравнении граней по S в нём ничего нет: эта тема намеренно исключена.

## 6. Сводка по пакету

### За W = 1 (косвенно)
1. **Детальное равновесие применяется как рабочий принцип на каждой грани.** Предэкспонент десорбции ν₂ выводится приравниванием скоростей адсорбции и десорбции в равновесии: ν₂ = s₀ c̄ f_g / f′_ad² (ур. 30–34, p. 84–85). Значит, K₁ каждой грани — термодинамическая величина: отношение статсумм × exp(E_ad/kT), а не независимая кинетическая. Различия s₀ между гранями (от 10⁻⁴ до 1) скомпенсированы различиями ν_d.
2. Угловые распределения десорбции на каждой грани согласуются с детальным равновесием: косинус на «открытых» гранях, cos⁵θ на Ni(111) с барьером во входном канале (p. 15, 19). Это прямая экспериментальная проверка микрообратимости для пары адсорбция/десорбция на каждой грани отдельно.
3. Равновесный состав H₂ : HD : D₂ = 1 : 2 : 1 одинаков на активных поверхностях (p. 121). Ступени и грань меняют только предэкспонент скорости (Pt(997) vs Pt(111): та же E* = 10.5 кДж/моль, скорость ×10, p. 125). Поверхность ведёт себя как катализатор: меняет скорость, но не равновесие.
4. На диаграмме подповерхностного водорода (Pd(110), рис. 54, p. 116) объёмные уровни (E_sol, E_diff) берутся из объёмных данных (Auer & Grabke; Lynch & Flanagan), а состояние поверхности (реконструкция) меняет **только барьер** входа в подповерхностные узлы. Это ровно та структура, при которой ΔG «газ → объём» не зависит от грани, а значит W = 1. Явного утверждения «растворимость не зависит от грани» в обзоре **нет**.
5. Для Pd(111) равновесие поверхность ⇄ подповерхность устанавливается быстро (EAM, p. 105), то есть промежуточные состояния термализованы.
6. Заселённость ступеней и террас подчиняется Больцману (ур. 10, p. 41): поверхностные неоднородности учитываются равновесной статистикой.

### Против W = 1 (то, что можно было бы так истолковать)
1. **Kubiak et al. (p. 89):** после проникновения H через Cu(110) и Cu(111) десорбирующие молекулы имеют избыток v = 1 примерно в 100 раз против Больцмана. Вывод в обзоре: аргументы детального равновесия «not suited to predict correctly the dissociative adsorption probability as a function of vibration», процессы «dynamically different». Это единственное в обзоре место, где детальное равновесие формально ставится под сомнение, и это эксперимент с проникновением через объём. Однако: (а) это неравновесный опыт (только десорбция при T_s = 850–960 K, без измерения адсорбции в тех же условиях); (б) ни потока, ни перепада давлений при равных давлениях не измеряли; (в) неравновесное колебательное распределение само по себе совместимо с детальным равновесием, если адсорбция сильно колебательно-активирована. *Примечание вне пакета, по памяти — проверить:* по-моему, позднее Rettner, Michelsen, Auerbach и соавт. (≈1992–1995, H₂/Cu(111)) показали, что состояние-разрешённые адсорбция и десорбция согласуются с детальным равновесием.
2. **Comsa et al., Pd(100) (p. 19):** быстрые D₂ «absorbed in the Pd bulk … recombine there directly, without equilibrating in the chemisorption well». Путь объём → газ не проходит через хемосорбционную яму. Это пример того, что путь через грань может отличаться от простой цепочки газ ⇄ поверхность ⇄ объём. Но для W важны константы равновесия, а не пути. По микрообратимости обратный канал (прямое поглощение горячих молекул) должен существовать. Кроме того, автор указывает на роль примесей серы.
3. Константы K₁ сильно различаются между гранями: E_ad,0 у Pd 88–102 кДж/моль, у Pt(111) 42–77 по разным работам; s₀ от 10⁻⁴ до 1; ν_d от 10⁻⁹ до 3. Это **не** аргумент против W = 1 (W = 1 требует постоянства произведения K₁K₂², а не K₁), но показывает, что «поверхностная» часть асимметрична. Из обзора не видно, компенсируется ли она K₂, — данных о K₂ нет.

### Нейтрально / нет данных
- Растворимость и проницаемость, S, K₂, k_r (в мембранном смысле), β, сравнение граней по S, проницаемость в двух направлениях, наблюдения потока или перепада давлений при равном давлении H₂ с обеих сторон — **отсутствуют**. Объём металла исключён из обзора (p. 4).
- Существование подповерхностных состояний: Pd(110), Pd(111), Ru(0001), Cu(110)/(111)/(100) — предложены или обсуждаются; для Ni(111) и Pt(111) не подтверждены (p. 116–119, 144). Отношение энергий поверхность/подповерхность различается по системам: Pd(111) 2.67 vs 2.65 эВ (EAM); у Ru подповерхностные узлы эндотермичны. Это данные о K₂ на отдельных гранях, но без количественного сравнения граней с объёмом.
- Реконструкция и мягкие фононы меняют f_s,ad ≠ f_s,0 (p. 85), то есть энтропийную часть K₁. Вывода о W отсюда нет.
- Ni(100) и Pt(111) и др.: таблицы s₀, E_ad, ν_d (разделы 4.2–4.4) полезны как входные параметры для модели W, но сами по себе W не определяют.

**Итог по пакету:** обзор Christmann 1988 не содержит ни прямого измерения W, ни прямого опровержения W = 1. Всё, что в нём есть о равновесии, построено на детальном равновесии на каждой грани, а поверхность в нём меняет кинетику (барьеры, предэкспоненты), но не равновесный состав. Единственное «тревожное» место (Kubiak, p. 89) касается состояние-разрешённой динамики в неравновесном опыте с проникновением и к асимметрии констант Сивертса отношения не имеет.


---

## Дополнение (2 нояб. 2026): проверка вывода Kubiak et al.


**Цитата сверена с PDF обзора (с. 89) дословно.** Ссылка [274]: G.D. Kubiak, G.O. Sitz, R.N. Zare, *J. Chem. Phys.* **83** (1985) 2538.

Опыт:
* атомарный водород проникал через кристаллы Cu(110) и Cu(111), десорбция в вакуум при $T_s$ = 850 и 960 K;
* вращательная энергия молекул — 80–90 % от $kT_s$;
* отношение заселённостей $v{=}1$ : $v{=}0$ примерно в 100 раз больше больцмановского при $T_s$.

По аннотации к Kubiak et al., *J. Chem. Phys.* **81**, 6397 (1984): превышение примерно в 50 раз для Cu(110) и примерно в 100 раз для Cu(111). Полные тексты 1984 и 1985 гг. не прочитаны (закрыты подпиской).

**Чем закончилось.** Колебательный «перегрев» десорбирующих молекул по микрообратимости требует, чтобы колебательно возбуждённые H$_2$ прилипали лучше. Это потом подтвердили прямые пучковые опыты (Rettner, Michelsen, Auerbach и др., 1990-е; ссылки по Gross & Scheffler, arXiv:mtrl-th/9607006). Для H$_2$/Cu(111) получены пороги адсорбции: $E_0\approx0.6$ эВ ($v{=}0$), 0.3 эВ ($v{=}1$), 0.1 эВ ($v{=}2$) (*J. Chem. Phys.* 102, 4625, 1995). Десорбционные данные, обработанные через детальное равновесие, воспроизводят тренды адсорбционных измерений (*Surf. Sci.*, 1993). Всё это — по аннотациям, полные тексты не прочитаны.

**Почему это не довод против $W=1$:**
1. Вывод Kubiak 1985 сделан до измерения адсорбции по отдельным колебательным состояниям. Когда её измерили, детальное равновесие для H$_2$/Cu подтвердилось.
2. Даже в исходной форме речь идёт о *распределении* потока десорбции по состояниям, а не о равновесных константах. Неравновесное распределение десорбирующего потока **требуется** детальным равновесием, если адсорбция активирована по колебательной координате.
3. Опыт односторонний: десорбция в вакуум, без обратного потока и без равных давлений. Для $W$ нужны интегральные $K_1$, $K_2$.

**Что нужно для дословной сверки:**
* Kubiak, Sitz, Zare, JCP 83, 2538 (1985), doi:10.1063/1.449300;
* Kubiak, Sitz, Zare, JCP 81, 6397 (1984);
* Michelsen, Auerbach, JCP 94, 7502 (1991), «A critical examination of data on the dissociative adsorption and associative desorption of hydrogen at copper surfaces».
