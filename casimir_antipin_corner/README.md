# Уголок Антипина: работа заново (с 20 ноября 2026)

Новая постановка задачи о силе Казимира на V-образный уголок (и соты): уравнения Максвелла целиком,
совместные ТМ+ТЕ, проницаемые материалы ($\varepsilon(\omega)$, $\mu(\omega)$), идеальный проводник — как предел.

## Соглашения
- Ось $z$ — вдоль образующей тел (как в `rezonatory_i_volnovody/MenDrive_qnm_theory.ipynb`); сечение — плоскость $(x,y)$.
- Гауссова система; поля $\propto e^{ik_z z - i\omega t}$; $\mathrm{rot}\,\mathbf E=\frac{i\omega}{c}\mu\mathbf H$, $\mathrm{rot}\,\mathbf H=-\frac{i\omega}{c}\varepsilon\mathbf E$.
- $s$ — длина дуги вдоль контура сечения, $\mathbf n$ — внешняя нормаль, $\boldsymbol\tau=\hat{\mathbf z}\times\mathbf n$.

## Файлы
- `Casimir_T_normalization_P1.ipynb` — нормировка $\langle\Theta_{ij}\rangle$ через функцию Грина и проверка П1 (две пластины) тремя способами: изображения, формула Лифшица, решатель на панелях.
- `Antipin_corner_problem_statement.ipynb` — **постановка задачи** (стиль MenDrive): утверждение Антипина и его модельные допущения; геометрия (варианты A и B, рис. `antipin_corner_geometry.png`); физическая модель; уравнения поля; определение силы через поток тензора; что следует до расчёта; проверки П1–П7; что сделано; план.
- `Casimir_cylinders_Maxwell_general.ipynb` — постановка: поперечные поля через $E_z,H_z$; граничные условия и связь ТМ–ТЕ;
  заряды и векторные токи; мнимая частота; интегральные уравнения; проверка решателя на точном решении для цилиндра.
  Тензор Максвелла и силы (§6–8) — ещё не сделаны.

## Предыдущая работа (не трогать, для справки)
- `casimir/Casimir_corner_Green_tensor_explained.ipynb` — идеальный проводник, скалярные ТМ/ТЕ, уголок с зазором и сплошной.
- `casimir/Antipin_angle_finite_DC.ipynb` — сверка с Дойчем–Канделасом, PFA, энергетический метод.
- `casimir/corner_lib.py` — функции из ноутбука выше одним файлом.
- `references/Casimir/Antipin/` — патент RU2610018C2 и статья Антипина 1702.0328v1.
- `casimir/ref/` — Deutsch & Candelas 1979, Fulling 2012, Milton 2013, Saharian 2008.

## Литература (`ref/`)
| файл | что | зачем |
|---|---|---|
| `Johnson2010_numerical_methods_Casimir_1007.0966.pdf` | S. G. Johnson, *Numerical methods for computing Casimir interactions*, Lect. Notes Phys. 834, ch. 6 (2011); arXiv:1007.0966 | обзор: функции Грина, мнимая частота, тензор напряжений, метод граничных элементов |
| `Reid2009_BEM_Casimir_3D_PRL_0904.0741.pdf` | Reid, Rodriguez, White, Johnson, PRL 103, 040401 (2009); arXiv:0904.0741 | граничные элементы для произвольных 3D-тел с острыми кромками |
| `Reid2011_BEM_Casimir_arbitrary_materials_1010.5539.pdf` | Reid, White, Johnson; arXiv:1010.5539 | то же для произвольных материалов: электрические **и магнитные** поверхностные токи |
| `Rahi2009_scattering_Casimir_PRD_0908.2649.pdf` | Rahi, Emig, Graham, Jaffe, Kardar, PRD 80, 085021 (2009); arXiv:0908.2649 | формула рассеяния (TGTG) для любых тел и материалов, вывод Лифшица |
| `Rodriguez2007_Casimir_arbitrary_geometries_PFA_failure_0704.1890.pdf` | Rodriguez, Ibanescu, Iannuzzi, Capasso, Joannopoulos, Johnson; arXiv:0704.1890 | тензор напряжений в сложной геометрии; где PFA не работает |

Не скачано (нет на arXiv / платное):

### Невзаимные (гиротропные, ферромагнитные) материалы — `ref/`
| файл | что | вывод для уголка |
|---|---|---|
| `GelbwaserKlimovsky2022_equilibrium_forces_nonreciprocal_2205.13470.pdf` | Gelbwaser-Klimovsky, Graham, Kardar, Krüger, PRB 106, 115106 (2022) | в равновесии (все T равны) сила и для невзаимных тел — градиент свободной энергии |
| `Milton2023_vacuum_torque_propulsion_nonreciprocal_2306.02197.pdf` | Milton и др., PRA 108, 022809 (2023) | вращающий момент и самодвижение — только вне равновесия |
| `ZhuFan2016_persistent_heat_current_nonreciprocal_1602.08454.pdf` | Zhu, Fan (2016) | равновесный незатухающий поток тепла в невзаимной системе многих тел |
| `Gangaraj2018_lateral_recoil_force_PTI_1711.04939.pdf` | Gangaraj, Hanson, Antezza, Silveirinha, PRB 97, 201108 (2018) | латеральная сила на *возбуждённый* атом у топологической поверхности |
| `Feigel2004_vacuum_momentum_dielectric_physics0304100.pdf` | Feigel, PRL 92, 020404 (2004), препринт | предсказан импульс вакуума у магнитоэлектрической жидкости в $\mathbf E\times\mathbf B$ |
| `RikkenVanTiggelen2011_magnetoelectric_momentum_experiment_1107.5284.pdf` | Rikken, van Tiggelen (2011) | эксперимент: предсказание о вакуумном вкладе опровергнуто |
| `NefedovRubi2017_lateral_drag_anisotropy_SciRep_1608.03200.pdf` | Nefedov, Rubi, Sci. Rep. 7, 6155 (2017); arXiv:1608.03200 | латеральная сила на частицу у наклонно-анизотропной поглощающей пластины; оспорена Gelbwaser-Klimovsky и др. 2022 |
| `KhandekarJacob2019_thermal_spin_photonics_NJP_1905.02745.pdf` | Khandekar, Jacob, New J. Phys. 21, 103030 (2019); arXiv:1905.02745 | равновесный незатухающий тепловой поток и спин у невзаимной пластины; предложено броуновское измерение латеральной силы; оспорено там же |

- K. A. Milton, *The Casimir Effect: Physical Manifestations of Zero-Point Energy* (World Scientific, 2001).
- A. W. Rodriguez, F. Capasso, S. G. Johnson, *The Casimir effect in microstructured geometries*, Nature Photonics 5, 211 (2011), doi:10.1038/nphoton.2011.39.
- Е. М. Лифшиц, Л. П. Питаевский, *Статистическая физика, ч. 2* (Ландау–Лифшиц, т. IX), главы о силах Ван-дер-Ваальса.
