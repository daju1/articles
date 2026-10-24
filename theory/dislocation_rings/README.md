# dislocation_rings — версии расчёта

Модель частицы как дислокационного кольца (колец) в калибровочной теории
дислокаций Лазара–Брагинского. Поле строится из двух каналов:
`central` (кулоно-аналог) и `vortex` (аналог B, сила Пича–Кёлера).
Поправка Ферми — три режима, включая мгновенную нелокальную плёнку.

## Хронология версий

### v1–v2 (dislocation_rings.c/.h, dislocation_rings_bind.py)
- Заряд кольца СКАЛЯРНЫЙ: `total_charge` × `charge_sign`.
- Структура RingConfig в .h: ...,`total_charge`,`n_elem`,`charge_sign`.
- Биндинг `dislocation_rings_bind.py` СОГЛАСОВАН с .h (правильный).
- Три режима Ферми: PLAIN / FILM_RETARDED / FILM_INSTANT (мгновенная плёнка).

### v3 (dislocation_rings_v3.c, dislocation_rings_bind_v3.py)
- Та же скалярная физика (`total_charge`), улучшенный расчёт поля.
- ⚠️ БАГ: `dislocation_rings_bind_v3.py` рассогласован с C-кодом →
  **segmentation fault**. Две причины:
    1. структура биндинга содержит `burgers_angle`+`burgers_magnitude`
       (2 double) вместо `total_charge` (1 double) → сдвиг полей →
       мусорный `n_elem` → выход за массив;
    2. лишний аргумент `glide` (8 параметров вместо 7).
  Этот биндинг был ЗАДЕЛОМ под будущую v4 (векторный Бюргерс), ошибочно
  наложенным на скалярный C-код v3.
- ✅ ИСПРАВЛЕНО: `dislocation_rings_bind_v3_fixed.py` — структура приведена
  к .h (`total_charge`, 7 параметров, без glide). Рабочая v3.

### v4 (dislocation_rings_v4.c/.h, dislocation_rings_bind_v4.py)  ← НОВОЕ
- Заряд кольца ВЕКТОРНЫЙ — вектор Бюргерса b^a = ∫_S T^a:
  поля `burgers_magnitude` (|b|) и `burgers_angle` (θ_b).
- МЁБИУСОВО расслоение знака (флаг `mobius_enable`): знак ~cos(θ/2),
  период 4π (спинорная двузначность). Из раздела mobius-compatibility v3.tex.
- Строго обобщает v3: при `mobius_enable=0`, `burgers_angle=0` сводится к
  скалярной физике v3 (проверено Python-портом).
- РАЗДЕЛЕНИЕ ЧАСТИЦ одним флагом:
    * `mobius=0` + два кольца ±  → ЭЛЕКТРОН (диполь, период 2π);
    * `mobius=1` + одно кольцо   → МЁБИУС-ЧАСТИЦА (Уровень 1, НЕ электрон,
      период 4π). См. NOTE_mobius_particle.md.

## Сборка

    gcc -O3 -shared -fPIC -o dislocation_rings_v3.so dislocation_rings_v3.c -lm
    gcc -O3 -shared -fPIC -o dislocation_rings_v4.so dislocation_rings_v4.c -lm

Затем для v3 использовать `dislocation_rings_bind_v3_fixed.py`,
для v4 — `dislocation_rings_bind_v4.py`.

## Ключевые численные результаты (Python-порт, см. ANALYSIS_* и NOTE_*)

- Диполь-диполь: F·R² → плато (хвост 1/R²) + мультиполь 1/R⁴ на малых R.
- Мгновенный закон: F ∝ cosΔφ/R², K→4 при больших R/a, без осцилляций.
- Мёбиус-частица: поле убывает как 1/R (из-за разрыва cos(θ/2) в точке
  склейки θ=0~2π — "особая точка"), период поля 4π vs 2π у электрона.
- Тождество масштаба G: (l_Planck/a_compton)² = α_G (переформулировка,
  не вывод — см. insert_burgers_scale_v3.tex).


## Python-порты (python_ports/)
Чистые Python/numpy реализации C-кода для верификации логики без компиляции
(on-device песочница блокирует ctypes/компилятор):
- `port_v4_field.py` — ядро поля v4 (каналы, режимы Ферми, мёбиус).
- `port_v4_stability.py` — натяжение и равновесный радиус.
- `port_mobius_analysis.py` — асимптотика поля, мультиполи, кинк.
- `port_tests.py` — воспроизведение всех результатов.
- `REFERENCE_LOG.txt` — эталонный лог прогона тестов.
- `README.md` — описание портов.
Запуск: `cd python_ports && python3 port_tests.py`

## Обновление по массе мёбиуса (23 окт, вечер)
Из характерных длин Lazar (arXiv:0802.0670): ℓ₁²=a₁/(2μ), ядро r_core≃6ℓ₁,
масса m=b²√(μa₁/2)/c² — ФОРМА выведена. Но ЧИСЛО массы НЕ предсказывается:
система (μ,a₁) недоопределена одним кулоновским якорем (см.
../NOTE_mobius_mass_calibration.md и ../NOTE_mobius_mass_lazar_lengths.md).

## Полный список связанных файлов в ../
- ANALYSIS_coulomb_gravity_rings.md — кулон+гравитация, 4 проверки.
- NOTE_mobius_particle.md — мёбиусова частица как НЕ электрон.
- NOTE_mobius_mass_lazar_lengths.md — формула массы и длины Lazar.
- NOTE_mobius_mass_calibration.md — статус предсказания массы.
- OVERVIEW_session_2026-10-23.md — сквозной обзор сессии.
- insert_burgers_scale_v3.tex, insert_planck_sector_v3.tex — LaTeX-вставки.
- Drozdov_unified_model_article.tex — полная статья (с приложениями, библиографией).
