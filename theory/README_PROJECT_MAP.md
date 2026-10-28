# Единая геометродинамическая модель — карта проекта
### Папка theory/ · обновлено 24 октября 2026

Навигация по проекту: частицы как дефекты упругого эфира (дислокационные кольца
в T(3)-калибровке Лазара-Кадича-Эделена); заряд, масса, спин и гравитация — из
геометрии «плёнки Ферми» (граница интеграла действия). Файлы сгруппированы по
назначению. ⭐ = актуальное/эталонное; ⚠️ = содержит подтасовки (см. AUDIT).

---

## 1. С ЧЕГО НАЧАТЬ (главные документы)
- ⭐ **UNIFIED_THEORY_consistent_2026-10-24.md** — непротиворечивое ядро теории
  (упругий базис, 6 каналов, иерархия=m_e/m_P). Главный сводный документ.
- ⭐ **OPEN_PROBLEMS_2026-10-24.md** — реестр: доказано / принято / открыто /
  закрыто честно. Карта того, что сделано и что осталось.
- ⭐ **SESSION_SUMMARY_2026-10-24.md** — хронология всей работы (части 1-5).
- **OVERVIEW_session_2026-10-23.md** — обзор предыдущей сессии.

## 2. СТАТЬИ ДЛЯ ПУБЛИКАЦИИ (arXiv-пакет)
- ⭐ **Drozdov_unified_theory_intro_ru.tex/.pdf** — обзорная версия (RU).
- ⭐ **Drozdov_unified_theory_detailed_ru.tex/.pdf** — с полными доказательствами (RU).
- ⭐ **Drozdov_unified_theory_intro_en.tex** — английская обзорная (для arXiv).
- ⭐ **Drozdov_unified_theory_detailed_en.tex** — английская с доказательствами
  (нужна папка figures_unified/ с 3 .png).
- **arXiv_submission_materials.md** — title/abstract/категории/чеклист (2 варианта:
  cond-mat.other / physics.gen-ph).
- **endorsement_request_template.md** — шаблон письма-запроса endorsement (EN+RU).

## 3. ЗАМЕТКИ-ИССЛЕДОВАНИЯ (по темам, все честные, с делением доказано/открыто)
Ядро теории:
- NOTE_core_three_franck_shear_atom_2026-10-24.md — вектор Франка, J=2, E_atom=c²q_m.
- NOTE_Qij_field_equation_2026-10-24.md — тензорный J=2 канал (короткодействие 1/R⁵).
- NOTE_charge_gravity_unity_2026-10-24.md — единство заряд-гравитация (форма, не величина).
- NOTE_hierarchy_10e84_2026-10-24.md — иерархия 10⁸⁴ = проблема m_e/m_P.
- NOTE_so3_vortex_ether_revision_2026-10-24.md — SO(3)-сектор, дисклинация=вихрь.
- ⭐ FERMI_FILM_PARAMETERS_fixed_2026-10-24.md — эталон параметров плёнки (два масштаба).
Гравитация и 4/3:
- NOTE_gravity_from_Fermi_film_dislocation_dipoles.md — вывод гравитации через плёнку.
- NOTE_four_thirds_aether.md — парадокс 4/3 (k_edge анизотропна; с самокоррекциями).
- ANALYSIS_coulomb_gravity_rings.md — ранний анализ кулон/гравитация.
Мёбиус/масса/заряд:
- NOTE_mobius_particle.md, NOTE_mobius_mass_calibration.md,
  NOTE_mobius_mass_lazar_lengths.md — мёбиус-частица и её масса.
- NOTE_v3_glide_magnetic_charge_dimensional.md — размерность магнитного заряда.
- NOTE_v3_calibration_electron.md — калибровка по электрону.
- NOTE_twist_plane_decoupled_from_mu.md — развязка плоскости твиста и μ.
- NOTE_weak_sector_and_dirac_g2.md — слабый сектор (хемитропные константы) и g=2.
Сонолюминесценция (в Schwinger/, не здесь):
- см. ../Schwinger/NOTE_sonoluminescence_honest_2026-10-24.md
- NOTE_kuzmenko_analysis_2026-10-24.md — анализ работ Кузьменко (E=ω²R³).

## 4. АУДИТ (ОБЯЗАТЕЛЬНО прочитать перед правкой старых .tex)
- ⭐ **AUDIT_fabricated_vs_honest_2026-10-24.md** — каталог: что честное, что
  сфальсифицировано восторженным ИИ (завышения на 51-135 порядков, подмены τ).
- ⚠️ **searl.tex/.pdf/.ipynb, searlNotes.bib** — ПОДТАСОВКА (антигравитация, over-unity).
- ⚠️ Старые .tex ниже могут содержать «closure»-тавтологии — сверять с AUDIT.

## 5. СТАРЫЕ ВЕРСИИ СТАТЕЙ (историч., частично с подтасовками)
- Drozdov_theory_full_exposition_v2/v3.tex/.pdf
- Drozdov_three_variants_calculation.tex/.pdf
- Drozdov_unified_model_article.tex/.pdf
- Lazar_Braginsky_Drozdov_full_theory.tex/.pdf
- insert_burgers_scale_v3.tex, insert_planck_sector_v3.tex — вставки для v3.
- classic_spin_model.pdf

## 6. КОД (C + ctypes-биндинги + независимые Python-реализации)
- dislocation_rings/ — основной код (v4.c/.h, биндинги, python_implementation/).
- dislocation_rings_v3_vector/ — v3 с векторными зарядами Бюргерса.
- disclination_fermion_v6/, disclination_fermion_v7/ — диспирационный фермион.
- mobius_dipole_v3B/ — пара мёбиус-колец (ветка B).
- positronium_v5/ — позитроний с нуля.
- unified_model_channels/ — все 6 каналов × 4 типа частиц.
  Сборка в каждой папке: gcc -O3 -shared -fPIC X.c -o X.so -lm; затем
  PYTHONIOENCODING=utf-8 python3 X_tests.py

## 7. РИСУНКИ И АНИМАЦИИ
- figures_unified/ — channel_scaling, stability, twist_period (для статей).
- disclination_fermion.png, electron_positron_glide_loops.gif,
  fermi_cap_circular_motion.gif, fixed_rings_precessing_moment.gif,
  fermi_local_perpendicular_thick_tube.png, film_3d_time_embedding.png,
  film_time_field.png, force_near_far_zone.png.

## 8. ИСТОЧНИКИ И МАТЕРИАЛЫ
- ref/ — Article1-3 Кузьменко (E=ω²R³), прочитаны полностью.
- ../references/gauge_theory_defects/ — Лазар, Осипов, Hehl и др. (вне theory/).
- планковская жесткость вакуум-3.txt — длинный чат про планковскую жёсткость (Gemini).
- cui2019.pdf.txt — справочный текст.
- chat-*.txt/.md/.pdf — архивные логи прежних ИИ-диалогов (сырьё, не результат).

## 9. СБОРКА PDF
- tex.sh — скрипт сборки (xelatex→bibtex→xelatex→pdflatex). Запуск:
  cd ~/articles/theory && sh tex.sh
- Для en-версий: pdflatex дважды каждую (detailed требует figures_unified/).

---
## Статус проекта (на 24 окт 2026)
ФОРМА выводится (законы 1/R², статистика, устойчивость, вектор Франка, единство
формы заряд-гравитация). ВЕЛИЧИНА (G, α, иерархия) упирается в ОДНУ стену —
малость m_e/m_P (стандартная проблема иерархии). Это честная, единственная
граница теории. arXiv-пакет готов; осталось внешнее (компиляция, endorsement).
