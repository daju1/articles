# Материалы для загрузки на arXiv
### Статьи: Drozdov_unified_theory_intro_en.tex + Drozdov_unified_theory_detailed_en.tex (24 окт 2026)

Справочный комплект для submission. Всё, что просит форма arXiv, собрано ниже.
Готовы ДВЕ английские версии: intro (обзорная, ~11 стр) и detailed (с полными
доказательствами, ~16 стр + 3 рисунка). Рекомендация: грузить как отдельные
препринты со взаимными ссылками, либо detailed как основную, intro — краткую.

---

## 1. TITLE (поле Title)
```
A Gauge Theory of Line Defects in an Elastic Continuum: Dislocation Rings,
the Peach-Koehler Force, and Emergent Interaction Laws
```
(Заголовок сдвинут под cond-mat акцент — калибровочная теория дефектов, а не
«unified model». Подзаголовок «A Geometrodynamic Model of Particles as Elastic
Defects» стоит в самих .tex как \large-строка. Для Варианта B (gen-ph) можно
вернуть прежний «unified» заголовок.)

## 2. AUTHORS (поле Authors)
```
A. Yu. Drozdov
```
(добавьте аффилиацию и ORCID, если есть — arXiv просит отдельным полем)

## 3. ABSTRACT (поле Abstract, ≤1920 знаков — у нас ~1180, в лимите)
```
We present a consistent account of a geometrodynamic model in which elementary
particles are treated as topological defects (dislocation and disclination
rings) of an elastic aether, and the fundamental interactions appear as
different projections of a single elastic response of the medium. The model is
built on the Kadic-Edelen and Lazar gauge theory of dislocations, in which the
magnetic field is identified with the dislocation density and the Peach-Koehler
force plays the role of the Lorentz force. Gravity is introduced through a
"Fermi film" mechanism: a deformable boundary of the action integral held
orthogonal to the curved world tubes of accelerating charges, generalizing
Fermi's 1923 treatment of the 4/3 paradox. We show that the inverse-square law
for gravity is a theorem (a consequence of the Laplacian Green function), not a
fit. The model is formulated in an elastic unit basis (L,T,M), without invoking
electric charge as a primary entity. The limits are drawn honestly: the
structure of the interactions (power laws, signs, hierarchy) is derived,
together with relations between scales, whereas the absolute values of the
fundamental constants (alpha, G) remain input parameters. The hierarchy of two
stiffness scales is shown to be equivalent to the standard hierarchy problem
m_e << m_Planck.
```
(1180 знаков с пробелами — в пределах лимита 1920. ASCII-safe: греческие буквы
заменены словами, тире обычное — так arXiv не ругается на поле ввода.)

## 4. CATEGORIES (поле Subject class)
ДВА ВАРИАНТА ПОЗИЦИОНИРОВАНИЯ (выбрать до загрузки):

**Вариант A (рекомендуется) — cond-mat.other как основная.**
Позиционирование: не «теория всего», а КАЛИБРОВОЧНАЯ ТЕОРИЯ ДЕФЕКТОВ упругого
континуума. Дислокации, вектор Бюргерса, сила Пича-Кёлера, теория Лазара —
стандартный аппарат cond-mat. Преимущества: (1) endorsement-путь обычно проще,
т.к. тема законная для cond-mat; (2) рецензент-профильный. Cross-list:
physics.gen-ph (вторичная). (gr-qc не рекомендуется: плёнка Ферми — НЕ GHY-член,
она не входит в лагранжиан.)
Для этого варианта — использовать ASCII-abstract из §3a (переформулирован).

**Вариант B — physics.gen-ph как основная.**
Позиционирование: оригинальная объединительная модель. Честно, но gen-ph
требует endorsement (см. §7) и воспринимается как «fringe». Cross-list:
cond-mat.other. Abstract из §3 (общая формулировка).

РЕКОМЕНДАЦИЯ: Вариант A (cond-mat.other). Он точнее отражает реальное
содержание (калибровочная теория дефектов), а не претензию на «всё», и снижает
барьер endorsement. Претензии на объединение — внутри текста, честно
разделённые на proved/assumed/open, а не в позиционировании категории.

## 3a. ABSTRACT для Варианта A (cond-mat.other, акцент на теорию дефектов)
```
We develop a gauge-theoretic model in which elementary excitations are treated
as topological line defects (dislocation and disclination rings) of an elastic
continuum, and their interactions arise as projections of a single elastic
response. Building on the Kadic-Edelen and Lazar translational (T(3)) gauge
theory of dislocations, we identify the magnetic field with the dislocation
density and the Peach-Koehler force with the Lorentz force, and formulate the
whole construction in a purely elastic unit basis (L,T,M) in which the magnetic
charge is the Burgers vector and has the dimension of length. We prove several
structural results: the dimensional closure of mass; a stability theorem for a
precessing defect pair (angular-momentum barrier); the Frank vector Omega = pi
of a cos(phi/2) Burgers twist (a wedge half-disclination, spin 1/2); and an
inverse-square interaction law as a Green-function theorem in three dimensions.
The approach places a long-interaction (gravity-like) channel on a Fermi-film
boundary term generalizing Fermi's 1923 treatment of the 4/3 problem. We are
explicit about limits: power laws, signs and scale relations are derived,
whereas the absolute coupling constants remain inputs, and the required ratio
of two elastic stiffness scales is shown to be equivalent to the m_e << m_P
hierarchy.
```
(1430 знаков — в пределах 1920. ASCII-safe. Акцент: теория дефектов, а не
«теория всего»; слова gravity/unification приглушены до «long-interaction
channel» и «elastic response».)

## 5. COMMENTS (поле Comments, опционально но полезно)
```
11 pages, 2 tables. Working preprint; introductory exposition. Companion
detailed version with full proofs in preparation.
```

## 6. ФАЙЛЫ ДЛЯ ЗАГРУЗКИ
ДВА препринта (взаимные ссылки в тексте уже есть — «introductory/detailed»):

**intro-препринт:**
- `Drozdov_unified_theory_intro_en.tex` (самодостаточный, библиография внутри).
- Рисунков НЕТ. Грузить один .tex.

**detailed-препринт:**
- `Drozdov_unified_theory_detailed_en.tex` (самодостаточный, библиография внутри).
- Рисунки ЕСТЬ (3 шт.): нужна папка `figures_unified/` с файлами
  `channel_scaling.png`, `stability.png`, `twist_period.png`. Грузить .tex +
  три .png (arXiv принимает архив; пути в \includegraphics уже относительные
  `figures_unified/...`). ПРОВЕРЬТЕ, что все три .png на месте перед сборкой.

Общее:
- Библиография внутри через \begin{thebibliography} — отдельный .bbl НЕ нужен.
- НЕ грузить .pdf — arXiv собирает сам из исходника.
- Проверьте локально pdflatex ДВАЖДЫ (для \tableofcontents и ссылок), для
  каждого файла:
    pdflatex Drozdov_unified_theory_intro_en.tex   (×2)
    pdflatex Drozdov_unified_theory_detailed_en.tex (×2)
  detailed требует наличия figures_unified/*.png, иначе \includegraphics упадёт.

## 7. ВАЖНО про endorsement (physics.gen-ph)
Первая загрузка в physics.gen-ph требует endorsement от автора, уже
публиковавшегося в этой категории. Без него arXiv не примет. Варианты:
(а) попросить поручительство у коллеги с arXiv-историей в physics.gen-ph;
(б) если есть co-author с публикациями — грузить от него;
(в) рассмотреть cond-mat.other (теория дефектов) — там endorsement может
    пройти проще, т.к. дислокации — законная тема cond-mat.
Кузьменко (ваш коллега) грузил на ResearchGate, не arXiv — у него endorsement
может не быть. Уточните, кто в вашем окружении имеет arXiv-историю.

## 8. ЛИЦЕНЗИЯ
arXiv попросит выбрать лицензию. Рекомендация для препринта: arXiv perpetual
non-exclusive license (стандарт, не мешает дальнейшей журнальной публикации).

## 9. ЧЕК-ЛИСТ ПЕРЕД SUBMIT
- [ ] pdflatex собирается без ошибок (прогнать дважды)
- [ ] abstract вставлен в поле (ASCII-версия из §3)
- [ ] категория выбрана (physics.gen-ph + endorsement решён)
- [ ] заголовок/авторы совпадают с .tex
- [ ] проверено, что все \cite разрешаются (нет [?] в PDF)
- [ ] метаданные (ORCID, аффилиация) добавлены

---
## Открытый вопрос к автору
arXiv для physics.gen-ph требует endorsement. Есть ли у вас или соавтора
arXiv-история? Если нет — это главное препятствие, решаем до технической
загрузки. Альтернатива — cond-mat.other (дислокации как законная тема).
