# python_binding — НАСТОЯЩИЙ Python-биндинг к C-коду

В отличие от `../python_implementation/` (где формулы переписаны на Python
С НУЛЯ), здесь Python ЛИНКУЕТСЯ к СКОМПИЛИРОВАННОМУ C-коду через ctypes.
То есть проверяется ИМЕННО ваш C, а не параллельная реализация.

## Как пользоваться
1. Собрать .so из C-файлов:
       ./build.sh
   (или вручную:
       gcc -O3 -shared -fPIC dislocation_rings_v4.c -o dislocation_rings_v4.so -lm
    в папке dislocation_rings/)
2. Прогнать тесты:
       python3 binding_tests.py

## Файлы
- `binding_v4.py` — ctypes-обёртка dislocation_rings_v4.so (векторный Бюргерс
  + мёбиус + устойчивость). Структура RingConfig ТОЧНО по dislocation_rings_v4.h:
  center,e1,e2,radius,omega,phase0, burgers_angle, burgers_magnitude,
  n_elem, charge_sign, mobius_enable. Функции: force_v4/force_multi (7 арг.,
  mode без glide), burgers_of, equilibrium_radius, total_energy.
- `binding_v3.py` — ctypes-обёртка для dislocation_rings_v3.so И базовой
  dislocation_rings.so (одинаковая структура из dislocation_rings.h:
  скалярный total_charge). Классы-экземпляры: B3.v3, B3.base.
- `build.sh` — сборка всех трёх .so.
- `binding_tests.py` — тесты через C.

## ВАЖНО: соответствие сигнатур (почему падали ноутбуки)
Ноутбуки v3/v4 звали:
  - force_v3(...)/force_v4(...) с 8 аргументами (лишний `glide`) — в C только 7;
  - force_multi(...)/make_dipole(...) — которых не было в старых биндингах.
Здесь это исправлено: force_v4/force_multi принимают 7 аргументов (mode без
glide), make_dipole и make_ring_b присутствуют. Если хотите glide-вариант —
его нет в текущем C API (compute_ring_dipole_force имеет 7 параметров).

## Три слоя проверки модели (итого)
1. C-код (эталон)              — dislocation_rings*.c/.h
2. python_binding/ (этот)      — линкует C через ctypes (проверка самого C)
3. python_implementation/      — независимая реализация формул на Python
Совпадение (2) и (3) = перекрёстная валидация C-кода и физики.


## Версия Python
Код совместим с Python 3.5+ (использован .format() вместо f-строк).
ВНИМАНИЕ: системный `python3` может быть старым (напр. 3.5.2). Если есть
новее — лучше им:
    python3.12 binding_tests.py      # или python3.11 / python3.10 / python3.9
Проверить, что установлено:
    ls /usr/local/bin/python3*  /opt/homebrew/bin/python3*  2>/dev/null
    which -a python3 python3.12 python3.11 python3.10 python3.9
numpy нужен только для python_implementation/, не для python_binding/.
