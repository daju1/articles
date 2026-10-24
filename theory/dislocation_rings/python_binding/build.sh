#!/usr/bin/env bash
# build.sh — компиляция всех C-файлов в .so для ctypes-биндингов.
# Запускать из папки dislocation_rings/ (на уровень выше python_binding/),
# ИЛИ из python_binding/ — скрипт сам перейдёт в родительскую папку с .c.
set -e
cd "$(dirname "$0")/.."   # -> dislocation_rings/
CC="${CC:-gcc}"
FLAGS="-O3 -shared -fPIC"
echo "Компилирую в $(pwd) компилятором $CC ..."
$CC $FLAGS dislocation_rings.c    -o dislocation_rings.so    -lm && echo "  OK dislocation_rings.so"
$CC $FLAGS dislocation_rings_v3.c -o dislocation_rings_v3.so -lm && echo "  OK dislocation_rings_v3.so"
$CC $FLAGS dislocation_rings_v4.c -o dislocation_rings_v4.so -lm && echo "  OK dislocation_rings_v4.so"
echo "Готово. Теперь: cd python_binding && python3 binding_tests.py"
