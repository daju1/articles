#!/bin/sh
# Сборка разделяемой библиотеки модели каналов
gcc -O3 -shared -fPIC unified_channels.c -o unified_channels.so -lm
echo "собрано: unified_channels.so"
PYTHONIOENCODING=utf-8 python3 unified_channels_tests.py
