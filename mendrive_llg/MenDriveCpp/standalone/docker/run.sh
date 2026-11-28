#!/usr/bin/env bash
# Использование:
#   ./run.sh build                 -- собрать образ
#   ./run.sh                       -- интерактивная оболочка в контейнере
#   ./run.sh make test             -- любая команда внутри контейнера
#   ./run.sh nohup ./mendrive -c examples/preisach2d_magnetic.cfg > out/p2d.log 2>&1
set -e
cd "$(dirname "$0")"
IMG=mendrive-standalone:latest
if [ "$1" = "build" ]; then docker build -t $IMG .; exit; fi
docker image inspect $IMG >/dev/null 2>&1 || docker build -t $IMG .
# монтируем MenDriveCpp целиком, чтобы работал ../mendrive_core.cpp
ROOT="$(cd .. && cd .. && pwd)"
docker run --rm -it -u "$(id -u):$(id -g)" -e HOME=/tmp \
  -v "$ROOT":/work/MenDriveCpp -w /work/MenDriveCpp/standalone \
  $IMG ${@:-bash}
