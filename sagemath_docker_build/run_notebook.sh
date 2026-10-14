#!/bin/bash
set -e
set -x

CONTAINER_NAME="sage_build_container"
# PROJECT_ROOT=$(dirname ${PWD})
CONTAINER_JUPYTER_BIN="${PWD}/sage/local/bin/jupyter"
HOST_MOUNT_SRC="${PWD}"
CONTAINER_MOUNT_DST="${PWD}/work"
ARTICLES_ROOT="$(dirname "$(realpath "$0")")"
SEARCH_ROOTS=(${PWD})
INPUT="$1"

if [ -z "$INPUT" ]; then
  echo "Usage: $0 <notebook.ipynb>"
  exit 1
fi

# --- Определяем, где выполняется скрипт: внутри контейнера или на хосте ---
if [ -f /.dockerenv ]; then
  RUN_MODE="local"
else
  RUN_MODE="docker"
fi

# --- Разрешение пути к ноутбуку (как раньше) ---
if [ -f "$INPUT" ]; then
  NOTEBOOK="$(realpath "$INPUT")"
elif [ -f "$ARTICLES_ROOT/$INPUT" ]; then
  NOTEBOOK="$(realpath "$ARTICLES_ROOT/$INPUT")"
  echo "Найден относительно $ARTICLES_ROOT"
else
  BASENAME_ONLY="$(basename "$INPUT")"
  echo "Файл не найден по прямому пути, ищу по имени: $BASENAME_ONLY ..."
  MATCHES=$(find "${SEARCH_ROOTS[@]}" -name "$BASENAME_ONLY" -type f 2>/dev/null)

  if [ -z "$MATCHES" ]; then
    echo "Ошибка: файл '$BASENAME_ONLY' не найден ни в одном из известных корней:"
    printf '  %s\n' "${SEARCH_ROOTS[@]}"
    exit 1
  fi

  MATCH_COUNT=$(echo "$MATCHES" | wc -l)
  if [ "$MATCH_COUNT" -gt 1 ]; then
    echo "Найдено несколько файлов с таким именем, укажите нужный путь явно:"
    echo "$MATCHES"
    exit 1
  fi

  NOTEBOOK="$MATCHES"
  echo "Найден единственный вариант: $NOTEBOOK"
fi

NOTEBOOK_DIR="$(dirname "$NOTEBOOK")"
NOTEBOOK_NAME="$(basename "$NOTEBOOK")"
BASENAME="${NOTEBOOK_NAME%.ipynb}"
RUNNING="${BASENAME}_running.ipynb"
LOG="${BASENAME}_run.log"

# ============================================================
# ВАРИАНТ А — выполняется внутри контейнера напрямую
# ============================================================
if [ "$RUN_MODE" = "local" ]; then
  echo "Режим: локальный запуск внутри контейнера"

  if [ ! -x "$CONTAINER_JUPYTER_BIN" ]; then
    echo "Ошибка: jupyter не найден или не исполняемый: $CONTAINER_JUPYTER_BIN"
    exit 1
  fi

  cd "$NOTEBOOK_DIR"

  if pgrep -f "nbconvert.*${RUNNING}" > /dev/null; then
    echo "Ошибка: расчёт для $RUNNING уже запущен"
    exit 1
  fi

  # export SAGE_ROOT=${HOST_MOUNT_SRC}/sage

  # cp "$NOTEBOOK_NAME" "$RUNNING"

  # nohup "$CONTAINER_JUPYTER_BIN" nbconvert --to notebook --execute \
  #   --ExecutePreprocessor.timeout=-1 \
  #   --output "$RUNNING" \
  #   "$RUNNING" > "$LOG" 2>&1 &

  # Запускаем Python-скрипт в фоне через nohup
  nohup ${HOST_MOUNT_SRC}/sage${SAGE_VER}/local/bin/python3 "${HOST_MOUNT_SRC}/nb_run.py" "$NOTEBOOK" "$RUNNING" > "$LOG" 2>&1 &

  PID=$!
  disown

  echo "Запущено. PID: $PID"
  echo "Лог: $NOTEBOOK_DIR/$LOG"

  sleep 2
  if ! kill -0 "$PID" 2>/dev/null; then
    echo "Внимание: процесс завершился почти сразу — проверьте $LOG"
    tail -20 "$LOG"
    exit 1
  fi

  echo "Процесс жив, расчёт идёт."
  echo "После успешного завершения выполните: mv $NOTEBOOK_DIR/$RUNNING $NOTEBOOK"

# ============================================================
# ВАРИАНТ Б — выполняется на хосте, запускает расчёт через docker exec -d
# ============================================================
else
  echo "Режим: запуск на хосте через docker exec в $CONTAINER_NAME"

  # Проверяем, что файл вообще виден контейнеру (лежит под смонтированным путём)
  case "$NOTEBOOK" in
    "$HOST_MOUNT_SRC"/*)
      CONTAINER_NOTEBOOK_DIR="${NOTEBOOK_DIR/#$HOST_MOUNT_SRC/$CONTAINER_MOUNT_DST}"
      ;;
    *)
      echo "Ошибка: файл лежит вне смонтированной в контейнер директории ($HOST_MOUNT_SRC)"
      echo "Контейнер его не увидит. Переместите файл в $HOST_MOUNT_SRC или используйте Вариант А внутри контейнера."
      exit 1
      ;;
  esac

  if ! docker ps --format '{{.Names}}' | grep -qx "$CONTAINER_NAME"; then
    echo "Ошибка: контейнер $CONTAINER_NAME не запущен"
    exit 1
  fi

  if docker exec "$CONTAINER_NAME" pgrep -f "nbconvert.*${RUNNING}" > /dev/null 2>&1; then
    echo "Ошибка: расчёт для $RUNNING уже запущен в контейнере"
    exit 1
  fi

  docker exec "$CONTAINER_NAME" bash -c "cd '$CONTAINER_NOTEBOOK_DIR' && cp '$NOTEBOOK_NAME' '$RUNNING'"

  docker exec -d "$CONTAINER_NAME" bash -c \
    "cd '$CONTAINER_NOTEBOOK_DIR' && '$CONTAINER_JUPYTER_BIN' nbconvert --to notebook --execute --ExecutePreprocessor.timeout=-1 --output '$RUNNING' '$RUNNING' > '$LOG' 2>&1"

  echo "Запущено внутри контейнера (detached через docker exec -d)."
  echo "Лог (виден и с хоста): $NOTEBOOK_DIR/$LOG"

  sleep 2
  if docker exec "$CONTAINER_NAME" pgrep -f "nbconvert.*${RUNNING}" > /dev/null 2>&1; then
    echo "Процесс жив, расчёт идёт."
  else
    echo "Внимание: процесс не найден через 2 секунды — проверьте лог:"
    tail -20 "$NOTEBOOK_DIR/$LOG" 2>/dev/null || echo "(лог ещё не появился)"
    exit 1
  fi

  echo "После успешного завершения выполните: mv $NOTEBOOK_DIR/$RUNNING $NOTEBOOK"
fi