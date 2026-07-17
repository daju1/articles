#!/bin/bash
# mendrive_build.sh — Удобный скрипт сборки библиотеки MenDrive
# Поддерживает выбор точности, уровня логгирования и режимов работы

set -e  # Автоматический выход при ошибке

# Цвета для вывода
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
BLUE='\033[0;34m'
NC='\033[0m' # No Color

# Параметры по умолчанию
PRECISION="long_double"
LOGLEVEL="info"
OPTIMIZE="3"
DEBUG=0
CLEAN=0
TARGET="mendrive"
GSL_PATH="../vimaniks/gsl/local"
KY_MODE=1
QNM_MODE=1
GYRO_TENSOR=0
JOBS=$(nproc 2>/dev/null || echo 2)
STATIC=0

# Справка
show_help() {
    cat << EOF
${BLUE}MenDrive Build Script v1.0${NC}

${GREEN}Использование:${NC} $0 [опции]

${GREEN}Основные опции:${NC}
  ${YELLOW}-p, --precision <режим>${NC}   Точность арифметики:
                                ${GREEN}long_double${NC}  — 80 бит (по умолчанию)
                                ${GREEN}float128${NC}     — 128 бит (__float128, требуется libquadmath)
                                ${GREEN}mpfr_512${NC}     — 512 бит (требует libmpfr + libgmp)

  ${YELLOW}-l, --loglevel <уровень>${NC}  Уровень логгирования:
                                ${GREEN}off${NC}    — отключено
                                ${GREEN}error${NC}  — только ошибки
                                ${GREEN}warn${NC}   — предупреждения + ошибки
                                ${GREEN}info${NC}   — информационные сообщения (по умолчанию)
                                ${GREEN}debug${NC}  — отладочная информация
                                ${GREEN}trace${NC}  — подробное трассирование

${GREEN}Дополнительные опции:${NC}
  ${YELLOW}-O, --optimize <0-3>${NC}      Уровень оптимизации GCC (по умолчанию: 3)
  ${YELLOW}-g, --debug${NC}              Включить отладочную информацию (-g)
  ${YELLOW}-c, --clean${NC}              Очистить объектные файлы перед сборкой
  ${YELLOW}-j, --jobs <N>${NC}           Параллельная компиляция (по умолчанию: $JOBS ядер)
  ${YELLOW}-s, --static${NC}             Собрать статическую библиотеку (.a) вместо динамической (.so)
  ${YELLOW}-t, --target <имя>${NC}       Имя целевой библиотеки (по умолчанию: "mendrive")
  ${YELLOW}-G, --gsl-path <путь>${NC}    Путь к GSL (по умолчанию: "../vimaniks/gsl/local")
  ${YELLOW}-K, --no-ky${NC}              Отключить режим KY (по умолчанию включён)
  ${YELLOW}-Q, --qnm${NC}                Включить режим QNM (по умолчанию выключен)
  ${YELLOW}-T, --gyro-tensor${NC}        Включить режим гиротензора (по умолчанию выключен)
  ${YELLOW}-h, --help${NC}               Показать эту справку

${GREEN}Примеры использования:${NC}
  $0 -p float128 -l debug
      Сборка с 128-битной точностью и подробным логгированием

  $0 -p mpfr_512 -O 2 -g -j 4
      Сборка с 512-битной точностью, оптимизация -O2, отладка, 4 потока

  $0 -c && $0 -p long_double -l info
      Полная пересборка с параметрами по умолчанию

  $0 -s -p float128
      Сборка статической библиотеки с 128-битной точностью
EOF
    exit 0
}

# Парсинг аргументов
while [[ $# -gt 0 ]]; do
    case $1 in
        -p|--precision)     PRECISION="$2"; shift 2 ;;
        -l|--loglevel)      LOGLEVEL="$2"; shift 2 ;;
        -O|--optimize)      OPTIMIZE="$2"; shift 2 ;;
        -g|--debug)         DEBUG=1; shift ;;
        -c|--clean)         CLEAN=1; shift ;;
        -j|--jobs)          JOBS="$2"; shift 2 ;;
        -s|--static)        STATIC=1; shift ;;
        -t|--target)        TARGET="$2"; shift 2 ;;
        -G|--gsl-path)      GSL_PATH="$2"; shift 2 ;;
        -K|--no-ky)         KY_MODE=0; shift ;;
        -Q|--qnm)           QNM_MODE=1; shift ;;
        -T|--gyro-tensor)   GYRO_TENSOR=1; shift ;;
        -h|--help)          show_help ;;
        *)                  echo -e "${RED}Ошибка:${NC} Неизвестный параметр: $1"; show_help ;;
    esac
done

# Валидация параметров
case "$PRECISION" in
    long_double|float128|mpfr_512) ;;
    *)
        echo -e "${RED}Ошибка:${NC} Неизвестный режим точности '$PRECISION'"
        echo "Допустимые значения: long_double, float128, mpfr_512"
        exit 1
        ;;
esac

case "$LOGLEVEL" in
    off|error|warn|info|debug|trace) ;;
    *)
        echo -e "${RED}Ошибка:${NC} Неизвестный уровень логгирования '$LOGLEVEL'"
        echo "Допустимые значения: off, error, warn, info, debug, trace"
        exit 1
        ;;
esac

if ! [[ "$OPTIMIZE" =~ ^[0-3]$ ]]; then
    echo -e "${RED}Ошибка:${NC} Уровень оптимизации должен быть 0-3"
    exit 1
fi

if ! [[ "$JOBS" =~ ^[0-9]+$ ]] || [ "$JOBS" -lt 1 ]; then
    echo -e "${RED}Ошибка:${NC} Некорректное число потоков: $JOBS"
    exit 1
fi

# Проверка зависимостей
check_dependency() {
    local dep=$1
    local pkg=$2
    if ! command -v "$dep" &> /dev/null; then
        echo -e "${YELLOW}⚠️  Предупреждение:${NC} $pkg не найден в PATH"
        return 1
    fi
    return 0
}

echo -e "${BLUE}🔍 Проверка зависимостей...${NC}"
check_dependency gcc "GCC compiler" || exit 1
if [ "$PRECISION" = "float128" ]; then
    echo -e "   Проверка libquadmath... \c"
    if ! gcc -lquadmath -x c -o /dev/null - <<<'' 2>/dev/null; then
        echo -e "${RED}НЕ НАЙДЕНА${NC}"
        echo -e "${RED}Ошибка:${NC} Для режима float128 требуется библиотека libquadmath"
        echo "Установите: sudo apt install libquadmath0 libquadmath-dev (Debian/Ubuntu)"
        exit 1
    fi
    echo -e "${GREEN}OK${NC}"
fi
if [ "$PRECISION" = "mpfr_512" ]; then
    echo -e "   Проверка libmpfr/libgmp... \c"
    if ! gcc -lmpfr -lgmp -x c -o /dev/null - <<<'' 2>/dev/null; then
        echo -e "${RED}НЕ НАЙДЕНЫ${NC}"
        echo -e "${RED}Ошибка:${NC} Для режима mpfr_512 требуются библиотеки libmpfr и libgmp"
        echo "Установите: sudo apt install libmpfr-dev libgmp-dev (Debian/Ubuntu)"
        exit 1
    fi
    echo -e "${GREEN}OK${NC}"
fi

# Очистка
if [ $CLEAN -eq 1 ]; then
    echo -e "${BLUE}🧹 Очистка объектных файлов...${NC}"
    rm -f *.o *.so lib${TARGET}.a ${TARGET}.a 2>/dev/null || true
fi

# Настройка флагов
CFLAGS="-O${OPTIMIZE} -fPIC -Wall -Wextra -Wno-unused-parameter"
LDFLAGS="-shared -fPIC"
ARITH_DEFINE=""
LOG_DEFINE=""
LIBS="-lm"

# Точность
case "$PRECISION" in
    float128)
        ARITH_DEFINE="-DARITH_FLOAT128"
        LIBS="$LIBS -lquadmath"
        PREC_LABEL="__float128 (128 бит, ~34 знака)"
        ;;
    mpfr_512)
        ARITH_DEFINE="-DARITH_MPFR_512 -DMPFR_PREC=512"
        LIBS="$LIBS -lmpfr -lgmp"
        PREC_LABEL="MPFR 512-bit (~154 знака)"
        ;;
    *)
        ARITH_DEFINE="-DARITH_LONG_DOUBLE"
        PREC_LABEL="long double (80 бит, ~19 знаков)"
        ;;
esac

# Логгирование
case "$LOGLEVEL" in
    off)   LOG_DEFINE="-DMPREC_LOG_LEVEL=0" ;;
    error) LOG_DEFINE="-DMPREC_LOG_LEVEL=1" ;;
    warn)  LOG_DEFINE="-DMPREC_LOG_LEVEL=2" ;;
    info)  LOG_DEFINE="-DMPREC_LOG_LEVEL=3" ;;
    debug) LOG_DEFINE="-DMPREC_LOG_LEVEL=4" ;;
    trace) LOG_DEFINE="-DMPREC_LOG_LEVEL=5" ;;
esac

# Отладка
if [ $DEBUG -eq 1 ]; then
    CFLAGS="$CFLAGS -g -DDEBUG"
    DEBUG_LABEL="ВКЛЮЧЕНА (-g)"
else
    CFLAGS="$CFLAGS -DNDEBUG"
    DEBUG_LABEL="отключена"
fi

# Режимы
MODE_FLAGS=""
[ $KY_MODE -eq 1 ] && MODE_FLAGS="$MODE_FLAGS -DKY"
[ $QNM_MODE -eq 1 ] && MODE_FLAGS="$MODE_FLAGS -DQNM"
[ $GYRO_TENSOR -eq 1 ] && MODE_FLAGS="$MODE_FLAGS -DGYRO_TENSOR"

# GSL
if [ -d "$GSL_PATH" ]; then
    CFLAGS="$CFLAGS -I$GSL_PATH/include"
    LDFLAGS="$LDFLAGS -L$GSL_PATH/lib -Wl,-rpath,'\$ORIGIN/../vimaniks/gsl/local/lib'"
    LIBS="$LIBS -lgsl -lgslcblas"
    GSL_STATUS="${GREEN}найден${NC} ($GSL_PATH)"
else
    GSL_STATUS="${YELLOW}не найден${NC} (некоторые функции могут быть недоступны)"
fi

# Финальные флаги
CFLAGS="$CFLAGS $ARITH_DEFINE $LOG_DEFINE $MODE_FLAGS"
[ $STATIC -eq 1 ] && LDFLAGS="" && TARGET_LIB="lib${TARGET}.a" || TARGET_LIB="${TARGET}.so"

# Вывод конфигурации
echo -e "${BLUE}╔════════════════════════════════════════════════════════════════════╗${NC}"
echo -e "${BLUE}║${NC} ${GREEN}Конфигурация сборки MenDrive${NC}                                     ${BLUE}║${NC}"
echo -e "${BLUE}╠════════════════════════════════════════════════════════════════════╣${NC}"
echo -e "${BLUE}║${NC} Точность:       ${GREEN}$PREC_LABEL${NC}                           ${BLUE}║${NC}"
echo -e "${BLUE}║${NC} Логгирование:   ${GREEN}$LOGLEVEL${NC}                                      ${BLUE}║${NC}"
echo -e "${BLUE}║${NC} Оптимизация:    ${GREEN}-O$OPTIMIZE${NC}                                    ${BLUE}║${NC}"
echo -e "${BLUE}║${NC} Отладка:        ${GREEN}$DEBUG_LABEL${NC}                              ${BLUE}║${NC}"
echo -e "${BLUE}║${NC} GSL:            $GSL_STATUS                              ${BLUE}║${NC}"
echo -e "${BLUE}║${NC} Режимы:         ${GREEN}$( [ $KY_MODE -eq 1 ] && echo "KY " )$( [ $QNM_MODE -eq 1 ] && echo "QNM " )$( [ $GYRO_TENSOR -eq 1 ] && echo "GYRO_TENSOR" )${NC}                 ${BLUE}║${NC}"
echo -e "${BLUE}║${NC} Потоков:        ${GREEN}$JOBS${NC}                                          ${BLUE}║${NC}"
echo -e "${BLUE}║${NC} Цель:           ${GREEN}$TARGET_LIB${NC}                                 ${BLUE}║${NC}"
echo -e "${BLUE}╚════════════════════════════════════════════════════════════════════╝${NC}"
echo ""

# Список исходников
SOURCES=(
    mendrive_det_qnm.c
    mendrive_isolines.c
    mendrive_isolines_traced.c
    mendrive_contour_intersections.c
    mendrive_contour_intersections_sharp.c
    mendrive_characteristic_roots.c
    mendrive_root.c
    mendrive_newton.c
    mendrive_newton_classic/mendrive_newton_classic.c
    mendrive_linsolve.c
)

# Проверка наличия файлов
MISSING=()
for src in "${SOURCES[@]}"; do
    [ ! -f "$src" ] && MISSING+=("$src")
done
if [ ${#MISSING[@]} -gt 0 ]; then
    echo -e "${RED}Ошибка:${NC} Отсутствуют исходные файлы:"
    for f in "${MISSING[@]}"; do echo "  - $f"; done
    exit 1
fi

# Компиляция объектных файлов
echo -e "${BLUE}🔨 Компиляция исходных файлов ($JOBS потоков)...${NC}"
COMPILE_JOBS=()
for src in "${SOURCES[@]}"; do
    obj="${src%.c}.o"
    echo -e "   ${GREEN}→${NC} $src"
    gcc $CFLAGS -c "$src" -o "$obj" &
    COMPILE_JOBS+=($!)
done

# Ожидание завершения всех задач
for job in "${COMPILE_JOBS[@]}"; do
    wait $job || { echo -e "${RED}Ошибка при компиляции${NC}"; exit 1; }
done

# Линковка
echo -e "${BLUE}🔗 Линковка библиотеки...${NC}"
if [ $STATIC -eq 1 ]; then
    ar rcs "$TARGET_LIB" *.o || { echo -e "${RED}Ошибка линковки${NC}"; exit 1; }
    echo -e "${GREEN}✅ Статическая библиотека создана:${NC} $TARGET_LIB"
else
    gcc *.o $LDFLAGS -o "$TARGET_LIB" $LIBS || { echo -e "${RED}Ошибка линковки${NC}"; exit 1; }
    echo -e "${GREEN}✅ Динамическая библиотека создана:${NC} $TARGET_LIB"
fi

# Финальное сообщение
echo ""
echo -e "${GREEN}✨ Сборка успешно завершена!${NC}"
echo ""
echo -e "${BLUE}💡 Советы по использованию:${NC}"
echo -e "   • Запуск тестов:"
echo -e "     ${YELLOW}gcc mendrive_libtest.c -o mendrive_libtest $CFLAGS $LIBS -ldl && ./mendrive_libtest${NC}"
echo -e ""
echo -e "   • Просмотр логов в реальном времени (для уровней debug/trace):"
echo -e "     ${YELLOW}./mendrive_libtest 2>&1 | grep '\\[DEBUG\\]\\|\\[TRACE\\]'${NC}"
echo -e ""
echo -e "   • Проверка символов в библиотеке:"
echo -e "     ${YELLOW}nm -D ${TARGET_LIB} | grep det_${NC}"