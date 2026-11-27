# Makefile -- автономная версия МенДрайв (MenDriveCpp/standalone)
#
#   make            -> ./mendrive          (std::thread, без внешних зависимостей)
#   make omp        -> ./mendrive_omp      (OpenMP; на macOS нужен brew install libomp)
#   make lib        -> libmendrive_new.*   (новое ядро для старой ctypes-обёртки)
#   make reflib     -> libmendrive_ref.*   (ИСХОДНОЕ ../mendrive_core.cpp -- для сравнения)
#   make test       -> сравнение с исходным ctypes-путём + проверка новых мод и экспорта полей (python3 + numpy)
#   make test-fields-> только новые моды (MODE) и экспорт полей (FIELDS): tests/test_fields_linear.py
#   make smoke      -> smoke-тест Preisach2D, замер памяти, проверка --resume
#   make plot       -> прогон с экспортом полей + графики (CFG=..., OUT=..., ARGS="key=val ...")
#   make plot-linear-> линейный эталон TE + сравнение с точным решением
#   make clean

CXX      ?= c++
UNAME_S  := $(shell uname -s)

# -march=native (x86, g++/clang) или -mcpu=native (Apple Silicon); если компилятор
# не понимает ни один из флагов -- собираем без них. Можно задать вручную: make ARCH=
ARCH ?= $(shell echo 'int main(){}' | $(CXX) -march=native -x c++ - -o /dev/null >/dev/null 2>&1 && echo -march=native || \
          (echo 'int main(){}' | $(CXX) -mcpu=native -x c++ - -o /dev/null >/dev/null 2>&1 && echo -mcpu=native))

BASEFLAGS = -O3 $(ARCH) -std=c++17 -Wall -Wextra
FLAGS     = $(BASEFLAGS) $(CXXFLAGS)
LDLIBS   += -pthread -lstdc++fs

ifeq ($(UNAME_S),Darwin)
  SHLIB_EXT = dylib
  SHARED    = -dynamiclib -fPIC
  LIBOMP   ?= $(shell brew --prefix libomp 2>/dev/null || echo /opt/homebrew/opt/libomp)
  OMPFLAGS  = -Xpreprocessor -fopenmp -I$(LIBOMP)/include
  OMPLIBS   = -L$(LIBOMP)/lib -lomp
else
  SHLIB_EXT = so
  SHARED    = -shared -fPIC
  OMPFLAGS  = -fopenmp
  OMPLIBS   = -fopenmp
endif

.PHONY: all omp lib reflib test test-fields smoke plot plot-linear clean

all: mendrive

mendrive: main.cpp mendrive_core.cpp
	$(CXX) $(FLAGS) main.cpp -o $@ $(LDLIBS)

omp: main.cpp mendrive_core.cpp
	$(CXX) $(FLAGS) -DUSE_OPENMP $(OMPFLAGS) main.cpp -o mendrive_omp $(OMPLIBS) $(LDLIBS)

lib: libmendrive_new.$(SHLIB_EXT)
libmendrive_new.$(SHLIB_EXT): mendrive_core.cpp
	$(CXX) $(FLAGS) $(SHARED) mendrive_core.cpp -o $@ $(LDLIBS)

reflib: libmendrive_ref.$(SHLIB_EXT)
libmendrive_ref.$(SHLIB_EXT): ../mendrive_core.cpp
	$(CXX) $(FLAGS) $(SHARED) ../mendrive_core.cpp -o $@

test: mendrive lib reflib
	python3 tests/compare_ctypes.py
	python3 tests/test_fields_linear.py

test-fields: mendrive
	python3 tests/test_fields_linear.py

smoke: mendrive
	sh tests/smoke.sh

CFG ?=
OUT ?= out/fields
ARGS ?=
plot: mendrive
	./mendrive $(if $(CFG),-c $(CFG),) out_dir=$(OUT) $(ARGS)
	python3 plot.py $(OUT)

plot-linear: mendrive
	./mendrive -c examples/linear_check_TE.cfg out_dir=out/linear_check_TE
	python3 linear1d.py out/linear_check_TE
	python3 plot.py out/linear_check_TE --compare out/linear_check_TE/linear_analytic_n1.csv --fit none --norm none

clean:
	rm -f mendrive mendrive_omp libmendrive_new.* libmendrive_ref.*
	rm -rf tests/out_* tests/time_*.txt
