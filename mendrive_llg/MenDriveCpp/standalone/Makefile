# Makefile -- автономная версия МенДрайв (MenDriveCpp/standalone)
#
#   make            -> ./mendrive          (std::thread, без внешних зависимостей)
#   make omp        -> ./mendrive_omp      (OpenMP; на macOS нужен brew install libomp)
#   make lib        -> libmendrive_new.*   (новое ядро для старой ctypes-обёртки)
#   make reflib     -> libmendrive_ref.*   (ИСХОДНОЕ ../mendrive_core.cpp -- для сравнения)
#   make test       -> сравнение с исходным ctypes-путём (нужны python3 + numpy)
#   make smoke      -> smoke-тест Preisach2D, замер памяти, проверка --resume
#   make clean

CXX      ?= c++
UNAME_S  := $(shell uname -s)

# -march=native (x86, g++/clang) или -mcpu=native (Apple Silicon); если компилятор
# не понимает ни один из флагов -- собираем без них. Можно задать вручную: make ARCH=
ARCH ?= $(shell echo 'int main(){}' | $(CXX) -march=native -x c++ - -o /dev/null >/dev/null 2>&1 && echo -march=native || \
          (echo 'int main(){}' | $(CXX) -mcpu=native -x c++ - -o /dev/null >/dev/null 2>&1 && echo -mcpu=native))

BASEFLAGS = -O3 $(ARCH) -std=c++17 -Wall -Wextra
FLAGS     = $(BASEFLAGS) $(CXXFLAGS)
LDLIBS   += -pthread

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

.PHONY: all omp lib reflib test smoke clean

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

smoke: mendrive
	sh tests/smoke.sh

clean:
	rm -f mendrive mendrive_omp libmendrive_new.* libmendrive_ref.*
	rm -rf tests/out_* tests/time_*.txt
