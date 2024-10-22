# CC=gcc
CC=c99
CFLAGS=-I. -I./local/include  -fPIC -O0 -g -Werror -Wformat # -DUSE_DEBUG
LDLIBS= -Wl,-Bstatic -lgsl -lgslcblas -Wl,-Bdynamic -lm
#LDLIBS=-lm
LDFLAGS=#-L./local/lib

# OBJ = vanderpol.o
OBJ = pendulum.o
LIB = pendulum_lib.o pendulum_tensor.o
LIBTEST = $(LIB) pendulum_lib_test.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

test: $(OBJ)
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)

libtest: $(LIBTEST)
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)

mylib: $(LIB)
	$(CC)  -shared  -o $@.so $^ $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(OBJ)
	rm -f test
	rm -f libtest
	rm -f mylib.so
	rm *.o