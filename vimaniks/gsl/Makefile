# CC=gcc
CC=c99
CFLAGS=-I. -I./local/include  -fPIC -O0 -g -Werror -Wformat # -DUSE_DEBUG
LDLIBS= -Wl,-Bstatic -lgsl -lgslcblas -Wl,-Bdynamic -lm
#LDLIBS=-lm
LDFLAGS=#-L./local/lib

VANDERPOL = vanderpol.o
PENDULUM = pendulum.o
PENDULUMLIB = pendulum_lib.o pendulum_tensor.o
PENDULUMLIBTEST = $(PENDULUMLIB) pendulum_lib_test.o

test: vanderpol pendulum pendulumlib pendulumlibtest

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

vanderpol: $(VANDERPOL)
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)

pendulum: $(PENDULUM)
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)

pendulumlib: $(PENDULUMLIB)
	$(CC)  -shared  -o $@.so $^ $(LDFLAGS) $(LDLIBS)

pendulumlibtest: $(PENDULUMLIBTEST)
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(OBJ)
	rm -f vanderpol
	rm -f pendulum
	rm -f libtest
	rm -f libtest
	rm -f pendulumlib.so
	rm *.o