CC=gcc
CFLAGS=-I. -I./local/include -O0 -g -Werror -Wformat # -DUSE_DEBUG
LDLIBS=-lm -lgsl -lgslcblas
LDFLAGS=-L./local/lib

# OBJ = vanderpol.o
OBJ = pendulum.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)  $(LDFLAGS) $(LDLIBS)

test: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS)  $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(OBJ)
	rm -f test
