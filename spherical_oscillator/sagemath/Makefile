CC=gcc
CFLAGS=-I. -O0 -g -Werror -Wformat # -DUSE_DEBUG
DEPS = hellomake.h
OBJ = tzap.o integrate.o calc.o backtrace.o main.o test_v0.o test_v1.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

testtzap: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -lm

clean:
	rm -f $(OBJ)
	rm -f testtzap
