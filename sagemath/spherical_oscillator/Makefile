CC=gcc
CFLAGS=-I. -O0 -g # -DUSE_DEBUG
DEPS = hellomake.h
OBJ = testtzap.o tzap.o backtrace.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

testtzap: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -lm

clean:
	rm -f $(OBJ)
	rm -f testtzap
