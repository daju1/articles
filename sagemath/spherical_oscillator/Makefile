CC=gcc
CFLAGS=-I. -DUSE_DEBUG
DEPS = hellomake.h
OBJ = testtzap.o tzap.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

testtzap: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -lm

clean:
	rm -f $(OBJ)
	rm -f testtzap
