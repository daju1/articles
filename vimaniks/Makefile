CC=gcc
CFLAGS=-I. -O0 -g -Werror -Wformat # -DUSE_DEBUG
DEPS = hellomake.h
OBJ = lw_rotate.o lw_vimanic.o main.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

test: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -lm

clean:
	rm -f $(OBJ)
	rm -f testtzap
