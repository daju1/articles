CC=gcc
CFLAGS=-I. -fPIC -O0 -g -Werror -Wformat # -DUSE_DEBUG
DEPS = hellomake.h
OBJ = lw_rotate.o lw_tensor.o lw_vimanic.o main.o
DROBTASK = lw_drobyshev_task.o lw_drobyshev_task_tensor.o
DROBTEST = lw_drobyshev_task.o lw_drobyshev_task_tensor.o lw_drobyshev_main.o

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

test: $(OBJ)
	$(CC) -o $@ $^ $(CFLAGS) -lm

libtest: $(DROBTEST)
	$(CC) -o $@ $^ $(LDFLAGS) $(LDLIBS)  -lm

mylib: $(DROBTASK)
	$(CC) -shared -o $@.so $^ $(LDFLAGS) $(LDLIBS)

clean:
	rm -f $(OBJ)
	rm -f testtzap
