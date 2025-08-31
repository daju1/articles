IDIR =$(PWD)/local/include
$(info IDIR=$(IDIR))

CC=gcc
CFLAGS=-I$(IDIR) -g -O0 -DINCLUDE_CUBA_H
CXXFLAGS=-I$(IDIR) -g -O0 -DINCLUDE_CUBA_H


LDIR =../lib

LIBS=-lm -lcuba -L/$(PWD)/local/lib

_DEPS = cuba.h
DEPS = $(patsubst %,$(IDIR)/%,$(_DEPS))
$(info DEPS=$(DEPS))

%.o: %.c $(DEPS)
	$(CC) -c -o $@ $< $(CFLAGS)

.PHONY: all clean

clean:
	rm -f *.o calc_sphere_mass calc_proton_mass calc_longitudinal_sphere_mass calc_transverse_sphere_mass calc_gradient_sphere_mass
