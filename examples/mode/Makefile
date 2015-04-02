
.PHONY: all

# Shared flags etc.
include ../../build/Makefile.conf

all: example

.PHONY: example

igl_lib=../../
ifndef EIGEN3_INC
	EIGEN3_INC=$(DEFAULT_PREFIX)/include/eigen3/
endif

CFLAGS+=-g
inc=-I$(igl_lib)/include -I$(EIGEN3_INC)
lib=-L$(igl_lib)/lib -ligl

example: example.o
	g++ $(CFLAGS) -o example example.o $(lib)

example.o: example.cpp
	g++ $(CFLAGS) -c example.cpp -o example.o $(inc)
clean:
	rm -f example.o
	rm -f example
