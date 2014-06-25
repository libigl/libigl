
all: example1

# Shared flags etc.
include ../../build/Makefile.conf

igl_lib=../../
eigen_lib=$(DEFAULT_PREFIX)/include/eigen3/

CFLAGS+=-g
inc=-I$(igl_lib)/include -I$(eigen_lib)
lib=-L$(igl_lib)/lib -ligl

example1: example1.o
	g++ $(CFLAGS) -o example1 example1.o $(lib)

example1.o: example1.cpp
	g++ $(CFLAGS) -c example1.cpp -o example1.o $(inc)


clean:
	rm -f example1.o
	rm -f example1
	rm -f bunny_out.off
