
all: example1

# Shared flags etc.
include ../../build/Makefile.conf

igl_lib=../../
ifndef EIGEN3_INC
	EIGEN3_INC=-I$(DEFAULT_PREFIX)/include/eigen3/
endif

CFLAGS+=-g -fopenmp
inc=-I$(igl_lib)/include $(EIGEN3_INC)
lib=-L$(igl_lib)/lib -ligl

example1: example1.o
	g++ $(CFLAGS) -o example1 example1.o $(lib)

example1.o: example1.cpp
	g++ $(CFLAGS) -c example1.cpp -o example1.o $(inc)


clean:
	rm -f example1.o
	rm -f example1
	rm -f bunny_out.off
