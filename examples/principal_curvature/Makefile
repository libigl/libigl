
.PHONY: all

# Shared flags etc.
include ../../Makefile.conf

all: curvature

.PHONY: curvature

igl_lib=../../
eigen=-I$(DEFAULT_PREFIX)/include/eigen3

CFLAGS=-g
inc=-I$(igl_lib)/include $(eigen)

curvature: curvature.o
	g++ $(CFLAGS) -o curvature curvature.o $(lib)

curvature.o: curvature.cpp
	g++ $(CFLAGS) -c curvature.cpp -o curvature.o $(inc)
clean:
	rm -f curvature.o
	rm -f curvature
