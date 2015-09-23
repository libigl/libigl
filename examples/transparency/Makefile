.PHONY: all

# Shared flags etc.
include ../Makefile.conf

all: example

.PHONY: example

CFLAGS+=-O3 

INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB)

example: example.o
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -o example example.o $(LIB)

example.o: example.cpp
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -c example.cpp -o example.o $(INC)
clean:
	rm -f example.o
	rm -f example
