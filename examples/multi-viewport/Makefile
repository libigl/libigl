.PHONY: all

# Shared flags etc.

all: example 

.PHONY:  example

include ../Makefile.conf

CFLAGS+=-std=c++11 -g -Wno-deprecated-declarations

INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB)

example: example.o
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -o example example.o $(LIB)

example.o: example.cpp
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -c example.cpp -o example.o $(INC)
clean:
	rm -f example.o
	rm -f example
