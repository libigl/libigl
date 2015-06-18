.PHONY: all

# Shared flags etc.
include ../Makefile.conf
LIBIGL_LIB+=-liglpng

all: example

.PHONY: example

INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC) $(YIMG_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB) $(YIMG_LIB)

example: example.o
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) $(LIB) -o example example.o 

example.o: example.cpp
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -c example.cpp -o example.o $(INC)
clean:
	rm -f example.o
	rm -f example
