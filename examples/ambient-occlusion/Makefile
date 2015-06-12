.PHONY: all

# Shared flags etc.
include ../Makefile.conf
LIBIGL_LIB+=-liglembree

all: example

.PHONY: example

EMBREE=$(LIBIGL)/external/embree
EMBREE_INC=-I$(EMBREE)/ -I$(EMBREE)/include
EMBREE_LIB=-L$(EMBREE)/build -lembree -lsys

INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC) $(EMBREE_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB) $(EMBREE_LIB)

example: example.o
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) $(LIB) -o example example.o 

example.o: example.cpp
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -c example.cpp -o example.o $(INC)
clean:
	rm -f example.o
	rm -f example
