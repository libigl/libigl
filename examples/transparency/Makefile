
.PHONY: all

# Shared flags etc.
include ../../build/Makefile.conf

all: example

.PHONY: example

LIBIGL=../../
LIBIGL_INC=-I$(LIBIGL)/include
LIBIGL_LIB=-L$(LIBIGL)/lib -ligl -liglmatlab

EIGEN3_INC=-I/opt/local/include/eigen3 -I/opt/local/include/eigen3/unsupported
CFLAGS+=-O3 -fopenmp

ANTTWEAKBAR_INC=-I$(LIBIGL)/external/AntTweakBar/include
ANTTWEAKBAR_LIB=-L$(LIBIGL)/external/AntTweakBar/lib -lAntTweakBar -framework AppKit

MATLAB_INC=-I$(MATLAB)/extern/include/
MATLAB_LIB=-L$(MATLAB)/bin/maci64 -lmx -lmat -lmex -lstdc++

INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC) $(MATLAB_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB) $(MATLAB_LIB)

example: example.o
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -o example example.o $(LIB)

example.o: example.cpp
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -c example.cpp -o example.o $(INC)
clean:
	rm -f example.o
	rm -f example
