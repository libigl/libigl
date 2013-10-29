
.PHONY: all

# Shared flags etc.
include ../../Makefile.conf

all: example

.PHONY: example

LIBIGL=../../
LIBIGL_INC=-I$(LIBIGL)/include
LIBIGL_LIB=-L$(LIBIGL)/lib -ligl -liglembree

EIGEN3_INC=-I/opt/local/include/eigen3 -I/opt/local/include/eigen3/unsupported

#EMBREE=$(LIBIGL)/external/embree
#EMBREE_INC=-I$(EMBREE)/rtcore -I$(EMBREE)/common
#EMBREE_LIB=-L$(EMBREE)/bin -lrtcore -lsys
EMBREE=$(LIBIGL)/external/embree
EMBREE_INC=-I$(EMBREE)/ -I$(EMBREE)/embree
EMBREE_LIB=-L$(EMBREE)/bin -lembree -lsys

ANTTWEAKBAR_INC=-I$(LIBIGL)/external/AntTweakBar/include
ANTTWEAKBAR_LIB=-L$(LIBIGL)/external/AntTweakBar/lib -lAntTweakBar -framework AppKit
INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC) $(EMBREE_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB) $(EMBREE_LIB)

example: example.o
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -o example example.o $(LIB)

example.o: example.cpp
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -c example.cpp -o example.o $(INC)
clean:
	rm -f example.o
	rm -f example
