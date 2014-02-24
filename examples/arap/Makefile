
.PHONY: all

# Shared flags etc.
include ../../Makefile.conf

all: example

.PHONY: example

LIBIGL=../../
LIBIGL_INC=-I$(LIBIGL)/include
LIBIGL_LIB=-L$(LIBIGL)/lib -ligl -liglpng -liglsvd3x3

EIGEN3_INC=-I/opt/local/include/eigen3 -I/opt/local/include/eigen3/unsupported

#EMBREE=$(LIBIGL)/external/embree
#EMBREE_INC=-I$(EMBREE)/rtcore -I$(EMBREE)/common
#EMBREE_LIB=-L$(EMBREE)/bin -lrtcore -lsys
EMBREE=$(LIBIGL)/external/embree
EMBREE_INC=-I$(EMBREE)/ -I$(EMBREE)/embree
EMBREE_LIB=-L$(EMBREE)/bin -lembree -lsys

# YIMAGE Library
YIMG=$(LIBIGL)/external/yimg/
YIMG_LIB=-L$(YIMG) -lyimg -lz -L/usr/X11/lib -lpng -bind_at_load
YIMG_INC=-I/usr/X11/include -I$(YIMG)

ANTTWEAKBAR_INC=-I$(LIBIGL)/external/AntTweakBar/include
ANTTWEAKBAR_LIB=-L$(LIBIGL)/external/AntTweakBar/lib -lAntTweakBar -framework AppKit

INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC) $(YIMG_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB) $(YIMG_LIB)

example: example.o
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) $(LIB) -o example $^

%.o: %.cpp
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -c $< -o $@ $(INC)

%.o: %.cpp %.h
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -c $< -o $@ $(INC)

clean:
	rm -f example.o
	rm -f example
