.PHONY: all

# Shared flags etc.
include ../Makefile.conf
LIBIGL_LIB+=-liglcgal

all: example

.PHONY: example

CGAL=/opt/local/
CGAL_LIB=-L$(CGAL)/lib -lCGAL -lCGAL_Core -lgmp -lmpfr -lboost_thread-mt -lboost_system-mt
CGAL_INC=-I$(CGAL)/include -I/usr/include/
# This is absolutely necessary for Exact Construction
CGAL_FLAGS=-frounding-math -fsignaling-nans 

INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC) $(EMBREE_INC) $(CGAL_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB) $(EMBREE_LIB) $(CGAL_LIB)
CFLAGS+=$(CGAL_FLAGS)

example: example.o
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) $(LIB) -o example example.o 

example.o: example.cpp
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -c example.cpp -o example.o $(INC)
clean:
	rm -f example.o
	rm -f example
