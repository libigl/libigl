.PHONY: all

# Shared flags etc.

all: upright

.PHONY: upright

include ../Makefile.conf

INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB)

CFLAGS+=-g $(OPTFLAGS)

upright: upright.o
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -o upright upright.o $(LIB)

upright.o: example.cpp
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -c example.cpp -o upright.o $(INC)
clean:
	rm -f upright.o
	rm -f upright
