
.PHONY: all

# Shared flags etc.

all: example 

.PHONY:  example

LIBIGL=/usr/local/igl/libigl/
include $(LIBIGL)/Makefile.conf
LIBIGL_INC=-I$(LIBIGL)/include
LIBIGL_LIB=-L$(LIBIGL)/lib -ligl -liglembree

EIGEN3_INC=-I/opt/local/include/eigen3 -I/opt/local/include/eigen3/unsupported
CFLAGS+=-std=c++11 -g

ANTTWEAKBAR_INC=-I$(LIBIGL)/external/AntTweakBar/include
ANTTWEAKBAR_LIB=-L$(LIBIGL)/external/AntTweakBar/lib -lAntTweakBar -framework AppKit
INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB)

example: example.o
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -o example example.o $(LIB)

example.o: example.cpp
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -c example.cpp -o example.o $(INC)
clean:
	rm -f example.o
	rm -f example
