
.PHONY: all

# Shared flags etc.
include ../../build/Makefile.conf

all: example

.PHONY: example

LIBIGL=../../
LIBIGL_INC=-I$(LIBIGL)/include
LIBIGL_LIB=-L$(LIBIGL)/lib -ligl -liglpng

EIGEN3_INC=-I/opt/local/include/eigen3 -I/opt/local/include/eigen3/unsupported

# YIMAGE Library
YIMG=$(LIBIGL)/external/yimg/
YIMG_LIB=-L$(YIMG) -lyimg -lz -L/opt/local/lib -lpng -bind_at_load
YIMG_INC=-I/usr/X11/include -I$(YIMG) 

ANTTWEAKBAR_INC=-I$(LIBIGL)/external/AntTweakBar/include
ANTTWEAKBAR_LIB=-L$(LIBIGL)/external/AntTweakBar/lib -lAntTweakBar -framework AppKit
INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC) $(YIMG_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB) $(YIMG_LIB)

example: example.o
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) $(LIB) -o example example.o 

example.o: example.cpp
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -c example.cpp -o example.o $(INC)
clean:
	rm -f example.o
	rm -f example
