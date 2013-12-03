
.PHONY: all

# Shared flags etc.

all: upright

.PHONY: upright

LIBIGL=/usr/local/igl/libigl/
include $(LIBIGL)/Makefile.conf
LIBIGL_INC=-I$(LIBIGL)/include
LIBIGL_LIB=-L$(LIBIGL)/lib -ligl

EIGEN3_INC=-I/opt/local/include/eigen3 -I/opt/local/include/eigen3/unsupported

CARBON_LIB=-framework Carbon

ANTTWEAKBAR_INC=-I$(LIBIGL)/external/AntTweakBar/include
ANTTWEAKBAR_LIB=-L$(LIBIGL)/external/AntTweakBar/lib -lAntTweakBar -framework AppKit
INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB) $(CARBON_LIB)

upright: upright.o
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -o upright upright.o $(LIB)

upright.o: example.cpp
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -c example.cpp -o upright.o $(INC)
clean:
	rm -f upright.o
	rm -f upright
