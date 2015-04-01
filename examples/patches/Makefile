
.PHONY: all

# Shared flags etc.
include ../../build/Makefile.conf

all: obj example

.PHONY: example

CFLAGS+=-g -std=c++11

LIBIGL=../../
LIBIGL_INC=-I$(LIBIGL)/include
LIBIGL_LIB=-L$(LIBIGL)/lib -ligl -liglmatlab -liglembree -liglboost

EIGEN3_INC=-I/opt/local/include/eigen3 -I/opt/local/include/eigen3/unsupported

ANTTWEAKBAR_INC=-I$(LIBIGL)/external/AntTweakBar/include
ANTTWEAKBAR_LIB=-L$(LIBIGL)/external/AntTweakBar/lib -lAntTweakBar -framework AppKit

EMBREE=$(LIBIGL)/external/embree
EMBREE_INC=-I$(EMBREE)/ -I$(EMBREE)/include
EMBREE_LIB=-L$(EMBREE)/build -lembree -lsys

MATLAB_INC=-I$(MATLAB)/extern/include/
MATLAB_LIB=-L$(MATLAB)/bin/maci64 -lmx -lmat -lmex -lstdc++

# Use free glut for mouse scrolling
GLUT_LIB=-framework GLUT
GLUT_INC=-framework GLUT

INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC) $(MATLAB_INC) $(GLUT_INC) $(EMBREE_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB) $(MATLAB_LIB) $(EMBREE_LIB)

CPP_FILES=$(wildcard ./*.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o))) 

example: obj $(OBJ_FILES)
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -o example $(OBJ_FILES) $(LIB)

obj:
	mkdir -p obj

obj/%.o: %.cpp
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -c $< -o $@ $(INC)

obj/%.o: %.cpp %.h
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -c $< -o $@ $(INC)

clean:
	rm -f $(OBJ_FILES)
	rm -f example
