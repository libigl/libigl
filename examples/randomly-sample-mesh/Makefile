
.PHONY: all

# Shared flags etc.
include ../../build/Makefile.conf

all: obj example

.PHONY: example

LIBIGL=../../
LIBIGL_INC=-I$(LIBIGL)/include
LIBIGL_LIB=-L$(LIBIGL)/lib -ligl -liglmatlab

EIGEN3_INC=-I/opt/local/include/eigen3 -I/opt/local/include/eigen3/unsupported

ANTTWEAKBAR_INC=-I$(LIBIGL)/external/AntTweakBar/include
ANTTWEAKBAR_LIB=-L$(LIBIGL)/external/AntTweakBar/lib -lAntTweakBar -framework AppKit

MATLAB_INC=-I$(MATLAB)/extern/include/
MATLAB_LIB=-L$(MATLAB)/bin/maci64 -lmx -lmat -lmex -lstdc++

CARBON_LIB=-framework Carbon

# Use free glut for mouse scrolling
#FREE_GLUT=/opt/local/
#FREE_GLUT_INC=-I$(FREE_GLUT)/include
#FREE_GLUT_LIB=-L$(FREE_GLUT)/lib -lglut
GLUT_LIB=-framework GLUT
GLUT_INC=

INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC) $(MATLAB_INC) $(GLUT_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB) $(MATLAB_LIB) $(CARBON_LIB)

CPP_FILES=$(wildcard ./*.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o))) 

CFLAGS+=-std=c++11

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
