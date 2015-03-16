
.PHONY: all

# Shared flags etc.
include ../../build/Makefile.conf

all: obj example

.PHONY: example

LIBIGL=../../
LIBIGL_INC=-I$(LIBIGL)/include
LIBIGL_LIB=-L$(LIBIGL)/lib -ligl -liglcgal -ligltetgen -liglbbw -liglmosek

EIGEN3_INC=-I/opt/local/include/eigen3 -I/opt/local/include/eigen3/unsupported

ANTTWEAKBAR_INC=-I$(LIBIGL)/external/AntTweakBar/include
ANTTWEAKBAR_LIB=-L$(LIBIGL)/external/AntTweakBar/lib -lAntTweakBar -framework AppKit

TETGEN=$(LIBIGL)/external/tetgen
TETGEN_LIB=-L$(TETGEN) -ltet 
TETGEN_INC=-I$(TETGEN)

EMBREE=$(LIBIGL)/external/embree
EMBREE_INC=-I$(EMBREE)/ -I$(EMBREE)/embree
EMBREE_LIB=-L$(EMBREE)/build -lembree -lsys

CARBON_LIB=-framework Carbon

# Use free glut for mouse scrolling
#FREE_GLUT=/opt/local/
#FREE_GLUT_INC=-I$(FREE_GLUT)/include
#FREE_GLUT_LIB=-L$(FREE_GLUT)/lib -lglut
GLUT_LIB=-framework GLUT
GLUT_INC=

ifdef IGL_NO_MOSEK
CFLAGS+=-DIGL_NO_MOSEK
else
# Adjust your mosek paths etc. accordingly
ifndef MOSEKPLATFORM
  MOSEKPLATFORM=osx64x86
endif
ifndef MOSEKVERSION
  MOSEKVERSION=7
endif
IGLMOSEK=../mosek/
IGLMOSEK_INC=-I$(IGLMOSEK)/
MOSEK=/usr/local/mosek
MOSEK_INC=-I$(MOSEK)/$(MOSEKVERSION)/tools/platform/$(MOSEKPLATFORM)/h
MOSEK_LIB=-L$(MOSEK)/$(MOSEKVERSION)/tools/platform/$(MOSEKPLATFORM)/bin -lmosek64
endif

CGAL=/opt/local/
CGAL_LIB=-L$(CGAL)/lib -lCGAL -lCGAL_Core -lgmp -lmpfr -lboost_thread-mt -lboost_system-mt
CGAL_INC=-I$(CGAL)/include -I/usr/include/
# This is absolutely necessary for Exact Construction
CGAL_FLAGS=-frounding-math -fsignaling-nans 
CFLAGS+=$(CGAL_FLAGS)

INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC) $(MATLAB_INC) $(GLUT_INC) ${CGAL_INC} ${TETGEN_INC} $(MOSEK_INC) $(EMBREE_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB) $(MATLAB_LIB) ${CGAL_LIB} $(CARBON_LIB) $(TETGEN_LIB) $(MOSEK_LIB) $(EMBREE_LIB)

CPP_FILES=$(wildcard ./*.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o))) 

CFLAGS+=-std=c++11 -g

example: obj $(OBJ_FILES)
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -o example $(LIB) $(OBJ_FILES) 

obj:
	mkdir -p obj

obj/%.o: %.cpp
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -c $< -o $@ $(INC)

obj/%.o: %.cpp %.h
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -c $< -o $@ $(INC)

clean:
	rm -f $(OBJ_FILES)
	rm -f example
