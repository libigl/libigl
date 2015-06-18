.PHONY: all

# Shared flags etc.
include ../Makefile.conf
LIBIGL_LIB+=-liglbbw -liglcgal

all: obj skeleton-poser

.PHONY: skeleton-poser

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
MOSEK=/usr/local/mosek
MOSEK_INC=-I$(MOSEK)/$(MOSEKVERSION)/tools/platform/$(MOSEKPLATFORM)/h
MOSEK_LIB=-L$(MOSEK)/$(MOSEKVERSION)/tools/platform/$(MOSEKPLATFORM)/bin -lmosek64
LIBIGL_LIB+=-liglmosek 
endif

CGAL=/opt/local/
CGAL_LIB=-L$(CGAL)/lib -lCGAL -lCGAL_Core -lgmp -lmpfr -lboost_thread-mt -lboost_system-mt
CGAL_INC=-I$(CGAL)/include -I/usr/include/
# This is absolutely necessary for Exact Construction
CGAL_FLAGS=-frounding-math -fsignaling-nans 
CFLAGS+=$(CGAL_FLAGS)

INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC) $(GLUT_INC) ${CGAL_INC} ${TETGEN_INC} $(MOSEK_INC) $(EMBREE_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB) ${CGAL_LIB} $(TETGEN_LIB) $(MOSEK_LIB) $(EMBREE_LIB)

CPP_FILES=$(wildcard ./*.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o))) 

CFLAGS+=-std=c++11 -g

skeleton-poser: obj $(OBJ_FILES)
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -o skeleton-poser $(LIB) $(OBJ_FILES) 

obj:
	mkdir -p obj

obj/%.o: %.cpp
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -c $< -o $@ $(INC)

obj/%.o: %.cpp %.h
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -c $< -o $@ $(INC)

clean:
	rm -f $(OBJ_FILES)
	rm -f skeleton-poser
