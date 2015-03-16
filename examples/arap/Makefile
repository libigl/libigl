
.PHONY: all

# Shared flags etc.
include ../../build/Makefile.conf

all: example

.PHONY: example

LIBIGL=../../
LIBIGL_INC=-I$(LIBIGL)/include
LIBIGL_LIB=-L$(LIBIGL)/lib -ligl -liglpng -liglsvd3x3

EIGEN3_INC=-I/opt/local/include/eigen3 -I/opt/local/include/eigen3/unsupported
#CFLAGS+=-DIGL_HEADER_ONLY

EMBREE=$(LIBIGL)/external/embree
EMBREE_INC=-I$(EMBREE)/ -I$(EMBREE)/include
EMBREE_LIB=-L$(EMBREE)/build -lembree -lsys

# YIMAGE Library
YIMG=$(LIBIGL)/external/yimg/
YIMG_LIB=-L$(YIMG) -lyimg -lz -L/opt/local/lib -lpng
YIMG_INC=-I/usr/X11/include -I$(YIMG)

ANTTWEAKBAR_INC=-I$(LIBIGL)/external/AntTweakBar/include
ANTTWEAKBAR_LIB=-L$(LIBIGL)/external/AntTweakBar/lib -lAntTweakBar
UNAME := $(shell uname)
# Apple needs to load the AppKit framework for anttweakbar and maybe bind at
# load for png
ifeq ($(UNAME), Darwin)
	YIMG_LIB+=-bind_at_load
	ANTTWEAKBAR_LIB+=-framework AppKit
endif

# SVD 
SINGULAR_VALUE_DECOMPOSITION_INC=\
	-I$(LIBIGL)/external/Singular_Value_Decomposition/

INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC) $(YIMG_INC) $(SINGULAR_VALUE_DECOMPOSITION_INC)
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
