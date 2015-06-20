.PHONY: all

# Shared flags etc.
include ../Makefile.conf
LIBIGL_LIB+=-liglbbw -liglmosek

all: obj example

.PHONY: example

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

INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC) $(GLUT_INC) ${TETGEN_INC} $(MOSEK_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB) $(TETGEN_LIB) $(MOSEK_LIB)

CPP_FILES=$(wildcard ./*.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o))) 

CFLAGS+=-std=c++11 -g

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
