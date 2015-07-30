.PHONY: all

# Shared flags etc.
include ../Makefile.conf
ifdef LIBIGL_USE_STATIC_LIBRARY
	LIBIGL_LIB+=-liglembree
endif

all: obj example

.PHONY: example

CFLAGS+=-g -std=c++11

INC=$(LIBIGL_INC) $(ANTTWEAKBAR_INC) $(EIGEN3_INC) $(GLUT_INC) $(EMBREE_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB) $(LIBIGL_LIB) $(EMBREE_LIB)

CPP_FILES=$(wildcard ./*.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o))) 

example: obj $(OBJ_FILES)
	g++ $(OPTFLAGS) $(OPENMP) $(AFLAGS) $(CFLAGS) -o example $(OBJ_FILES) $(LIB)

obj:
	mkdir -p obj

obj/%.o: %.cpp
	g++ $(OPTFLAGS) $(OPENMP) $(AFLAGS) $(CFLAGS) -c $< -o $@ $(INC)

obj/%.o: %.cpp %.h
	g++ $(OPTFLAGS) $(OPENMP) $(AFLAGS) $(CFLAGS) -c $< -o $@ $(INC)

clean:
	rm -f $(OBJ_FILES)
	rm -f example
