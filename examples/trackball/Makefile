
.PHONY: all

all: example

# Shared flags etc.
include ../../Makefile.conf

.PHONY: example

igl_lib=../../

CFLAGS=-g -Wall 
#deps=-MMD -MF depends.txt
EIGEN3_INC=-I$(DEFAULT_PREFIX)/include/eigen3 -I$(DEFAULT_PREFIX)/include/eigen3/unsupported
inc=-I$(igl_lib)/include $(EIGEN3_INC)
lib=$(OPENGL_LIB) $(GLUT_LIB) -L$(igl_lib)/lib -ligl

example: example.o
	g++ $(CFLAGS) -o example example.o $(lib)

example.o: example.cpp
	g++ $(CFLAGS) $(deps) -c example.cpp -o example.o $(inc)
clean:
	rm -f example.o
	rm -f example
	rm -f depends.txt
