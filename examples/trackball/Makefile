.PHONY: all

all: example

# Shared flags etc.
include ../Makefile.conf

.PHONY: example

igl_lib=../../

CFLAGS+=-g -Wall 
#deps=-MMD -MF depends.txt
inc=$(LIBIGL_INC) $(EIGEN3_INC)
lib=$(LIBIGL_LIB) $(OPENGL_LIB) $(GLUT_LIB)

example: example.o
	g++ $(CFLAGS) -o example example.o $(lib)

example.o: example.cpp
	g++ $(CFLAGS) $(deps) -c example.cpp -o example.o $(inc)
clean:
	rm -f example.o
	rm -f example
	rm -f depends.txt
