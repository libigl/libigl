
.PHONY: all

all: example

# Shared flags etc.
include ../../Makefile.conf

.PHONY: example

igl_lib=../../

CFLAGS=-g -Wall 
inc=-I$(igl_lib)/include
lib=$(OPENGL_LIB) $(GLUT_LIB) -L$(igl_lib)/lib -ligl

example: example.o
	g++ $(CFLAGS) -o example example.o $(lib)

example.o: example.cpp
	g++ $(CFLAGS) $(deps) -c example.cpp -o example.o $(inc)
clean:
	rm -f example.o
	rm -f example
