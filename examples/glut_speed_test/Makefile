.PHONY: all

all: example

# Shared flags etc.
include ../Makefile.conf

.PHONY: example

igl_lib=../../

CFLAGS+=-g -Wall 
INC=$(LIBIGL_INC)
LIB=$(OPENGL_LIB) $(GLUT_LIB) $(LIBIGL_LIB)

example: example.o
	g++ $(CFLAGS) -o example example.o $(LIB)

example.o: example.cpp
	g++ $(CFLAGS) $(deps) -c example.cpp -o example.o $(INC)

clean:
	rm -f example.o
	rm -f example
