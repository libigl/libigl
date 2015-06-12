.PHONY: all

# Shared flags etc.
include ../Makefile.conf
LIBIGL_LIB+=-liglpng

# YIMG dependency

all: example

.PHONY: example

IGL=../../

inc=$(LIBIGL_INC) $(YIMG_INC) $(EIGEN3_INC)
lib=$(LIBIGL_LIB) $(OPENGL_LIB) $(GLUT_LIB) $(YIMG_LIB)

example: example.o
	g++ $(CFLAGS) -o example example.o $(lib)

example.o: example.cpp
	g++ $(CFLAGS) -c example.cpp -o example.o $(inc)
clean:
	rm -f example.o
	rm -f example
