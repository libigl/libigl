.PHONY: all

# Shared flags etc.
include ../Makefile.conf

all: example

.PHONY: example

inc=$(LIBIGL_INC) $(ANTTWEAKBAR_INC)
lib=$(LIBIGL_LIB) $(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB)
CFLAGS+=-g

example: example.o
	g++ $(CFLAGS) -o example example.o $(lib)

example.o: example.cpp
	g++ $(CFLAGS) -c example.cpp -o example.o $(inc)
clean:
	rm -f example.o
	rm -f example
