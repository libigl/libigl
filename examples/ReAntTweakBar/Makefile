
.PHONY: all

all: example

.PHONY: example

STATIC_ANTTWEAKBAR=""

igl_lib=../../

ifdef STATIC_ANTTWEAKBAR
CFLAGS=-g -DSTATIC_ANTTWEAKBAR
inc=-I$(igl_lib) -I./static
lib=-Lstatic -lAntTweakBar -framework AppKit
else
CFLAGS=-g
inc=-I$(igl_lib)
lib=-lAntTweakBar
endif

example: example.o
	g++ $(CFLAGS) -o example example.o -framework OpenGL -framework GLUT $(lib)
	rm example.o

example.o: example.cpp
	g++ $(CFLAGS) -c example.cpp -o example.o $(inc)
clean:
	rm -f example.o
	rm -f example
