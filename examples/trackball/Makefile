
.PHONY: all

all: example

.PHONY: example

igl_lib=../../

CFLAGS=-g
inc=-I$(igl_lib)

example: example.o
	g++ $(CFLAGS) -o example example.o -framework OpenGL -framework GLUT $(lib)
	rm example.o

example.o: example.cpp
	g++ $(CFLAGS) -c example.cpp -o example.o $(inc)
clean:
	rm -f example.o
	rm -f example
