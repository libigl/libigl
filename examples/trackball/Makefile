
.PHONY: all

all: example

.PHONY: example

igl_lib=../../

CFLAGS=-g -Wall 
#deps=-MMD -MF depends.txt
inc=-I$(igl_lib)/include
lib=-L$(igl_lib)/lib -ligl

example: example.o
	g++ $(CFLAGS) -o example example.o -framework OpenGL -framework GLUT $(lib)
	rm example.o

example.o: example.cpp
	g++ $(CFLAGS) $(deps) -c example.cpp -o example.o $(inc)
clean:
	rm -f example.o
	rm -f example
	rm -f depends.txt
