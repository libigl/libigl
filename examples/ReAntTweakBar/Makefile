
.PHONY: all

all: example

.PHONY: example

igl_lib=../../

lib=-lAntTweakBar -framework AppKit -L$(igl_lib)/lib -ligl -L/opt/local/lib
inc=-I$(igl_lib)/include -I/opt/local/include

example: example.o
	g++ $(CFLAGS) -o example example.o -framework OpenGL -framework GLUT $(lib)

example.o: example.cpp
	g++ $(CFLAGS) -c example.cpp -o example.o $(inc)
clean:
	rm -f example.o
	rm -f example
