
.PHONY: all

all: example

.PHONY: example

igl_lib=../../
eigen=/opt/local/include/eigen3/

CFLAGS=-g
inc=-I$(igl_lib)/include -I$(eigen)
lib=-L$(igl_lib)/lib -ligl

example: example.o
	g++ $(CFLAGS) -o example example.o $(lib)
	rm example.o

example.o: example.cpp
	g++ $(CFLAGS) -c example.cpp -o example.o $(inc)
clean:
	rm -f example.o
	rm -f example
