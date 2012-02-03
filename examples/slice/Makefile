
.PHONY: all

all: example

.PHONY: example

igl_lib=../../
eigen=-I/usr/local/include/eigen3 -I/usr/local/include/eigen3/unsupported

#CFLAGS=-g
CFLAGS=-Os -DNDEBUG
inc=-I$(igl_lib) $(eigen)
lib=-L$(igl_lib) -ligl

example: example.o
	g++ $(CFLAGS) -o example example.o $(lib)
	rm example.o

example.o: example.cpp
	g++ $(CFLAGS) -c example.cpp -o example.o $(inc)
clean:
	rm -f example.o
	rm -f example
