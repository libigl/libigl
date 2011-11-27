
.PHONY: all

all: example

.PHONY: example

igl_lib=-I../../
eigen=-I/usr/local/include/eigen3 -I/usr/local/include/eigen3/unsupported

#CFLAGS=-g
CFLAGS=-Os -DNDEBUG
inc=$(igl_lib) $(eigen)
lib=

example: example.o
	g++ $(CFLAGS) -o example example.o $(lib)
	rm example.o

example.o: example.cpp
	g++ $(CFLAGS) -c example.cpp -o example.o $(inc)
clean:
	rm -f example.o
	rm -f example
