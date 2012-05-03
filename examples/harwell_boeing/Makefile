
.PHONY: all

# Shared flags etc.
include ../../Makefile.conf

all: example

.PHONY: example

igl_lib=../../
eigen=-I$(OPT)/local/include/eigen3 -I$(OPT)/local/include/eigen3/unsupported

CFLAGS=-g
inc=-I$(igl_lib)/include $(eigen)
lib=-L$(igl_lib)/lib -ligl

example: example.o
	g++ $(CFLAGS) -o example example.o $(lib)

example.o: example.cpp
	g++ $(CFLAGS) -c example.cpp -o example.o $(inc)
clean:
	rm -f example.o
	rm -f example
