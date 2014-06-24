.PHONY: all

all: example_static example_header_only

# Shared flags etc.
include ../../build/Makefile.conf

igl_lib=../../

CFLAGS=-g
inc=-I$(igl_lib)/include

STATIC_LIB=-ligl -L../../lib
HEADER_ONLY_DEF=-DIGL_HEADER_ONLY

example_static: example_static.o
	#g++ $(CFLAGS) -o example_static example_static.o $(lib) $(STATIC_LIB) ../../obj/example_fun.o
	g++ $(CFLAGS) -o example_static example_static.o ../../obj/example_fun.o

example_static.o: example.cpp
	g++ $(CFLAGS) -c example.cpp -o example_static.o $(inc)

example_header_only: example_header_only.o
	g++ $(CFLAGS) -o example_header_only example_header_only.o $(lib)

example_header_only.o: example.cpp
	g++ $(CFLAGS) -c example.cpp -o example_header_only.o $(inc) $(HEADER_ONLY_DEF)

clean:
	rm -f example_static.o
	rm -f example_static
	rm -f example_header_only.o
	rm -f example_header_only
