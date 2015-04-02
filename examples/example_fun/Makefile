.PHONY: all

all: example_static example_header_only

# Shared flags etc.
include ../../build/Makefile.conf

CFLAGS+=-g
INC=$(LIBIGL_INC)

STATIC_LIB=-ligl -L../../lib -DIGL_STATIC_LIBRARY

example_static: example_static.o
	#g++ $(CFLAGS) -o example_static example_static.o $(lib) $(STATIC_LIB) ../../obj/example_fun.o
	g++ $(CFLAGS) -o example_static example_static.o ../../obj/example_fun.o

example_static.o: example.cpp
	g++ $(CFLAGS) -c example.cpp -o example_static.o $(INC)

example_header_only: example_header_only.o
	g++ $(CFLAGS) -o example_header_only example_header_only.o $(lib)

example_header_only.o: example.cpp
	g++ $(CFLAGS) -c example.cpp -o example_header_only.o $(INC)

clean:
	rm -f example_static.o
	rm -f example_static
	rm -f example_header_only.o
	rm -f example_header_only
