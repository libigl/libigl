.PHONY: all

# Shared flags etc.
include ../Makefile.conf

all: example

.PHONY: example

#CFLAGS+=-g
CFLAGS+=-Os -DNDEBUG
inc=$(LIBIGL_INC) $(EIGEN3_INC)
lib=$(LIBIGL_LIB)

example: example.o
	g++ $(CFLAGS) -o example example.o $(lib)

example.o: example.cpp
	g++ $(CFLAGS) -c example.cpp -o example.o $(inc)
clean:
	rm -f example.o
	rm -f example
