.PHONY: all

# Shared flags etc.
include ../Makefile.conf

all: example

.PHONY: example

CFLAGS+=-g
INC=$(LIBIGL_INC) $(EIGEN3_INC)
LIB=$(LIBIGL_LIB)

example: example.o
	g++ $(CFLAGS) -o example example.o $(LIB)

example.o: example.cpp
	g++ $(CFLAGS) -c example.cpp -o example.o $(INC)
clean:
	rm -f example.o
	rm -f example
