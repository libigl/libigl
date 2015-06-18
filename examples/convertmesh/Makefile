.PHONY: all

# Shared flags etc.
include ../Makefile.conf

all: convertmesh

.PHONY: convertmesh

INC=$(LIBIGL_INC) $(EIGEN3_INC) 
LIB=$(LIBIGL_LIB) 
OPTFLAGS+=-O3 -DNDEBUG $(OPENMP)

convertmesh: convertmesh.o
	g++ $(OPTFLAGS) $(AFLAGS) $(CFLAGS) $(LIB) -o convertmesh convertmesh.o 

convertmesh.o: convertmesh.cpp
	g++ $(OPTFLAGS) $(AFLAGS) $(CFLAGS) -c convertmesh.cpp -o convertmesh.o $(INC)
clean:
	rm -f convertmesh.o
	rm -f convertmesh
