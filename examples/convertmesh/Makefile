.PHONY: all

# Shared flags etc.
include ../../build/Makefile.conf

all: convertmesh

.PHONY: convertmesh

INC=$(LIBIGL_INC) $(EIGEN3_INC) 
LIB=$(LIBIGL_LIB) 

convertmesh: convertmesh.o
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) $(LIB) -o convertmesh convertmesh.o 

convertmesh.o: convertmesh.cpp
	g++ $(OPENMP) $(AFLAGS) $(CFLAGS) -c convertmesh.cpp -o convertmesh.o $(INC)
clean:
	rm -f convertmesh.o
	rm -f convertmesh
