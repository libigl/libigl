.PHONY: all

# Shared flags etc.
include ../Makefile.conf
LIBIGL_LIB+=-liglmatlab

MATLAB_INC=-I$(MATLAB)/extern/include/
MATLAB_LIB=-L$(MATLAB)/bin/maci64 -lmx -lmat

LIB+=$(LIBIGL_LIB) $(MATLAB_LIB)
INC+=$(LIBIGL_INC) $(EIGEN3_INC) $(MATLAB_INC)

CFLAGS += -g

ifeq ($(IGL_WITH_MATLAB),1)
	EXAMPLE=example
endif

all: $(EXAMPLE)

example: example.o
	g++ -o example example.o $(LIB)

example.o: example.cpp
	g++ $(CFLAGS) -o $@ -c $< $(INC)

clean:
	rm -f example.o
	rm -f example
