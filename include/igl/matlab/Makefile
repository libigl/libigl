include ../../../Makefile.conf
all: CFLAGS += -O3 -DNDEBUG
debug: CFLAGS += -g -Wall -Werror

.PHONY: all
all: libiglmatlab
debug: libiglmatlab

.PHONY: libmatlab
libiglmatlab: obj ../../../lib/libiglmatlab.a

CPP_FILES=$(wildcard *.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

# include igl headers
INC+=-I../../../include/

# EXPECTS THAT CFLAGS IS ALREADY SET APPROPRIATELY 

# Eigen dependency
EIGEN3_INC=-I$(DEFAULT_PREFIX)/include/eigen3 -I$(DEFAULT_PREFIX)/include/eigen3/unsupported
INC+=$(EIGEN3_INC)

# Matlab dependency
ifndef MATLAB
  MATLAB=/Applications/MATLAB_R2013b.app/
  $(warning MATLAB undefined. Setting MATLAB=${MATLAB})
endif
MATLAB_INC=-I$(MATLAB)/extern/include/
MATLAB_LIB=-L$(MATLAB)/bin/maci64 -lmx -leng
INC+=$(MATLAB_INC)

obj: 
	mkdir -p obj

../../../lib/libiglmatlab.a: $(OBJ_FILES)
	rm -f $@
	ar cqs $@ $(OBJ_FILES)

obj/%.o: %.cpp %.h
	g++ $(AFLAGS) $(CFLAGS) -c -o $@ $< $(INC)

clean:
	rm -f obj/*.o
	rm -f ../../../lib/libiglmatlab.a
