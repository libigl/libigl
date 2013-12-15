include ../../../Makefile.conf

.PHONY: all
all: libiglsvd3x3
debug: libiglsvd3x3

include ../../../Makefile.conf
all: OPTFLAGS += -O3 -DNDEBUG $(OPENMP)
debug: OPTFLAGS += -g -Wall
CFLAGS += $(OPTFLAGS) -std=c++11

.PHONY: libiglsvd3x3
libiglsvd3x3: obj ../../../lib/libiglsvd3x3.a

CPP_FILES=$(wildcard *.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

# include igl headers
INC+=-I../../../include/

# EXPECTS THAT CFLAGS IS ALREADY SET APPROPRIATELY 


# SVD 
SINGULAR_VALUE_DECOMPOSITION_INC=\
	-I../../../external/Singular_Value_Decomposition/
# Eigen dependency
EIGEN3_INC=-I$(DEFAULT_PREFIX)/include/eigen3 -I$(DEFAULT_PREFIX)/include/eigen3/unsupported
INC+=$(EIGEN3_INC) $(SINGULAR_VALUE_DECOMPOSITION_INC)

obj: 
	mkdir -p obj

../../../lib/libiglsvd3x3.a: $(OBJ_FILES)
	rm -f $@
	ar cqs $@ $(OBJ_FILES)

obj/%.o: %.cpp %.h
	$(GG) $(AFLAGS) $(CFLAGS) -c -o $@ $< $(INC)

clean:
	rm -f obj/*.o
	rm -f ../../../lib/libiglsvd3x3.a
