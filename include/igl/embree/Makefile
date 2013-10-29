include ../../../Makefile.conf

.PHONY: all
all: libiglembree
debug: libiglembree

include ../../../Makefile.conf
all: CFLAGS += -O3 -DNDEBUG -std=c++11
debug: CFLAGS += -g -Wall -std=c++11

.PHONY: libiglembree
libiglembree: obj ../../../lib/libiglembree.a

CPP_FILES=$(wildcard *.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

# include igl headers
INC+=-I../../../include/

# EXPECTS THAT CFLAGS IS ALREADY SET APPROPRIATELY 

# Eigen dependency
EIGEN3_INC=-I$(DEFAULT_PREFIX)/include/eigen3 -I$(DEFAULT_PREFIX)/include/eigen3/unsupported
INC+=$(EIGEN3_INC)

# embree dependency
#EMBREE=../../../external/embree
#EMBREE_INC=-I$(EMBREE)/common -I$(EMBREE)/rtcore/
EMBREE=../../../external/embree/
EMBREE_INC=-I$(EMBREE)/ -I$(EMBREE)/embree
#EMBREE_LIB=-L$(EMBREE)/bin -lrtcore -lsys
INC+=$(EMBREE_INC)

obj: 
	mkdir -p obj

../../../lib/libiglembree.a: $(OBJ_FILES)
	rm -f $@
	ar cqs $@ $(OBJ_FILES)

obj/%.o: %.cpp %.h
	g++ $(AFLAGS) $(OPENMP) $(CFLAGS) -c -o $@ $< $(INC)

clean:
	rm -f obj/*.o
	rm -f ../../../lib/libiglembree.a
