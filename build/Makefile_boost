include ../../../Makefile.conf

.PHONY: all
all: libiglboost
debug: libiglboost

include ../../../Makefile.conf
all: CFLAGS += -O3 -DNDEBUG
debug: CFLAGS += -g -Wall 

.PHONY: libiglboost
libiglboost: obj ../../../lib/libiglboost.a

CPP_FILES=$(wildcard *.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

# include igl headers
INC+=-I../../../include/

# EXPECTS THAT CFLAGS IS ALREADY SET APPROPRIATELY 

# Eigen dependency
EIGEN3_INC=-I$(DEFAULT_PREFIX)/include/eigen3 -I$(DEFAULT_PREFIX)/include/eigen3/unsupported
INC+=$(EIGEN3_INC)

# boost dependency
BOOST=/opt/local/
BOOST_INC=-I$(BOOST)/include
INC+=$(BOOST_INC)

obj: 
	mkdir -p obj

../../../lib/libiglboost.a: $(OBJ_FILES)
	rm -f $@
	ar cqs $@ $(OBJ_FILES)

obj/%.o: %.cpp %.h
	g++ $(AFLAGS) $(OPENMP) $(CFLAGS) -c -o $@ $< $(INC)

clean:
	rm -f obj/*.o
	rm -f ../../../lib/libiglboost.a
