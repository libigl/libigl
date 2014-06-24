
.PHONY: all
all: libigltetgen
debug: libigltetgen

include ../../../Makefile.conf
all: CFLAGS += -O3 -DNDEBUG 
debug: CFLAGS += -g -Wall -Werror

.PHONY: libtetgen
libigltetgen: obj ../../../lib/libigltetgen.a

CPP_FILES=$(wildcard *.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

# include igl headers
INC+=-I../../../include/

# EXPECTS THAT CFLAGS IS ALREADY SET APPROPRIATELY 

# Eigen dependency
EIGEN3_INC=-I$(DEFAULT_PREFIX)/include/eigen3 -I$(DEFAULT_PREFIX)/include/eigen3/unsupported
INC+=$(EIGEN3_INC)

# Tetgen dependency
TETGEN=../../../external/tetgen
TETGEN_INC=-I$(TETGEN)
INC+=$(TETGEN_INC)
TETGEN_STATIC_LIB=$(TETGEN)/libtet.a

obj: 
	mkdir -p obj

../../../lib/libigltetgen.a: $(OBJ_FILES)
	rm -f $@
	ar cqs $@ $(OBJ_FILES)

obj/%.o: %.cpp %.h
	g++ $(AFLAGS) $(CFLAGS) -c -o $@ $< $(INC)

clean:
	rm -f obj/*.o
	rm -f ../../../lib/libigltetgen.a
