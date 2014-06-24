
.PHONY: all
all: libiglcgal
debug: libiglcgal

include ../../../Makefile.conf
all: CFLAGS += -O3 -DNDEBUG 
debug: CFLAGS += -g -Wall -Werror

.PHONY: libcgal
libiglcgal: obj ../../../lib/libiglcgal.a

CPP_FILES=$(wildcard *.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

# include igl headers
INC+=-I../../../include/

# EXPECTS THAT CFLAGS IS ALREADY SET APPROPRIATELY 

# Eigen dependency
EIGEN3_INC=-I$(DEFAULT_PREFIX)/include/eigen3 -I$(DEFAULT_PREFIX)/include/eigen3/unsupported
INC+=$(EIGEN3_INC)

# CGAL dependency
CGAL=$(DEFAULT_PREFIX)
CGAL_INC=-I$(CGAL)/include
CGAL_FLAGS=-frounding-math -fsignaling-nans 
CFLAGS+=$(CGAL_FLAGS)
INC+=$(CGAL_INC)

obj: 
	mkdir -p obj

../../../lib/libiglcgal.a: $(OBJ_FILES)
	rm -f $@
	ar cqs $@ $(OBJ_FILES)

obj/%.o: %.cpp %.h
	g++ $(AFLAGS) $(CFLAGS) -c -o $@ $< $(INC)

clean:
	rm -f obj/*.o
	rm -f ../../../lib/libiglcgal.a
