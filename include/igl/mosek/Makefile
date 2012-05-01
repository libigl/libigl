include ../../../Makefile.conf

.PHONY: all
all: libiglmosek
debug: libiglmosek

include ../../../Makefile.conf
all: CFLAGS += -O3 -DNDEBUG -j 
debug: CFLAGS += -g -Wall -Werror

.PHONY: libiglmosek
libiglmosek: obj ../../../lib/libiglmosek.a

CPP_FILES=$(wildcard *.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

# include igl headers
INC+=-I../../../include/

# EXPECTS THAT CFLAGS IS ALREADY SET APPROPRIATELY 

# Eigen dependency
EIGEN3_INC=-I/opt/local/include/eigen3 -I/opt/local/include/eigen3/unsupported
INC+=$(EIGEN3_INC)

# mosek dependency
# TODO: linux, 32 bit etc
MOSEKPLATFORM=osx64x86
MOSEK=/usr/local/mosek
MOSEK_INC=-I$(MOSEK)/6/tools/platform/$(MOSEKPLATFORM)/h
MOSEK_LIB=-L$(MOSEK)/6/tools/platform/$(MOSEKPLATFORM)/bin -lmosek64
INC+=$(MOSEK_INC)

obj: 
	mkdir -p obj

../../../lib/libiglmosek.a: $(OBJ_FILES)
	rm -f $@
	ar cqs $@ $(OBJ_FILES)

obj/%.o: %.cpp %.h
	g++ $(AFLAGS) $(CFLAGS) -c -o $@ $< $(INC)

clean:
	rm -f obj/*.o
	rm -f ../../../lib/libiglmosek.a
