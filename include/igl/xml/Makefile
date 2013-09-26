include ../../../Makefile.conf

.PHONY: all
all: libiglxml
debug: libiglxml

include ../../../Makefile.conf
all: CFLAGS += -O3 -DNDEBUG
debug: CFLAGS += -g -Wall

.PHONY: libiglxml
libiglxml: obj ../../../lib/libiglxml.a

CPP_FILES=$(wildcard *.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

# include igl headers
INC+=-I../../../include/

# EXPECTS THAT CFLAGS IS ALREADY SET APPROPRIATELY 

# Eigen dependency
EIGEN3_INC=-I$(DEFAULT_PREFIX)/include/eigen3 -I$(DEFAULT_PREFIX)/include/eigen3/unsupported
INC+=$(EIGEN3_INC)

#AntTweakbar dependency
ANTTWEAKBAR_INC=-I../../../external/AntTweakBar/include

# xml dependency
# TODO: linux, 32 bit etc
TINYXML2=../../../external/tinyxml2
TINYXML2_INC=-I$(TINYXML2)
#TINYXML2_LIB=-L$(TINYXML2) -ltinyxml2
INC+=$(TINYXML2_INC) $(ANTTWEAKBAR_INC)

# AntTweakBar dependency
ANTTWEAKBAR_INC=-I../../../external/AntTweakBar/include
INC+=$(ANTTWEAKBAR_INC)

#CFLAGS+=-std=c++11

obj: 
	mkdir -p obj

../../../lib/libiglxml.a: $(OBJ_FILES)
	rm -f $@
	ar cqs $@ $(OBJ_FILES)

obj/%.o: %.cpp %.h
	g++ $(AFLAGS) $(CFLAGS) -c -o $@ $< $(INC)

clean:
	rm -f obj/*.o
	rm -f ../../../lib/libiglxml.a
