include ../../../Makefile.conf
all: CFLAGS += -O3 -DNDEBUG 
debug: CFLAGS += -g -Wall -Werror

.PHONY: all
all: libiglpng
debug: libiglpng

.PHONY: libpng
libiglpng: obj ../../../lib/libiglpng.a

CPP_FILES=$(wildcard *.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

# include igl headers
INC+=-I../../../include/

# EXPECTS THAT CFLAGS IS ALREADY SET APPROPRIATELY 

# YIMG dependency
YIMG=../../../external/yimg
YIMG_INC=-I$(YIMG)
INC+=$(YIMG_INC)
YIMG_STATIC_LIB=$(YIMG)/libyimg.a

# BOOST dependency
INC+=-I/opt/local/include/

obj: 
	mkdir -p obj

../../../lib/libiglpng.a: $(OBJ_FILES)
	rm -f $@
	ar cqs $@ $(OBJ_FILES)

obj/%.o: %.cpp %.h
	g++ $(AFLAGS) $(CFLAGS) -c -o $@ $< $(INC)

clean:
	rm -f obj/*.o
	rm -f ../../../lib/libiglpng.a
