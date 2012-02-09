all: lib examples

.PHONY: all
.PHONY: examples

debug: lib

lib: obj libigl.a

examples:
	make -C examples

CPP_FILES=$(wildcard ./*.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

# optimized default settings
all: LFLAGS +=
all: CFLAGS += -O3 -DNDEBUG -Wall
debug: CFLAGS += -g -Wall -Werror

# Eigen dependency
EIGEN3_INC=-I/usr/local/include/eigen3 -I/usr/local/include/eigen3/unsupported
INC+=$(EIGEN3_INC)

# AntTweakBar dependency
MY_INC=-I/opt/local/include -I/opt/local/include/eigen3 -I/opt/local/include/eigen3/unsupported
INC+=$(MY_INC)


## OpenGL dependency
#LIB+=-framework OpenGL
#LIB+=-framework GLUT
#LIB+=-framework AppKit

obj: 
	mkdir -p obj

libigl.a: $(OBJ_FILES)
	rm -f $@
	ar cqs $@ $(OBJ_FILES)

obj/%.o: %.cpp %.h
	g++ $(CFLAGS) -c -o $@ $< $(INC)

clean:
	rm -f obj/*.o
	rm -f libigl.a
	make -C examples clean
