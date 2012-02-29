all: lib examples

.PHONY: all
.PHONY: examples

debug: lib

lib: obj lib/libigl.a

examples:
	make -C examples

CPP_FILES=$(wildcard include/igl/*.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

# include igl headers
INC+=-Iinclude/

CFLAGS += -Wall -arch i386 -arch x86_64
# optimized default settings
all: LFLAGS +=
all: CFLAGS += -O3 -DNDEBUG -j 
debug: CFLAGS += -g -Wall -Werror

# Eigen dependency
EIGEN3_INC=-I/opt/local/include/eigen3 -I/opt/local/include/eigen3/unsupported
INC+=$(EIGEN3_INC)

# AntTweakBar dependency
ANTTWEAKBAR_INC=-I/opt/local/include
INC+=$(ANTTWEAKBAR_INC)


## OpenGL dependency
#LIB+=-framework OpenGL
#LIB+=-framework GLUT
#LIB+=-framework AppKit

obj: 
	mkdir -p obj

lib/libigl.a: $(OBJ_FILES)
	rm -f $@
	ar cqs $@ $(OBJ_FILES)

obj/%.o: include/igl/%.cpp include/igl/%.h
	g++ $(CFLAGS) -c -o $@ $< $(INC)

clean:
	rm -f obj/*.o
	rm -f lib/libigl.a
	make -C examples clean
