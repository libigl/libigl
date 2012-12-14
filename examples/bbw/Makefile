include ../../Makefile.conf

IGL=../../
IGL_INC=-I${IGL}/include
IGL_LIB=-L${IGL}/lib -ligl

ifndef MOSEKPLATFORM
  MOSEKPLATFORM=osx64x86
endif
MOSEK=/usr/local/mosek
MOSEK_INC=-I$(MOSEK)/6/tools/platform/$(MOSEKPLATFORM)/h
MOSEK_LIB=-L$(MOSEK)/6/tools/platform/$(MOSEKPLATFORM)/bin -lmosek64 -liglmosek

EIGEN3_INC=-I$(DEFAULT_PREFIX)/include/eigen3 -I$(DEFAULT_PREFIX)/include/eigen3/unsupported

TETGEN=$(IGL)/external/tetgen
TETGEN_LIB=-L$(TETGEN) -ligltetgen -ltet 
TETGEN_INC=-I$(TETGEN)

INC=-I. ${EIGEN3_INC} ${IGL_INC} ${TETGEN_INC} ${MOSEK_INC} 
LIB=${TETGEN_LIB} ${MOSEK_LIB} ${OPENGL_LIB} ${GLUT_LIB} ${IGL_LIB} 

# Make file for bbw_demo
.PHONY: all
all: bbw_demo

CPP_FILES=$(wildcard ./*.cpp)
#C_FILES=$(wildcard ./*.c)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o))) 
#$(addprefix obj/,$(notdir $(C_FILES:.c=.o)))
CFLAGS+=-DSHADER_DIR="\"./GLSL\""

bbw_demo: obj $(OBJ_FILES)
	g++ -o bbw_demo $(CFLAGS) $(OBJ_FILES) $(LIB)

obj:
	mkdir -p obj

#obj/%.o: %.c %.h
#	gcc $(CFLAGS) -o $@ -c $< $(INC)
obj/%.o: %.cpp %.h
	g++ $(CFLAGS) -o $@ -c $< $(INC)
obj/%.o: %.cpp
	g++ $(CFLAGS) -o $@ -c $< $(INC)

.PHONY: clean
clean:
	rm -f $(OBJ_FILES)
	rm -f bbw_demo
