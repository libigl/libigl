include ../Makefile.conf

CXX=g++

LIBIGL_LIB+=-liglmosek -liglbbw -ligltetgen

ifdef IGL_NO_MOSEK
CFLAGS+=-DIGL_NO_MOSEK
else
# Adjust your mosek paths etc. accordingly
ifndef MOSEKPLATFORM
  MOSEKPLATFORM=osx64x86
endif
ifndef MOSEKVERSION
  MOSEKVERSION=6
endif
MOSEK=/usr/local/mosek
MOSEK_INC=-I$(MOSEK)/$(MOSEKVERSION)/tools/platform/$(MOSEKPLATFORM)/h
MOSEK_LIB=-L$(MOSEK)/$(MOSEKVERSION)/tools/platform/$(MOSEKPLATFORM)/bin -lmosek64
INC+=$(MOSEK_INC)
endif

INC=-I. ${EIGEN3_INC} ${LIBIGL_INC} ${TETGEN_INC} ${MOSEK_INC} 
LIB=${TETGEN_LIB} ${MOSEK_LIB} ${OPENGL_LIB} ${GLUT_LIB} ${LIBIGL_LIB} 

# Make file for bbw_demo
.PHONY: all
all: bbw_demo

CPP_FILES=$(wildcard ./*.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o))) 

bbw_demo: obj $(OBJ_FILES)
	$(CXX) -o bbw_demo $(CFLAGS) $(AFLAGS) $(OPENMP) $(OBJ_FILES) $(LIB)

obj:
	mkdir -p obj

obj/%.o: %.cpp %.h
	$(CXX) $(CFLAGS) $(AFLAGS) $(OPENMP) -o $@ -c $< $(INC)
obj/%.o: %.cpp
	$(CXX) $(CFLAGS) $(AFLAGS) $(OPENMP) -o $@ -c $< $(INC)

.PHONY: clean
clean:
	rm -f $(OBJ_FILES)
	rm -f bbw_demo
