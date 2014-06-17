.PHONY: all
all: lib extras examples
framework: lib extras lib/igl.framework/

# Shared flags etc.
include Makefile.conf
$(info Hello, $(IGL_USERNAME)!)

# optimized default settings
all: LFLAGS +=
OPTFLAGS+=-O3 -DNDEBUG $(OPENMP)
#debug: OPTFLAGS= -g -Wall -Werror
debug: OPTFLAGS= -g -Wall
debug: DEBUG=debug
CFLAGS += $(OPTFLAGS)
#CFLAGS += -DIGL_NO_OPENGL -DIGL_NO_ANTTWEAKBAR
# We use well-supported features of c++11

EXTRA_DIRS=
ifeq ($(IGL_WITH_BBW),1)
	EXTRA_DIRS+=include/igl/bbw
	EXTRAS += bbw
endif
ifeq ($(IGL_WITH_BOOST),1)
	EXTRA_DIRS+=include/igl/boost
	EXTRAS += boost
endif
ifeq ($(IGL_WITH_CGAL),1)
	EXTRA_DIRS+=include/igl/cgal
	EXTRAS += cgal
endif
ifeq ($(IGL_WITH_EMBREE),1)
	EXTRA_DIRS+=include/igl/embree
	EXTRAS += embree
endif
ifeq ($(IGL_WITH_MATLAB),1)
	EXTRA_DIRS+=include/igl/matlab
	EXTRAS += matlab
endif
ifeq ($(IGL_WITH_MOSEK),1)
	EXTRA_DIRS+=include/igl/mosek
	EXTRAS += mosek
endif
ifeq ($(IGL_WITH_PNG),1)
	EXTRA_DIRS+=include/igl/png
	EXTRAS += png
endif
ifeq ($(IGL_WITH_SVD3X3),1)
	EXTRA_DIRS+=include/igl/svd3x3
	EXTRAS += svd3x3
endif
ifeq ($(IGL_WITH_TETGEN),1)
	# append tetgen extra dir
	EXTRA_DIRS+=include/igl/tetgen
	EXTRAS += tetgen
endif
ifeq ($(IGL_WITH_VIEWER),1)
	EXTRA_DIRS+=include/igl/viewer
	EXTRAS += viewer
endif
ifeq ($(IGL_WITH_XML),1)
	EXTRA_DIRS+=include/igl/xml
	EXTRAS += xml
endif

.PHONY: examples
.PHONY: extras
debug: lib extras
lib: lib/libigl.a
examples: lib extras
	make -C examples
extras:
	for p in  $(EXTRA_DIRS); \
	do \
	echo "cd $$p" ; \
	$(MAKE) -C $$p $(DEBUG); \
	done


#############################################################################
# SOURCE 
#############################################################################
CPP_FILES=$(wildcard include/igl/*.cpp)
H_FILES=$(wildcard include/igl/*.h)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

# include igl headers
INC+=-Iinclude/

#############################################################################
# DEPENDENCIES
#############################################################################
INC+=$(OPENGL_INC)

# Eigen dependency
ifndef EIGEN3_INC
	EIGEN3_INC=-I$(DEFAULT_PREFIX)/include/eigen3 -I$(DEFAULT_PREFIX)/include/eigen3/unsupported
endif
INC+=$(EIGEN3_INC)

# AntTweakBar dependency
#ANTTWEAKBAR_INC=-I$(DEFAULT_PREFIX)/include
ANTTWEAKBAR_INC=-Iexternal/AntTweakBar/include
INC+=$(ANTTWEAKBAR_INC)

## OpenGL dependency
#LIB+=-framework OpenGL
#LIB+=-framework GLUT
#LIB+=-framework AppKit

.PHONY: obj
obj: 
	mkdir -p obj

lib/libigl.a: obj $(OBJ_FILES)
	mkdir -p lib
	rm -f $@
	ar cqs $@ $(OBJ_FILES)

obj/%.o: include/igl/%.cpp include/igl/%.h
	$(GG) $(CFLAGS) $(AFLAGS) -c -o $@ $< $(INC)

lib/igl.framework/:
	mkdir -p $@
	cp lib/*.a $@
	mv $@/libigl.a $@/igl
	mkdir -p $@/Libraries
	mv $@/*.a $@/Libraries
	mkdir -p $@/Headers
	cp $(H_FILES) $@/Headers
	for p in $(EXTRAS); \
	do \
	mkdir $@/Headers/$$p; \
	cp include/igl/$$p/*.h $@/Headers/$$p; \
	done


clean:
	rm -rf lib/igl.framework/
	rm -f obj/*.o
	rm -f lib/libigl.a
	make -C examples clean
	for p in  $(EXTRA_DIRS); \
	do \
	echo "cd $$p" ; \
	$(MAKE) -C $$p clean; \
	done
