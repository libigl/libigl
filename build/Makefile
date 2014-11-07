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

EXTRAS=
ifeq ($(IGL_WITH_BBW),1)
	EXTRAS += bbw
endif
ifeq ($(IGL_WITH_BOOLEAN),1)
	EXTRAS += boolean
endif
ifeq ($(IGL_WITH_BOOST),1)
	EXTRAS += boost
endif
ifeq ($(IGL_WITH_CGAL),1)
	EXTRAS += cgal
endif
ifeq ($(IGL_WITH_EMBREE),1)
	EXTRAS += embree
endif
ifeq ($(IGL_WITH_COMISO),1)
	EXTRAS += comiso
endif
ifeq ($(IGL_WITH_MATLAB),1)
	EXTRAS += matlab
endif
ifeq ($(IGL_WITH_MOSEK),1)
	EXTRAS += mosek
endif
ifeq ($(IGL_WITH_PNG),1)
	EXTRAS += png
endif
ifeq ($(IGL_WITH_SVD3X3),1)
	EXTRAS += svd3x3
endif
ifeq ($(IGL_WITH_TETGEN),1)
	# append tetgen extra dir
	EXTRAS += tetgen
endif
ifeq ($(IGL_WITH_VIEWER),1)
	EXTRAS += viewer
endif
ifeq ($(IGL_WITH_XML),1)
	EXTRAS += xml
endif

.PHONY: examples
.PHONY: extras
debug: lib extras
lib: ../lib/libigl.a
examples: lib extras
	make -C ../examples
extras:
	for p in  $(EXTRAS); \
	do \
	$(MAKE) -f Makefile_$$p $(DEBUG); \
	done


#############################################################################
# SOURCE
#############################################################################
CPP_FILES=$(wildcard ../include/igl/*.cpp)
H_FILES=$(wildcard ../include/igl/*.h)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

# include igl headers
INC+=-I../include/

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
ANTTWEAKBAR_INC=-I../external/AntTweakBar/include
INC+=$(ANTTWEAKBAR_INC)

## OpenGL dependency
#LIB+=-framework OpenGL
#LIB+=-framework GLUT
#LIB+=-framework AppKit

.PHONY: obj
obj:
	mkdir -p obj

../lib/libigl.a: obj $(OBJ_FILES)
	mkdir -p ../lib
	rm -f $@
	ar cqs $@ $(OBJ_FILES)

obj/%.o: ../include/igl/%.cpp ../include/igl/%.h
	$(GG) $(CFLAGS) $(AFLAGS) -c -o $@ $< $(INC)

#lib/igl.framework/:
#	mkdir -p $@
#	cp lib/*.a $@
#	mv $@/libigl.a $@/igl
#	mkdir -p $@/Libraries
#	mv $@/*.a $@/Libraries
#	mkdir -p $@/Headers
#	cp $(H_FILES) $@/Headers
#	for p in $(EXTRAS); \
#	do \
#	mkdir $@/Headers/$$p; \
#	cp include/igl/$$p/*.h $@/Headers/$$p; \
#	done


clean:
	rm -rf ../lib/igl.framework/
	rm -f $(OBJ_FILES)
	rm -f ../lib/libigl.a
	make -C ../examples clean
	for p in  $(EXTRA_DIRS); \
	do \
	echo "cd $$p" ; \
	$(MAKE) -C $$p clean; \
	done
