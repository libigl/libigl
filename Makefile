.PHONY: all
all: lib examples extras

# Shared flags etc.
include Makefile.conf

# optimized default settings
all: LFLAGS +=
all: CFLAGS += -O3 -DNDEBUG -j 
debug: CFLAGS += -g -Wall -Werror

EXTRA_DIRS=
ifeq ($(IGL_WITH_TETGEN),1)
	# append tetgen extra dir
	EXTRA_DIRS+=include/igl/tetgen
endif

.PHONY: examples
.PHONY: extras
debug: lib
lib: obj lib/libigl.a
examples:
	make -C examples
extras:
	for p in  $(EXTRA_DIRS); \
	do \
	echo "cd $$p" ; \
	$(MAKE) -C $$p ; \
	done


#############################################################################
# SOURCE 
#############################################################################
CPP_FILES=$(wildcard include/igl/*.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

# include igl headers
INC+=-Iinclude/

#############################################################################
# DEPENDENCIES
#############################################################################

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
	g++ $(CFLAGS) $(AFLAGS) -c -o $@ $< $(INC)

clean:
	rm -f obj/*.o
	rm -f lib/libigl.a
	make -C examples clean
	for p in  $(EXTRA_DIRS); \
	do \
	echo "cd $$p" ; \
	$(MAKE) -C $$p clean; \
	done
