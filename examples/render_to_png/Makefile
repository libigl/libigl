
.PHONY: all

# Shared flags etc.
include ../../Makefile.conf

# YIMG dependency
YIMG=../../external/yimg
YIMG_LIB=-L$(YIMG) -lyimg -lz -L/usr/X11/lib -lpng -bind_at_load
YIMG_INC=-I/usr/X11/include -I$(YIMG) 

EIGEN3_INC=-I/opt/local/include/eigen3 -I/opt/local/include/eigen3/unsupported

# LIBPNG dependency

all: example

.PHONY: example

IGL=../../

inc=-DIGL_HEADER_ONLY -I$(IGL)/include -I$(DEFAULT_PREFIX)/include $(YIMG_INC) $(EIGEN3_INC)
lib=-L$(IGL)/lib -ligl -liglpng -L$(DEFAULT_PREFIX)/lib $(OPENGL_LIB) $(GLUT_LIB) $(YIMG_LIB) $(LIBPNG_LIB)

example: example.o
	g++ $(CFLAGS) -o example example.o $(lib)

example.o: example.cpp
	g++ $(CFLAGS) -c example.cpp -o example.o $(inc)
clean:
	rm -f example.o
	rm -f example
