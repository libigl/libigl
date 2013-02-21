
.PHONY: all

# Shared flags etc.
include ../../Makefile.conf

all: example

.PHONY: example

IGL=../../

ANTTWEAKBAR_INC=-I$(IGL)/external/AntTweakBar/include
ANTTWEAKBAR_LIB=-L$(IGL)/external/AntTweakBar/lib -lAntTweakBar -framework AppKit
inc=-DIGL_HEADER_ONLY -I$(IGL)/include -I$(DEFAULT_PREFIX)/include $(ANTTWEAKBAR_INC)
lib=-L$(IGL)/lib -ligl -L$(DEFAULT_PREFIX)/lib $(OPENGL_LIB) $(GLUT_LIB) $(ANTTWEAKBAR_LIB)

example: example.o
	g++ $(CFLAGS) -o example example.o $(lib)

example.o: example.cpp
	g++ $(CFLAGS) -c example.cpp -o example.o $(inc)
clean:
	rm -f example.o
	rm -f example
