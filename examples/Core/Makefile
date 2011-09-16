
all: example1 example2

STATIC_ANTTWEAKBAR=""

igl_lib=../../
eigen_lib=/opt/local/include/eigen3/

ifdef STATIC_ANTTWEAKBAR
CFLAGS=-g -DSTATIC_ANTTWEAKBAR
inc=-I$(igl_lib) -I./static -I$(eigen_lib)
lib= -lAntTweakBar -framework AppKit
else
CFLAGS=-g
inc=-I$(igl_lib) -I$(eigen_lib)
lib=-lAntTweakBar
endif

matlablibs=-I/Applications/MATLAB_R2010b.app/extern/include -L/Applications/MATLAB_R2010b.app/bin/maci64/ -leng -lmx

example1: example1.o
	g++ $(CFLAGS) -o example1 example1.o -framework OpenGL -framework GLUT $(lib)
	rm example1.o

example1.o: example1.cpp
	g++ $(CFLAGS) -c example1.cpp -o example1.o $(inc)

example2: example2.o
	g++ $(CFLAGS) $(matlablibs) -o example2 example2.o -framework OpenGL -framework GLUT $(lib)
	rm example2.o

example2.o: example2.cpp
	g++ $(CFLAGS) $(matlablibs) -c example2.cpp -o example2.o $(inc)

clean:
	rm -f example1.o
	rm -f example1
	rm -f example2.o
	rm -f example2
	rm -f bunny_out.off
