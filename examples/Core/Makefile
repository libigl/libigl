
all: example1

igl_lib=../../
eigen_lib=/opt/local/include/eigen3/

CFLAGS=-g
inc=-I$(igl_lib)/include -I$(eigen_lib)
lib=-L$(igl_lib)/lib -ligl

matlab_inc=-I$(MATLAB)/extern/include 
matlab_lib=-L$(MATLAB)/bin/maci64/ -leng -lmx

example1: example1.o
	g++ $(CFLAGS) -o example1 example1.o -framework OpenGL -framework GLUT $(lib)

example1.o: example1.cpp
	g++ $(CFLAGS) -c example1.cpp -o example1.o $(inc)


clean:
	rm -f example1.o
	rm -f example1
	rm -f bunny_out.off
