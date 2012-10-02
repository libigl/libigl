LIBS  = -framework OpenGL -framework GLUT
CFLAGS = -Wall -O3

C_FILES=$(wildcard *.c)
OBJ_FILES=$(notdir $(C_FILES:.c=.o))

%.o : %.c
	gcc -c $(CFLAGS) $< -o $@ -I.

.PHONY: libmesh
libmesh: libmesh.a

libmesh.a: $(OBJ_FILES) 
	ar cqs $@ $(OBJ_FILES)
