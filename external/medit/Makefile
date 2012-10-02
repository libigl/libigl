LIBMESH = ./libmesh/
INC = -I${LIBMESH}
LIBS  = -framework OpenGL -framework GLUT -L${LIBMESH} -lmesh
CFLAGS = -Wall -O3

C_FILES=$(wildcard *.c)
OBJ_FILES=$(notdir $(C_FILES:.c=.o))

%.o : %.c
	gcc -c $(CFLAGS) $< -o $@ $(INC)

medit: $(OBJ_FILES)
	gcc -o $@ $^ $(CFLAGS) $(LIBS)
