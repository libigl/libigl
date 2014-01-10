NL=../OpenNL3.2.1/
NL_INC=-I${NL}/src/NL -I${NL}/src/

JMESH=../JMeshLib-1.2/
JMESH_INC=-I${JMESH}/include/

INC=${JMESH_INC} ${NL_INC} -I./include
CFLAGS+=-DIS64BITPLATFORM
OPTFLAGS+=-O3

CPP_FILES=$(wildcard ./src/*.cpp)
OBJ_FILES=$(addprefix obj/,$(notdir $(CPP_FILES:.cpp=.o)))

.PHONY:
all: lib/libjmeshext.a

obj:
	mkdir -p obj

lib/libjmeshext.a: obj $(OBJ_FILES) obj/predicates.o
	rm -f $@
	ar cqs $@ $(OBJ_FILES) obj/predicates.o

obj/predicates.o: src/JRS_Predicates/jrs_predicates.c
	gcc  -O0 $(CFLAGS) -o $@ -c $< $(INC)

obj/%.o: src/%.cpp
	g++ $(OPTFLAGS) $(CFLAGS) -o $@ -c $< $(INC)

clean:
	rm -f $(OBJ_FILES) obj/predicates.o
	rm -f lib/libjmeshext.a
