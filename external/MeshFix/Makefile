JMESH=./JMeshLib-1.2/
JMESH_INC=-I${JMESH}/include/
JMESH_LIB=-L${JMESH}/lib/ -ljmesh
JMESHEXT=./JMeshExt-1.0alpha_src/
JMESHEXT_INC=-I${JMESHEXT}/include/
JMESHEXT_LIB=-L${JMESHEXT}/lib -ljmeshext
NL=./OpenNL3.2.1/build/Darwin-Release
NL_LIB=-L${NL}/binaries/lib/ -lnl -L../SuperLU_4.3/lib -lsuperlu_4.3 -framework Accelerate
JMESHEXT_LIB=-L${JMESHEXT}/lib -ljmeshext

INC=${JMESH_INC} ${JMESHEXT_INC}
LIB=${JMESH_LIB} ${JMESHEXT_LIB} ${NL_LIB}

CFLAGS+=-DIS64BITPLATFORM
OPTFLAGS+=-O3

meshfix:
	g++ ${OPTFLAGS} -c meshfix.cpp -o meshfix.o ${INC} ${CFLAGS}
	g++ -o meshfix meshfix.o ${LIB}
clean:
	rm -f meshfix.o
	rm -f meshfix
