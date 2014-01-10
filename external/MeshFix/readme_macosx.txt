Short notes on porting meshfix to mac.

Follow readme.txt placing downloaded OpenNL and JMeshLib directories in this
directory [same as meshfix.cpp]

When compiling JMeshLib be sure in makeconf file to change the lines:

...
# On 64-bit machines you need to uncomment the following line
# -DIS64BITPLATFORM

MOREFLAGS = $(OPTM) $(STRICTALIAS)
...

to:

...
MOREFLAGS = $(OPTM) $(STRICTALIAS) -DIS64BITPLATFORM
...

When compiling JMeshExt I separated the predicates.cxx file from tetgen's
source to make jrs_predicates.h and jrs_predicates.c. Placing them according to
readme.txt. One note is that if jrs_predicates.c is to be a .c file rather than
a .cpp file then you need to extern "C"{...} encapsulate the jrs_predicates.h
header file or you will get Undefined symbol linker errors.

Finally in JMeshExt and meshfix.cpp there are several casts of pointers to
ints. This won't compile on 64-bit machines. This similar operation was fixed
in JMeshLib by the -DIS64BITPLATFORM. Thus, I've fixed the following files to
use the j_voidint typedef defined by IS64BITPLATFORM in j_mesh.h:
  JMeshExt-1.0alpha_src/src/holeFilling.cpp
  meshfix.cpp
