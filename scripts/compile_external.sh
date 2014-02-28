#!/bin/bash
git clone git@github.com:libigl/libigl.git
cd libigl/
make -C external/AntTweakBar/src -f Makefile.osx.igl
make -C external/yimg
make -C external/medit/libmesh
make -C external/medit/ -f Makefile.igl medit
cd external/tetgen
make clean
mkdir obj
make -f Makefile.igl tetgen
rm -f obj/*.o
make -f Makefile.igl tetlib
mkdir -p ../embree/build
cd ../embree/build
cmake .. -DCMAKE_C_COMPILER=/opt/local/bin/gcc -DCMAKE_CXX_COMPILER=/opt/local/bin/g++
make
cd ../../tinyxml2
cmake .
make
cd ../../
make -j12
make -j12
