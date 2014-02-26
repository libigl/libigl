#!/bin/bash
git clone git@github.com:libigl/libigl.git
cd libigl/
make -C external/AntTweakBar/src -f Makefile.osx.igl
make -C external/yimg
cd external/tetgen
make clean
make -f Makefile.igl tetlib
rm -f obj/*.o
rm tetgen
make -f Makefile.igl tetgen
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
