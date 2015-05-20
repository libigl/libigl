cd ../external/AntTweakBar/src/
rm ../lib/*.a
make -f makefile.osx.igl

cd ../../glfw/
cmake -DCMAKE_BUILD_TYPE=Release .
rm src/*.a
make

cd ../embree/
rm -fr build
mkdir build
cd build
cp ../ispc/ispc-v1.8.1-osx ispc
cmake -DCMAKE_BUILD_TYPE=Release ../
make -j 16

cd ../../../tutorial
