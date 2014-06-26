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
cmake -DCMAKE_BUILD_TYPE=Release ../
make

cd ../../../tutorial
