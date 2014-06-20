cd ../external/AntTweakBar/src/
rm ../lib/*.a
make -f makefile.osx.igl

cd ../../glfw/
cmake -DCMAKE_BUILD_TYPE=Release .
rm src/*.a
make

cd ../../tutorial
