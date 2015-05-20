sudo apt-get install git
sudo apt-get install build-essential
sudo apt-get install libeigen3-dev
sudo apt-get install cmake
sudo apt-get install libx11-dev
sudo apt-get install mesa-common-dev libgl1-mesa-dev libglu1-mesa-dev
sudo apt-get install libxrandr-dev
sudo apt-get install libxi-dev
sudo apt-get install freeglut3-dev
sudo apt-get install libxmu-dev
sudo apt-get install libblas-dev libsuitesparse-dev

cd ../external/AntTweakBar/src/
rm ../lib/*.a
make

cd ../../glfw/
cmake -DCMAKE_BUILD_TYPE=Release .
rm src/*.a
make

cd ../embree/
rm -fr build
mkdir build
cd build
cp ../ispc/ispc-v1.8.1-linux ispc
cmake -DCMAKE_BUILD_TYPE=Release ../
make

cd ../../../tutorial
