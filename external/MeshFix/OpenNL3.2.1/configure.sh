#!/bin/sh

# This file for Linux users, 
# launches CMake and creates configuration for
# Release and Debug modes.

echo ============= Checking for cmake ============
if (cmake --version)
then
   echo "found CMake"
else
   echo "cmake not found, please install it (see http://www.cmake.org/)" 
   exit
fi

os=`uname`

for config in Release Debug
do
   echo
   echo ============= Creating makefiles for $config mode ============
   mkdir -p build/$os-$config
   (cd build/$os-$config; cmake ../../ -DCMAKE_BUILD_TYPE:STRING=$config)
done

echo
echo ============== OpenNL build configured ==================


cat << EOF
to build the libBL library :
  go to build/$os-Release or build/$os-Debug
  and 'make'

Note: local configuration can be specified in CMakeOptions.txt
(see CMakeOptions.txt.sample for an example)
You'll need to re-run configure.sh if you create or modify CMakeOptions.txt

EOF

