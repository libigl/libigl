This directory contains external libraries which are either difficult to
obtain, difficult to compile or patched for libigl.

## AntTweakBar
Patched to support copying and pasting from text fields (see
http://www.alecjacobson.com/weblog/?p=2378)

Added makefiles for Mac OS X with modern GCC (AntTweakBar/src/Makefile.osx.igl)
and for building with Mesa (AntTweakBar/src/Makefile.mesa.igl)


## Embree

Install `ispc`  and `tbb` for example (`brew install ispc tbb`),
then

    mkdir build
    cd build
    cmake -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DCMAKE_BUILD_TYPE=Release ..
    make

## TinyXML2
double precision bug fixes:

tinyxml2.h line 1286 add function:

void SetAttribute( const char* name, float value ) {
	XMLAttribute* a = FindOrCreateAttribute( name );
	a->SetAttribute( value );
}

tinyxml2.cpp line 434 replace with:

TIXML_SNPRINTF( buffer, bufferSize, "%.15e", v );

