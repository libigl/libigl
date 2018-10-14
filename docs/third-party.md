The `external/` directory contains external libraries which are either difficult to obtain, difficult to compile or patched for libigl.

## AntTweakBar

Patched to support copying and pasting from text fields (see http://www.alecjacobson.com/weblog/?p=2378)

Added makefiles for Mac OS X with modern GCC (AntTweakBar/src/Makefile.osx.igl) and for building with Mesa (AntTweakBar/src/Makefile.mesa.igl)

## Embree

Install `ispc`  and `tbb` for example (`brew install ispc tbb`), then

```bash
mkdir build
cd build
cmake -DCMAKE_C_COMPILER=/usr/bin/gcc -DCMAKE_CXX_COMPILER=/usr/bin/g++ -DCMAKE_BUILD_TYPE=Release ..
make
```

## TinyXML2

double precision bug fixes:

tinyxml2.h line 1286 add function:

```cpp
void SetAttribute( const char* name, float value ) {
	XMLAttribute* a = FindOrCreateAttribute( name );
	a->SetAttribute( value );
}
```

tinyxml2.cpp line 434 replace with:

```cpp
TIXML_SNPRINTF( buffer, bufferSize, "%.15e", v );
```

## CGAL

CGAL can be built as a CMake external project thanks to the `CMakeLists.txt` provided in this folder. While this is mainly intended as a convenience for Windows users, and CI builds on AppVeyor, this should work on Linux/macOS as well. To build CGAL and Boost with the provided CMake script, build this folder as you would compile any CMake project (use CMake GUI and MSVC on Windows):

```bash
mkdir build
cd build
cmake ..
make -j4
```

Once this is done, just build the libigl tutorials, and it should properly detect CGAL and Boost.
