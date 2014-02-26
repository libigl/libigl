libigl - A simple c++ geometry processing library
=================================================

<http://igl.ethz.ch/projects/libigl/>
<https://github.com/alecjacobson/libigl/>

Copyright 2013 - Alec Jacobson, Daniele Panozzo, Olga Diamanti, Kenshi
Takayama, Leo Sacht 

This is first and foremost a *header* library. Each header file should contain
a single function.  The function may have multiple prototypes. All functions
should use the igl namespace and should adhere to the conventions and styles
listed below. 

## Dependencies ##
- Eigen3  Last tested with Eigen Version 3.2

### Optional ###
- OpenGL (`IGL_NO_OPENGL`)
- AntTweakBar  (`IGL_NO_ANTTWEAKBAR`) Last tested 1.16 (see
-   libigl/external/AntTweakBar)
- GLEW  Windows only
- OpenMP  
- libpng  libiglpng extra only
- Mosek  libiglmosek extra only
- Matlab  libiglmatlab extra only
- boost  libboost extra only

### Optional (included in external/) ###
- TetGen  libigltetgen extra only
- Embree  libiglembree extra only
- tinyxml2  libiglxml extra only
 
## Header only ##
libigl is designed to work "out-of-the-box" as a headers only library. To
include libigl in your project. You need only include the libigl/include/
directory in your include path and define the `IGL_HEADER_ONLY` macro. To 
compile a hello-word example.cpp:

    #include <Eigen/Dense>
    #include <igl/readOBJ.h>
    #include <iostream>
    int main(int argc, char * argv[])
    {
      if(argc>1)
      {
        Eigen::MatrixXd V;
        Eigen::MatrixXi F;
        igl::readOBJ(argv[1],V,F);
        std::cout<<"Hello, mesh with "<<V.rows()<<" vertices!"<<std::endl;
      }else{
        std::cout<<"Hello, world!"<<std::endl;
      }
      return 0;
    }

using gcc (replacing appropriate paths):

    g++ -DIGL_HEADER_ONLY -I/usr/local/igl/libigl/include \
      -I/opt/local/include/eigen3 example.cpp -o example

Then run this example with:

    ./example examples/shared/TinyTorus.obj

## Compile ##
libigl is developed most often on Mac OS X, though has current users in Linux and Windows.

### Linux/Mac OS X/Cygwin ###
  
libigl may also be compiled to a static library. This is advantageous when
building a project with libigl, since the header only directive can slow down
compile times.

To build the entire libigl library producing lib/libigl.a, issue:
  
    make lib
  
You may need to edit Makefile.conf accordingly. Best to give yourself an
`IGL_USERNAME` and add a custom install suite for yourself. Then you can enable
appropriate extras.
  
#### Extras ####
Once you've set up an `IGL_USERNAME` and enabled extras within Makefile.conf.
You can build the extra libraries (into lib/ligiglpng.a, lib/libiglmatlab.a,
lib/libigltetgen.a, lib/libiglmosek.a, etc.) by issuing:
  
    make extras
  
#### Examples ####
You can make a slew of examples by issuing:
  
    make examples
  
#### External ####
Finally there are a number of external libraries that we include in
./external/ because they are either difficult to obtain or they have been
patched for easier use with libigl. Please see the respective readmes in
those directories.


##### Installing AntTweakBar #####
To build the a static AntTweakBar library on Mac OS X issue:

    cd external/AntTweakBar/src
    make -f Makefile.osx.igl

##### Installing Tetgen #####
To build the tetgen library and executable on Mac OS X issue:

    cd external/tetgen
    make clean
    rm -f obj/*.o
    make -f Makefile.igl tetlib
    rm -f obj/*.o
    rm tetgen
    make -f Makefile.igl tetgen

##### Installing Embree 2.0 #####
To build the embree library and executables on Mac OS X issue:

    cd external/embree
    mkdir build
    cd build
    cmake ..
    # Or using a different compiler
    #cmake .. -DCMAKE_C_COMPILER=/opt/local/bin/gcc -DCMAKE_CXX_COMPILER=/opt/local/bin/g++
    make
    # Could also install embree to your root, but libigl examples don't expect
    # this
    #sudo make install

##### Installing tinyxml2 #####
To build the a static tinyxml2 library on Mac OS X issue:

    cd external/tinyxml2
    cmake .
    make


##### Installing YImg #####
To build the a static YImg library on Mac OS X issue:

    cd external/yimg
    make


### Windows (Experimental) ###
To build a static library (.lib) on windows, open Visual Studio 2010.

- New > Project ...
- Visual C++ > Win32
- Win32 Console Application
- Name: libiglVisualStudio
- Uncheck "Create directory for solution"
- Then hit OK, and then Next
- Check "Static Library"
- Uncheck "Precompiled headers"
- Add all include/igl/*.cpp to the sources directory
- Add all include/igl/*.h to the headers directory
- Open Project > libigl Properties...
- Add the path to eigen3 to the include paths
- Change the target name to libigl
- Build and pray (this should create libigl.lib

[Source](http://msdn.microsoft.com/en-us/library/ms235627(v=vs.80).aspx)

## Examples ##
To get started, we advise that you take a look at a few examples:

    ./examples/hello-world/

    ./examples/meshio/

    ./examples/basic-topology/

    ./examples/ReAntTweakBar/

## Development ##
Further documentation for developers is listed in tutorial.html,
style_guidelines.html

## License ##
See `LICENSE.txt`

## Zipping ##
Zip this directory without .git litter and binaries using:

    git archive —prefix=libigl/ -o libigl.zip master

## Contact ##
libigl is a group endeavor led by Alec Jacobson and Daniele Panozzo. Please
contact [alecjacobson@gmail.com](mailto:alecjacobson@gmail.com) if you have
questions or comments. We are happy to get feedback! Enjoy!
