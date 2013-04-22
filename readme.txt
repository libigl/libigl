libigl - A simple c++ geometry processing library

Copyright 2013 - Alec Jacobson, Daniele Panozzo, Olga Diamanti, Kenshi
Takayama, Leo Sacht, Interactive Geometry Lab - ETH Zurich

This is first and foremost a *header* library. Each header file should contain
a single function.  The function may have multiple prototypes. All functions
should use the igl namespace and should adhere to the conventions and styles
listed below. 

= Dependencies =
Eigen3  Last tested with Eigen Version 3.1.2
AntTweakBar  Last tested 1.14 (see External)
OpenGL
GLUT
GLEW  Windows only

  = Optional =
  OpenMP  
  libpng  libiglpng extra only
  Mosek  libiglmosek extra only
  Matlab  libiglmatlab extra only

  = Optional (included in external/ =
  TetGen  libigltetgen extra only

= Header only =
libigl is designed to work "out-of-the-box" as a headers only library. To
include libigl in your project. You need only include the libigl/include/
directory in your include path and define the IGL_HEADER_ONLY macro. To 
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

g++ -DIGL_HEADER_ONLY -I/usr/local/igl/igl_lib/include \
  -I/opt/local/include/eigen3 example.cpp -o example

Then run this example with:

./example examples/shared/TinyTorus.obj

= Compile =
libigl may also be compiled to a static library. This is advantageous when
building a project with libigl, since the header only directive can slow down
compile times.

To build the entire libigl library producing lib/libigl.a, issue:

make lib

You may need to edit Makefile.conf accordingly. Best to give yourself an
IGL_USERNAME and add a custom install suite for yourself. Then you can enable
appropriate extras.

  = Extras =
  Once you've set up an IGL_USERNAME and enabled extras within Makefile.conf.
  You can build the extra libraries (into lib/ligiglpng.a, lib/libiglmatlab.a,
  lib/libigltetgen.a and lib/libiglmosek.a) by issuing:

  make extras

  = Examples =
  You can make a slew of examples by issuing:

  make examples

  = External =
  Finally there are a number of external libraries that we include in
  ./external/ because they are either difficult to obtain or they have been
  patched for easier use with libigl. Please see the respective readmes in
  those directories.

= Examples =
To get started, we advise that you take a look at a few examples:

  ./examples/hello-world/

  ./examples/meshio/

  ./examples/basic-topology/

  ./examples/ReAntTweakBar/

= Development =
Further documentation is listed for developers in tutorial.html,
style_guidelines.html

= License =
For now, all files are Copyright 2013 - Alec Jacobson, Daniele Panozzo, Olga
Diamanti, Kenshi Takayama, Leo Sacht, Interactive Geometry Lab - ETH Zurich
unless otherwise noted. Soon we hope to upgrade to a more liberal license.

= Zipping =
Zip this directory without .hg litter using:

make clean
zip -9 -r --exclude=*.hg*  libigl.zip ../libigl


= Contact =
libigl is a group endeavor led by Alec Jacobson and Daniele Panozzo. Please
contact alecjacobson@gmail.com if you have questions or comments. We are happy
to get feedback! Enjoy!
