// ======================================================================== //
// Copyright 2009-2013 Intel Corporation                                    //
//                                                                          //
// Licensed under the Apache License, Version 2.0 (the "License");          //
// you may not use this file except in compliance with the License.         //
// You may obtain a copy of the License at                                  //
//                                                                          //
//     http://www.apache.org/licenses/LICENSE-2.0                           //
//                                                                          //
// Unless required by applicable law or agreed to in writing, software      //
// distributed under the License is distributed on an "AS IS" BASIS,        //
// WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied. //
// See the License for the specific language governing permissions and      //
// limitations under the License.                                           //
// ======================================================================== //

Embree is a collection of high-performance ray tracing kernels,
developed at Intel. The target user of Embree are graphics application
engineers that want to improve the performance of their application by
leveraging the optimized ray tracing kernels of Embree. The kernels
are optimized for photo-realistic rendering on the latest Intel(R)
processors with support for SSE, AVX, AVX2, and the 16-wide Xeon
Phi(TM) vector instructions. Embree supports runtime code selection to
choose the traversal and build algorithms that best matches the
instruction set of your CPU. We recommend using Embree through its API
to get the highest benefit from future improvements. Embree is
released as Open Source under the Apache 2.0 license.

Embree supports applications written with the Intel SPMD Programm
Compiler (ISPC, http://ispc.github.com) by also providing an ISPC
interface to the core ray tracing algorithms. This makes it possible
to write a renderer in ISPC that leverages SSE, AVX, AVX2, and Xeon
Phi(TM) instructions without any code change. ISPC also supports
runtime code selection, thus ISPC will select the best code path for
your application, while Embree selects the optimal code path for the
ray tracing algorithms.

Embree contains algorithms optimized for incoherent workloads (e.g.
Monte Carlo ray tracing algorithms) and coherent workloads
(e.g. primary visibility and hard shadow rays). For standard CPUs, the
single-ray traversal kernels in Embree provide the best performance
for incoherent workloads and are very easy to integrate into existing
rendering applications. For Xeon Phi(TM), a renderer written in ISPC
using the default hybrid ray/packet traversal algorithms have shown to
perform best, but requires writing the renderer in ISPC. In general
for coherent workloads, ISPC outperforms the single ray mode on each
platform. Embree also supports dynamic scenes by implementing high
performance two-level spatial index structure construction algorithms.
        
In addition to the ray tracing kernels, Embree provides some tutorials
to demonstrate how to use the Embree API. Documentation on the Embree
API can be found in embree/doc/embree_api.pdf. The example
photorealistic renderer that was originally included in the Embree
kernel package is now available in a separate GIT repository.

--- Supported Platforms ---

Embree supports Windows, Linux and MacOS, each in 32bit and 64bit
modes. The code compiles with the Intel Compiler, the Microsoft
Compiler, GCC and CLANG. Using the Intel Compiler improves performance
by approximately 10%. Performance also varies across different
operating systems. Embree is optimized for Intel CPUs supporting SSE,
AVX, and AVX2 instructions, and requires at least a CPU with support
for SSE2.

The Xeon Phi(TM) version of Embree only works under Linux in 64bit
mode. For compilation of the the Xeon Phi(TM) code the Intel Compiler
is required. The host side code compiles with GCC, CLANG, and the
Intel Compiler.

--- Folder Structure ---

Once you downloaded or checked out Embree you will see the following
folder structure:
      
embree                  Embree root folder
embree/include          User API to the ray tracing kernels
embree/kernels          Embree ray tracing kernels implementation
embree/kernels/xeon     Embree kernels for Intel(R) Xeon(R) CPUs
embree/kernels/xeonphi  Embree kernels for Intel(R) Xeon Phi(TM) Accelerators
embree/tutorials        Embree tutorials

--- Compiling Embree on Linux and MacOS ---

Embree requires the Intel SPMD Compiler (ISPC) to compile. We have
tested ISPC version 1.6.0, but more recent versions of ISPC should
also work. You can download and install the ISPC binaries from
http://ispc.github.com/downloads.html. After installation, put the
path to the ispc executable permanently into your PATH.
      
  export PATH=path-to-ispc:$PATH

You additionally have to install CMake and the developer version of
GLUT. Under MaxOS, these dependencies can be installed using MacPorts:

  sudo port install cmake freeglut

Under Linux you can install these dependencies using yum. Depending
on your Linux distribution, some of these packages might already be
installed or might have slightly different names.

  sudo yum install cmake.x86_64
  sudo yum install freeglut.x86_64 freeglut-devel.x86_64
  sudo yum install libXmu.x86_64 libXi.x86_64 
  sudo yum install libXmu-devel.x86_64 libXi-devel.x86_64
        
Finally you can compile Embree using CMake. Create a build directory
and execute "ccmake .." inside this directory.
        
  mkdir build
  cd build
  cmake ..

This will open a configuration dialog where you should set the
CMAKE_BUILD_TYPE to "Release" and the compiler to "GCC", "CLANG" or
"ICC". You should also select all targets that you want Embree to
generate optimized code for. We recommend to enable TARGET_SSE41,
TARGET_AVX, and TARGET_AVX2 if you want to use Embree on standard
CPUs, and you have to enable TARGET_XEON_PHI if you want to use Embree
on Xeon Phi(TM). You need at least Intel Compiler 11.1 or GCC 4.4 to
enable AVX and Intel Compiler 12.1 or GCC 4.7 to enable AVX2. Now
press c (for configure) and g (for generate) to generate a Makefile
and leave the configuration. The code can be compiled by executing
make.

  make

The executables will be generated inside the build folder. We
recommend to finally install the Embree library and header files on
your system:

  sudo make install

--- Compiling Embree on Windows ---

Embree requires the Intel SPMD Compiler (ISPC) to compile. We have
tested ISPC version 1.6.0, but more recent versions of ISPC should
also work. You can download and install the ISPC binaries from
http://ispc.github.com/downloads.html. After installation, put the
path to ispc.exe permanently into your PATH environment variable. You
have to restart Visual Studio for this change to take effect.
      
For compilation of Embree under Windows use the Visual Studio 2008
solution file embree_vs2008.sln or Visual Studio 2010 solution file
embree_vs2010.sln. The project compiles in 32 bit and 64 bit mode. The
solution is by default setup to use the Microsoft Compiler. You can
switch to the Intel Compiler by right clicking onto the solution in
the Solution Explorer and then selecting the Intel Compiler. We
recommend using 64 bit mode and the Intel Compiler for best
performance.
      
In Visual Studio, you will find 4 build configurations, Debug (for
SSE2 debug mode), Release (for SSE2 release mode), ReleaseAVX (for AVX
release mode), and ReleaseAVX2 (for AVX2 release mode). When using the
Microsoft Compiler you can only use the Debug and Release
configuration. For enabling the ReleaseAVX configuration you need at
least Intel Compiler 11.1 and for the ReleaseAVX2 configuration you
need at least Intel Compiler 12.1.

There is a known issue with compiling the ISPC files in Visual Studio
2010, resulting in link errors for the first build. Build the project
a second time (no rebuild) for it to link properly.

We recommend enabling syntax highlighting for the .ispc source and
.isph header files. To do so open Visual Studio 2008, go to Tools ->
Options -> Text Editor -> File Extension and add the isph and ispc
extension for the "Microsoft Visual C++" editor.

--- Running the Tutorials ---

Some tutorials come as C++ and ISPC version, e.g.:

  ./tutorial00
  ./tutorial00_ispc

You can select an initial camera using the -vp (camera position), -vi
(camera lookat point), -vu (camera up vector), and -fov (field of
view) command line parameters:

  ./tutorial00 -vp 10 10 10 -vi 0 0 0

You can select the initial windows size using the -size command line
parameter, or start the tutorials in fullscreen using the -fullscreen
parameter:

  ./tutorial00 -size 1024 1024
  ./tutorial00 -fullscreen

Implementation specific parameters can be passed to the ray tracing
core through the -rtcore command line parameter, e.g.:

  ./tutorial00 -rtcore verbose=2,threads=1,accel=bvh4.triangle1

The navigation in the interactive display mode follows the camera
orbit model, where the camera revolves around the current center of
interest. With the left mouse button you can rotate around the center
of interest (the point initially set with -vi). Holding Control
pressed while klicking the left mouse button rotates the camera around
its location. You can also use the arrow keys for navigation.

--- Contact ---

Please contact embree_support@intel.com if you have questions related to
Embree or if you want to report a bug.
