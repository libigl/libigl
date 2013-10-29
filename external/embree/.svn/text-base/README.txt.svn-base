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
developed at Intel. The target user of Embree are graphics
application engineers that want to improve the performance of their 
application by leveraging the optimized ray tracing kernels of Embree.
The kernels are optimized for photo-realistic rendering on the latest
Intel® processors with support for SSE, AVX, and 16 wide Xeon Phi
vector instructions. Embree supports applications written with the
Intel SPMD Programm Compiler (ISPC, http://ispc.github.com) by
providing an ISPC interface to the core ray tracing algorithms. This 
makes it possible to write a renderer in ISPC that leverages SSE, AVX, 
and Xeon Phi instructions without any code change. 

Embree contains algorithms optimized for incoherent workloads (e.g. 
Monte Carlo ray tracing algorithms) and coherent workloads (e.g. primary
visibility and hard shadow rays). For standard CPUs, the single-ray 
traversal kernels in Embree provide the best performance for
incoherent workloads and are very easy to integrate into existing
rendering applications. For Xeon Phi, a renderer written in ISPC using 
the optimized hybrid ray/packet traversal algorithms have shown to
perform best. In general for coherent workloads, ISPC outperforms the
single ray mode on each platform.

In addition to the ray tracing kernels, Embree provides some tutorials 
and an example photo-realistic rendering engine to demonstrate how the ray
tracing kernels are used in practice and to measure the performance of
the kernels in a realistic application scenario. The Embree example
renderer is not a full featured renderer and not designed to be used
for production renderering.

Embree is released as Open Source under the Apache 2.0 license.

--- Supported Platforms ---

Embree for SSE and AVX runs on Windows, Linux and MacOSX, each in 32bit and 64bit
modes. The code compiles with the Intel Compiler, the Microsoft
Compiler and with GCC. Using the Intel Compiler improves
performance by approximately 10%. Performance also varies across different 
operating systems. Embree is optimized for Intel CPUs supporting SSSE3, 
SSE4.1 and AVX instructions.

The Xeon Phi version of Embree only works under Linux in 64bit
mode. For compilation of the the Xeon Phi code the Intel Compiler is 
required. The host side code compiles with GCC and the Intel
Compiler.

--- Compiling Embree on Linux and MacOSX ---

For compilation under Linux and MacOSX you have to install CMake (for
compilation) the developer version of GLUT (for display) and we
recommend installing the ImageMagick and OpenEXR developer packages
(for reading and writing images). 

Under MacOSX you can install these dependencies using MacPorts:

   sudo port install cmake freeglut openexr ImageMagick

Under Linux you can install the dependencies using yum:

   sudo yum install cmake.x86_64
   sudo yum install freeglut.x86_64 freeglut-devel.x86_64
   sudo yum install libXmu.x86_64 libXi.x86_64 libXmu-devel.x86_64 libXi-devel.x86_64
   sudo yum install OpenEXR.x86_64 OpenEXR-devel.x86_64
   sudo yum install ImageMagick.x86_64 ImageMagick-c++.x86_64 ImageMagick-devel.x86_64 ImageMagick-c++-devel.x86_64 

When enabling the ISPC renderer or tutorials you also have to compile and install ISPC first (see below). To compile 
the code using CMake create a build directory and execute ccmake .. 
inside this directory. 

   mkdir build
   cd build
   ccmake ..

This will open a configuration dialog where you should set the build
mode to “Release”, the compiler target to either SSSE3, SSE4.1, SSE4.2, or
AVX, and possibly enable the ICC compiler for better performance. You
can also configure which parts of Embree to build:

   BUILD_SINGLE_RAY_DEVICE     : Single ray device for CPU operating on individual rays
   BUILD_SINGLE_RAY_DEVICE_KNC : Single ray device for Xeon Phi operating on individual rays
   BUILD_ISPC_DEVICE_SSE       : ISPC CPU device using SSE (ray packets of size 4)
   BUILD_ISPC_DEVICE_AVX       : ISPC CPU device using AVX (ray packets of size 8)
   BUILD_ISPC_DEVICE_KNC       : ISPC Xeon Phi device (ray packets of size 16)
   BUILD_TUTORIALS_SSE         : Compile ISPC tutorials for SSE
   BUILD_TUTORIALS_AVX         : Compile ISPC tutorials for AVX
   BUILD_TUTORIALS_KNC         : Compile ISPC tutorials for KNC

For building the AVX tutorials and AVX rendering device you have to
enable AVX as compilation target.

Press c (for configure) and g (for generate) to generate a Makefile
and leave the configuration. The code can now be compiled by executing
make. The executables will be generated in the build folder.

      make

--- Compiling Embree on Windows ---

For compilation under Windows use the Visual Studio
2008 solution files. Use embree_ispc.sln to compile 
Embree with ISPC support and embree.sln to compile
Embree without ISPC support. Inside Visual Studio you 
can switch between the Microsoft Compiler and the Intel 
Compiler by right clicking on the solution and then 
selecting the compiler. The project compiles with both 
compilers in 32 bit and 64 bit mode. We recommend using 
64 bit mode and the Intel Compiler for best performance. 

The solution file requires ISPC to be installed properly (see
below). For compiling the solution without ISPC, simply delete all tutorials,
device_ispc, and embree_ispc projects from the solution.

When using the Microsoft Compiler, SSE4 is enabled by default in the 
codebase. Disabling this default setting by removing the __SSE4_1__
and __SSE4_2__ define in common/sys/platform.h is necessary when SSE4 
is not supported on your system.

We recommend enabling syntax highlighting for the .ispc source 
and .isph header files. To do so open Visual Studio 2008, go to 
Tools -> Options -> Text Editor -> File Extension and add the isph
and ispc extension for the "Microsoft Visual C++" editor.

--- Installing ISPC ---

For the ISPC projects of Embree to work you have to install
ISPC from ispc.github.com. Best use ISPC v1.4.2 as we used that 
version for testing. You can download precompiled ISPC binaries or 
compile ISPC from sources. We recommend using the precompiled
binaries. After installing ISPC you have to set the ISPC_DIR
environment variable and put the ispc executable into your path:

  export ISPC_DIR=path-to-ispc
  export PATH=path-to-ispc:$PATH

Best set the ISPC_DIR variable and PATH permanently.

--- Folder structure ---

    embree                               embree ray tracing kernels
    embree/include/                      user API to the ray tracing kernels
    embree_ispc                          ISPC binding of the ray tracing kernels
    embree_ispc/include                  ISPC user API to the ray tracing kernels
    models                               Simple models for testing
    examples                             Example applications for Embree
    examples/renderer                    Photo realistic renderer building on Embree
    examples/renderer/device_singleray   Single ray implementation of renderer
    examples/renderer/device_ispc        ISPC implementation of renderer
    examples/renderer/viewer             Viewer frontend for the renderer.
    examples/tutorialXX                  Simple ISPC tutorials for embree.

--- Running the Embree example renderer ---

This section describes how to run the embree example renderer. Execute 
embree -help for a complete list of parameters. Embree ships with a few simple test
scenes, each consisting of a scene file (.xml or .obj) and an Embree
command script file (.ecs). The command script file contains command
line parameters that set the camera parameters, lights and render
settings. The following command line will render the Cornell box
scene with 16 samples per pixel and write the resulting image to the
file cb.tga in the current directory:

   renderer -c ../../models/cornell_box.ecs -spp 16 -o cb.tga

To interactively display the same scene, enter the following command:

   renderer -c ../../models/cornell_box.ecs

A window will open and you can control the camera using the mouse and
keyboard. Pressing c in interactive mode outputs the current camera
parameters, pressing r enables or disables the progressive refinement
mode.

By default the renderer uses the single ray device. For selecting a
different device use the -device command line parameter at the very 
beginning:

  renderer -device singleray -c ../../models/cornell_box.ecs
  renderer -device ispc -c ../../models/cornell_box.ecs

Under Linux and MacOSX the ISPC device also carries the instruction
set it was compiled for:

  renderer -device ispc_sse -c ../../models/cornell_box.ecs
  renderer -device ispc_avx -c ../../models/cornell_box.ecs
  renderer -device ispc_knc -c ../../models/cornell_box.ecs

For each of these configurations embree will use the best traversal 
algorithm automatically.

--- Using Embree renderer in network mode ---

For using the network device start the embree server on some machine:

  renderer_server

Make sure that port 8484 is not blocked by the firewall. Now you can 
connect from a second machine to the network server:

  renderer -connect ip_of_network_server -c ../../models/cornell_box.ecs
  
--- Navigation ---

The navigation in the interactive display mode follows the camera
orbit model, where the camera revolves around the current center of
interest. The camera navigation assumes the y-axis to point
upwards. If your scene is modelled using the z-axis as up axis we
recommend rotating the scene.

	LMB: Rotate around center of interest
	MMB: Pan
	RMB: Dolly (move camera closer or away from center of interest)
	Strg+LMB: Pick center of interest
	Strg+Shift+LMB: Pick focal distance
	Alt+LMB: Roll camera around view direction
	L: Decrease lens radius by one world space unit
	Shift+L: Increase lens radius by one world space unit

--- Contact ---

Please contact embree_support@intel.com if you have questions related to
Embree or if you want to report a bug.
