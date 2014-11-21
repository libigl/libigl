# Compiling libigl as a static library

> Warning: compiling libigl as a static library is considerably more difficult
> than using it as a header-only library (see `../README.md` instead). Do it
> only if you are experienced with C++ and you want to improve your compilation
> times.

Libigl is developed most often on Mac OS X, though has current users in Linux
and Windows.

### Linux/Mac OS X/Cygwin ###

Libigl may also be compiled to a static library. This is advantageous when
building a project with libigl, since when used as an header-only library can
slow down compile times.

To build the entire libigl library producing `lib/libigl.a`, issue:

    cd build
    make lib

You may need to edit `Makefile.conf` accordingly. Best to give yourself an
`IGL_USERNAME` and add a custom install suite for yourself. Then you can enable
appropriate extras.

#### Extras ####
Once you've set up an `IGL_USERNAME` and enabled extras within Makefile.conf.
You can build the extra libraries (into `lib/ligiglpng.a`, `lib/libiglmatlab.a`,
`lib/libigltetgen.a`, `lib/libiglmosek.a`, etc.) by issuing:

    cd build
    make extras

#### Examples ####
You can make a slew of examples by issuing:

    cd build
    make examples

#### External ####
Finally there are a number of external libraries that we include in
`./external/` because they are either difficult to obtain or they have been
patched for easier use with libigl. Please see the respective readmes in those
directories.


##### Installing AntTweakBar #####
To build the a static AntTweakBar library on Mac OS X issue:

    cd external/AntTweakBar/src
    make -f Makefile.osx.igl

##### Installing Tetgen #####
To build the tetgen library and executable on Mac OS X issue:

    cd external/tetgen
    make clean
    rm -f obj/*.o
    make -f Makefile.igl tetgen
    rm -f obj/*.o
    make -f Makefile.igl tetlib

##### Installing medit #####
To build the igl version of the medit executable on Mac OS X issue:

    cd external/medit
    make -C libmesh
    make -f Makefile.igl medit

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

You may need to install libpng. Systems with X11 might find this already
installed at `/usr/X11/lib`.


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

## Extras ##
Libigl compartmentalizes dependences via its organization into a _main_ libigl
library and "extras."


### bbw ###
This library extra contains functions for computing Bounded Biharmonic Weights, can
be used with and without the [mosek](#mosek) extra via the `IGL_NO_MOSEK`
macro.

### boost ###
This library extra utilizes the graph functions in the boost library for find
connected components and performing breadth-first traversals.

### cgal ###
This library extra utilizes CGAL's efficient and exact intersection and
proximity queries.

### embree ###
This library extra utilizes embree's efficient ray tracing queries.

### matlab ###
This library extra provides support for reading and writing `.mat` workspace
files, interfacing with Matlab at run time and compiling mex functions.

### mosek ###
This library extra utilizes mosek's efficient interior-point solver for
quadratic programs.

### png ###
This library extra uses `libpng` and `YImage` to read and write `.png` files.

### svd3x3 ###
This library extra implements "as-rigid-as-possible" (ARAP) deformation
techniques using the fast singular value decomposition routines
written specifically for 3x3 matrices to use `SSE` intrinsics. This extra can
still be compiled without sse support and support should be determined
automatically at compile time via the `__SSE__` macro.

### tetgen ###
This library extra provides a simplified wrapper to the tetgen 3d tetrahedral meshing
library.

### viewer ###
This library extra utilizes glfw and glew to open an opengl context and launch
a simple mesh viewer.

### xml ###
This library extra utilizes tinyxml2 to read and write serialized classes
containing Eigen matrices and other standard simple data-structures.

## Development ##
Further documentation for developers is listed in 
[style_guidelines.html](../style_guidelines.html).

## License ##
See `LICENSE.txt`

## Zipping ##
Zip this directory without .git litter and binaries using:

    git archive -prefix=libigl/ -o libigl.zip master

## Explicit specialization of templated functions

Special care must be taken by the developers of each function and
class in the libigl library that uses C++ templates. If this function
is intended to be compiled into the statically linked libigl library
then function is only compiled for each <i>explicitly</i> specialized
declaration. These should be added at the bottom of the corresponding
.cpp file surrounded by a

    #ifdef IGL_STATIC_LIBRARY

Of course, a developer may not know ahead of time which
specializations should be explicitly included in the igl static lib.
One way to find out is to add one explicit specialization for each
call in one's own project. This only ever needs to be done once for
each template.

The process is somewhat mechanical using a linker with reasonable error
output.

Supposed for example we have compiled the igl static lib, including the
cat.h and cat.cpp functions, without any explicit instanciation. Say
using the makefile in the `libigl` directory:

    cd $LIBIGL
    make

Now if we try to compile a project and link against it we may get
an error like:


    Undefined symbols for architecture x86_64:
    "Eigen::Matrix<int, -1, -1, 0, -1, -1> igl::cat<Eigen::Matrix<int, -1, -1, 0, -1, -1> >(int, Eigen::Matrix<int, -1, -1, 0, -1, -1> const&, Eigen::Matrix<int, -1, -1, 0, -1, -1> const&)", referenced from:
    uniform_sample(Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, Eigen::Matrix<int, -1, -1, 0, -1, -1> const&, int, double, Eigen::Matrix<double, -1, -1, 0, -1, -1>&)in Skinning.o
    "Eigen::SparseMatrix<double, 0, int> igl::cat<Eigen::SparseMatrix<double, 0, int> >(int, Eigen::SparseMatrix<double, 0, int> const&, Eigen::SparseMatrix<double, 0, int> const&)", referenced from:
    covariance_scatter_matrix(Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, Eigen::Matrix<int, -1, -1, 0, -1, -1> const&, ArapEnergy, Eigen::SparseMatrix<double, 0, int>&)in arap_dof.o
    arap_rhs(Eigen::Matrix<double, -1, -1, 0, -1, -1> const&, Eigen::Matrix<int, -1, -1, 0, -1, -1> const&, ArapEnergy, Eigen::SparseMatrix<double, 0, int>&)in arap_dof.o

This looks like a mess, but luckily we don't really need to read it
all. Just copy the first part in quotes

    Eigen::Matrix<int, -1, -1, 0, -1, -1> igl::cat<Eigen::Matrix<int, -1, -1, 0, -1, -1> >(int, Eigen::Matrix<int, -1, -1, 0, -1, -1> const&, Eigen::Matrix<int, -1, -1, 0, -1, -1> const&)

, then append it
to the list of explicit template specializations at the end of
`cat.cpp` after the word
**template** and followed by a semi-colon.
Like this:

    #ifdef IGL_STATIC_LIBRARY
    // Explicit template specialization
    template Eigen::Matrix<int, -1, -1, 0, -1, -1> igl::cat<Eigen::Matrix<int, -1, -1, 0, -1, -1> >(int, Eigen::Matrix<int, -1, -1, 0, -1, -1> const&, Eigen::Matrix<int, -1, -1, 0, -1, -1> const&);
    #endif

Then you must recompile the IGL static library.

    cd $LIBIGL
    make

And try to compile your project again, potentially repeating this
process until no more symbols are undefined.

`It may be useful to check that you code compiles with
no errors first using the headers-only version to be sure that all errors are from missing template
specializations.`

If you're using make then the following command will
reveal each missing symbol on its own line:

    make 2>&1 | grep "referenced from" | sed -e "s/, referenced from.*//"

Alternatively you can use the `autoexplicit.sh` function
which (for well organized .h/.cpp pairs in libigl) automatically
create explicit instanciations from your compiler's error messages.
Repeat this process until convergence:

    cd /to/your/project
    make 2>$LIBIGL/make.err
    cd $LIBIGL
    cat make.err | ./autoexplicit.sh
    make clean
    make


### Benefits of static library

* **Faster compile time**: Because the libigl library
    is already compiled, only the new code in ones project must be
    compiled and then linked to IGL. This means compile times are
    generally faster.
* **Debug or optimized**: The IGL static
    library may be compiled in debug mode or optimized release mode
    regardless of whether one's project is being optimized or
    debugged.

### Drawbacks of static library

*  **Hard to use templates**: Special
    care</a> (by the developers of the library) needs to be taken when
    exposing templated functions.

# Compressed .h/.cpp pair
Calling the script:

    scripts/compress.sh igl.h igl.cpp

will create a single header `igl.h` and a single cpp file `igl.cpp`.

Alternatively, you can also compress everything into a single header file:

    scripts/compress.sh igl.h

### Benefits of compressed .h/.cpp pair

* **Easy incorporation**: This can be easily incorporated
  into external projects.

### Drawbacks of compressed .h/.cpp pair

* **Hard to debug/edit**: The compressed files are
  automatically generated. They're huge and should not be edited. Thus
  debugging and editting are near impossible.

* **Compounded dependencies**:
  An immediate disadvantage of this
  seems to be that even to use a single function (e.g.
  `cotmatrix`), compiling and linking against
  `igl.cpp` will require linking to all of `libigl`'s
  dependencies (`OpenGL`, `GLUT`,
  `AntTweakBar`, `BLAS`). However, because all
  depencies other than Eigen should be encapsulated between
  `#ifndef` guards (e.g. `#ifndef IGL_NO_OPENGL`, it
  is possible to ignore certain functions that have such dependencies.</li>
  <li><strong>Long compile:</strong> Compiling `igl.cpp` takes a long time and isn't easily parallelized (no `make -j12` equivalent).</li>
  </ul>

Here's a tiny test example using `igl.h` and `igl.cpp`. Save the following in `test.cpp`:

    #include <igl.h>
    #include <Eigen/Core>

    int main(int argc, char * argv[])
    {
    Eigen::MatrixXd V;
    Eigen::MatrixXi F;
    return (argc>=2 &amp;&amp; igl::read_triangle_mesh(argv[1],V,F)?0:1);
    }

Then compile `igl.cpp` with:

    g++ -o igl.o -c igl.cpp -I/opt/local/include/eigen3 -DIGL_NO_OPENGL -DIGL_NO_ANTTWEAKBAR

Notice that we're using `-DIGL_NO_OPENGL -DIGL_NO_ANTTWEAKBAR` to disable any libigl dependencies on OpenGL and AntTweakBar.

Now compile `test.cpp` with:

    g++ -g -I/opt/local/include/eigen3/ -I/usr/local/igl/libigl/ -L/usr/local/igl/libigl/ -ligl -DIGL_NO_OPENGL -DIGL_NO_ANTTWEAKBAR -o test

Try running it with:

    ./test path/to/mesh.obj


The following bash one-liner will find all source files that contain the string `OpenGL` but don't contain and `IGL_NO_OPENGL` guard:

    grep OpenGL `grep -L IGL_NO_OPENGL include/igl/*`

### Optional ###
- OpenGL (disable with `IGL_NO_OPENGL`)
    * OpenGL >= 4 (enable with `IGL_OPENGL_4`)
- AntTweakBar  (disable with `IGL_NO_ANTTWEAKBAR`) Last tested 1.16 (see
  `libigl/external/AntTweakBar`)
- GLEW  Windows and Linux
- OpenMP
- libpng  libiglpng extra only
- Mosek  libiglmosek extra only
- Matlab  libiglmatlab extra only
- boost  libiglboost, libiglcgal extra only
- SSE/AVX  libiglsvd3x3 extra only
- CGAL  libiglcgal extra only
    * boost
    * gmp
    * mpfr
- CoMiSo libcomiso extra only

### Optional (included in external/) ###
- TetGen  libigltetgen extra only
- Embree  libiglembree extra only
- tinyxml2  libiglxml extra only
- glfw libviewer extra only
- LIM  liblim extra only
