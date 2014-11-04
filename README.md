libigl - A simple c++ geometry processing library
=================================================

<http://igl.ethz.ch/projects/libigl/>
<https://github.com/alecjacobson/libigl/>

Copyright 2014 - Alec Jacobson, Daniele Panozzo, Olga Diamanti, Kenshi
Takayama, Leo Sacht, Wenzel Jacob, etc.

This is first and foremost a *header* library. Each header file should contain
a single function.  The function may have multiple prototypes. All functions
should use the igl namespace and should adhere to the conventions and styles
listed below.

> **New:** As of 1 July 2014, we have released libigl as beta version 1.0.
> There are a number of changes we collected for this release to minimize
> confusion and changes to how you use libigl. See [Version 1.0
> Changes][version1.0changes].

## Installation ##
Libigl is a *header* library. You do **not** need to build anything to install.
Simple add `igl/` to your include path and include relevant headers. Here's a
small "Hello, World" program:

```cpp
#include <igl/cotmatrix.h>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <iostream>
int main()
{
  Eigen::MatrixXd V(4,2);
  V<<0,0,
     1,0,
     1,1,
     0,1;
  Eigen::MatrixXi F(2,3);
  F<<0,1,2,
     0,2,3;
  Eigen::SparseMatrix<double> L;
  igl::cotmatrix(V,F,L);
  std::cout<<"Hello, mesh: "<<std::endl<<L*V<<std::endl;
  return 0;
}
```

If you save this in `hello.cpp`, then on `gcc` with Eigen installed via
macports for example you could compile this with:

```bash
gcc -I/opt/local/include/eigen3 -I./igl/ hello.cpp -o hello
```

Running `./hello` would then produce

```
Hello, mesh:
 0.5  0.5
-0.5  0.5
-0.5 -0.5
 0.5 -0.5
```

## Tutorial ##

As of version 1.0, libigl includes an introductory
[tutorial](http://libigl.github.io/libigl/tutorial/tutorial.html) that covers
its basic functionalities.

## Dependencies ##
- Eigen3  Last tested with Eigen Version 3.2

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

## Header only ##
Libigl is designed to work "out-of-the-box" as a headers only library. To
include libigl in your project. You need only include the libigl/include/
directory in your include path. To
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

    g++ -I/usr/local/igl/libigl/include \
      -I/opt/local/include/eigen3 example.cpp -o example

Then run this example with:

    ./example examples/shared/TinyTorus.obj

## Compilation as a static library ##
Libigl is developed most often on Mac OS X, though has current users in Linux
and Windows.

### Linux/Mac OS X/Cygwin ###

Libigl may also be compiled to a static library. This is advantageous when
building a project with libigl, since when used as an header-only library can
slow down compile times.

To build the entire libigl library producing lib/libigl.a, issue:

    cd build
    make lib

You may need to edit Makefile.conf accordingly. Best to give yourself an
`IGL_USERNAME` and add a custom install suite for yourself. Then you can enable
appropriate extras.

#### Extras ####
Once you've set up an `IGL_USERNAME` and enabled extras within Makefile.conf.
You can build the extra libraries (into lib/ligiglpng.a, lib/libiglmatlab.a,
lib/libigltetgen.a, lib/libiglmosek.a, etc.) by issuing:

    cd build
    make extras

#### Examples ####
You can make a slew of examples by issuing:

    cd build
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
Further documentation for developers is listed in tutorial.html,
style_guidelines.html

## License ##
See `LICENSE.txt`

## Zipping ##
Zip this directory without .git litter and binaries using:

    git archive -prefix=libigl/ -o libigl.zip master

## Version 1.0 Changes ##
Our beta release marks our confidence that this library can be used outside of
casual experimenting. To maintain order, we have made a few changes which
current users should read and adapt their code accordingly.

### Renamed functions ###
The following table lists functions which have changed name as of version
1.0.0:

Old                              | New
-------------------------------- | -------------------------------------
`igl::add_barycenter`            | `igl::false_barycentric_subdivision`
`igl::areamatrix`                | `igl::vector_area_matrix`
`igl::barycentric2global`        | `igl::barycentric_to_global`
`igl::boundary_faces`            | `igl::boundary_facets`
`igl::boundary_vertices_sorted`  | `igl::boundary_loop`
`igl::cotangent`                 | `igl::cotmatrix_entries`
`igl::edgetopology`              | `igl::edge_topology`
`igl::gradMat`                   | `igl::grad`
`igl::is_manifold`               | `igl::is_edge_manifold`
`igl::mexStream`                 | `igl::MexStream`
`igl::moveFV`                    | `igl::average_onto_vertices`
`igl::moveVF`                    | `igl::average_onto_faces`
`igl::plot_vector`               | `igl::print_vector`
`igl::pos`                       | `igl::HalfEdgeIterator`
`igl::plane_project`             | `igl::project_isometrically_to_plane`
`igl::project_points_mesh`       | `igl::line_mesh_intersection`
`igl::read`                      | `igl::read_triangle_mesh`
`igl::removeDuplicates.cpp`      | `igl::remove_duplicates`
`igl::removeUnreferenced`        | `igl::remove_unreferenced`
`igl::tt`                        | `igl::triangle_triangle_adjacency`
`igl::vf`                        | `igl::vertex_triangle_adjacency`
`igl::write`                     | `igl::write_triangle_mesh`
`igl::manifold_patches`          | `igl::orientable_patches`
`igl::selfintersect`             | `igl::remesh_self_intersections`
`igl::project_mesh`              | `igl::line_mesh_intersection`
`igl::triangulate`               | `igl::polygon_mesh_to_triangle_mesh`
`igl::is_manifold`               | `igl::is_edge_manifold`
`igl::triangle_wrapper`          | `igl::triangulate`

### Miscellaneous ###
 - To match interfaces provided by (all) other quadratic optimization
   libraries, `igl::min_quad_with_fixed` and `igl::active_set` now expect as
   input twice the quadratic coefficients matrix, i.e. the Hessian. For
   example, `igl::min_quad_with_fixed(H,B,...)` minimizes $\frac{1}{2}x^T H
   x+x^T B$.
 - We have inverted the `IGL_HEADER_ONLY` macro to `IGL_STATIC_LIBRARY`. To
   compile using libigl as a header-only library, simply include headers and
   libigl in the header search path. To link to libigl, you must define the
   `IGL_STATIC_LIBRARY` macro at compile time and link to the `libigl*.a`
   libraries.
 - Building libigl as a static library is now more organized. There is a
   `build/` directory with Makefiles for the main library (`Makefile`) and each
   dependency (e.g. `Makefile_mosek` for `libiglmosek.a`)
 - `igl::polar_svd` now always returns a rotation in `R`, never a reflection.
   This mirrors the behavior of `igl::polar_svd3x3`.  Consequently the `T`
   part may have negative skews.
 - We have organized the static
 - The previous `igl::grad` function, which computed the per-triangle gradient
   of a per-vertex scalar function has been replaced. Now `igl::grad` computes
   the linear operator (previous computed using `igl::gradMat`). The gradient
   values can still be recovered by multiplying the operator against the scalar
   field as a vector and reshaping to have gradients per row.
 - `MASSMATRIX_*` has become `MASSMATRIX_TYPE_*`
 - The function `igl::project_normals`, which cast a line for each vertex of
   mesh _A_ in the normal direction and found the closest intersection along
   these lines with mesh _B_, has been removed.

## Contact ##
Libigl is a group endeavor led by Alec Jacobson and Daniele Panozzo. Please
contact [alecjacobson@gmail.com](mailto:alecjacobson@gmail.com) if you have
questions or comments. We are happy to get feedback! Enjoy!

If you're using libigl in your projects, quickly [drop us a
note](mailto:alecjacobson@gmail.com). Tell us who you are and what you're using
it for. This helps us apply for funding and justify spending time maintaining
this.

If you find bugs or have problems please use our [github issue tracking
page](https://github.com/libigl/libigl/issues).

## Academic citation ##
If you use libigl in your research projects, please cite the papers we
implement as appropriate. To cite the library in general, you could use this
BibTeX entry:

```bibtex
@misc{libigl,
  title = {{libigl}: A simple {C++} geometry processing library},
  author = {Alec Jacobson and Daniele Panozzo and others},
  note = {http://igl.ethz.ch/projects/libigl/},
  year = {2014},
}
```
