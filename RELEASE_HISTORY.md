title: libigl Tutorial
author: Alec Jacobson
date: 17 June 2015
css: tutorial/style.css
html header:   <script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
<link rel="stylesheet" href="http://yandex.st/highlightjs/7.3/styles/default.min.css">
<script src="http://yandex.st/highlightjs/7.3/highlight.min.js"></script>
<script>hljs.initHighlightingOnLoad();</script>

# Libigl version tracking

Version | Short description
--------|----------------------------------------------------------------------
1.2.1   | Reorganization opengl-dependent functions: opengl and opengl2 extras
1.2.0   | Reorganization of "extras", rm deprecated funcs, absorb boost & svd3x3
1.1.7   | Switch build for static library to cmake.
1.1.6   | Major boolean robustness fix, drop CGAL dependency for AABB/distances
1.1.5   | Bug fix in booleans
1.1.4   | Edge collapsing and linear program solving
1.1.3   | Bug fixes in active set and boundary_conditions
1.1.1   | PLY file format support
1.1.0   | Mesh boolean operations using CGAL and cork, implementing [Attene 14]
1.0.3   | Bone heat method
1.0.2   | Bug fix in winding number code
1.0.1   | Bug fixes and more CGAL support
1.0.0   | Major beta release: many renames, tutorial, triangle, org. build
0.4.6   | Generalized Winding Numbers
0.4.5   | CGAL extra: mesh selfintersection
0.4.4   | STL file format support
0.4.3   | ARAP implementation
0.4.1   | Migrated much of the FAST code including extra for Sifakis' 3x3 svd
0.4.0   | Release under MPL2 license
0.3.7   | Embree2.0 support
0.3.6   | boost extra, patches, mosek 7 support, libiglbbw (mosek optional)
0.3.5   | More examples, naive primitive sorting
0.3.3   | Many more examples, ambient occlusion with Embree.
0.3.1   | Linearly dependent constraints in min_quad_with_fixed, SparseQR buggy
0.3.0   | Better active set method support
0.2.3   | More explicits, active set method, opengl/anttweakbar guards
0.2.2   | More explicit instanciations, faster sorts and uniques
0.2.1   | Bug fixes in barycenter and doublearea found by Martin Bisson
0.2.0   | XML serializer more stable and fixed bug in remove_duplicate_vertices
0.1.8   | Embree and xml (windows only) extras
0.1.5   | Compilation on windows, bug fix for compilation with cygwin
0.1.1   | Alpha release with core functions, extras, examples

## Version 1.2 Changes ##
This change introduces better organization of dependencies and removes some
deprecated/repeated functions. The 3x3 svd code and dependent functions
(including ARAP) were absorbed into the main library. Similarly, the boost
dependency extra was absorbed.


### External libraries as git subrepos ###
The core functionality of libigl (still) just depends on stl, c++11 and Eigen.
There are additional _optional_ dependencies (e.g. CGAL, embree, glfw, tetgen,
triangle). Libigl functions using these are located (still) in sub-folders of
the include directory (e.g.  `include/igl/cgal/`, `include/igl/embree/`). Prior
to version 1.2 we included copies of the code for some of these dependencies in the
`external/` directory. As of
version 1.2, these have been replaced with git sub-repos. If you have cloned
libigl _before version 1.2_ then you should issue 

    git submodule update --init --recursive

### Deprecated/repeated functions ###

Old                                     | New
--------------------------------------- | -----------------------------------
`igl::angles`                           | `igl::internal_angles`
`igl::get_modifiers`                    | [deleted]
`igl::nchoosek(offset,K,N,std::vector)` | `igl::nchoosek(Eigen,K,Eigen)`
`#include <igl/boost/components.h>`     | `#include <igl/components.h>`
`#include <igl/boost/bfs_orient.h>`     | `#include <igl/bfs_orient.h>`
`#include <igl/boost/orientable_patches.h>` | `#include <igl/orientable_patches.h>`
`#include <igl/svd3x3/arap.h>`          | `#include <igl/arap.h>`
`#include <igl/svd3x3/arap_dof.h>`      | `#include <igl/arap_dof.h>`
`#include <igl/svd3x3/fit_rotations.h>` | `#include <igl/fit_rotations.h>`
`#include <igl/svd3x3/polar_svd3x3.h>`  | `#include <igl/polar_svd3x3.h>`
`#include <igl/svd3x3/svd3x3.h>`        | `#include <igl/svd3x3.h>`
`#include <igl/svd3x3/svd3x3_avx.h>`    | `#include <igl/svd3x3_avx.h>`
`#include <igl/svd3x3/svd3x3_sse.h>`    | `#include <igl/svd3x3_sse.h>`


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
 - We have organized the static library build
 - The previous `igl::grad` function, which computed the per-triangle gradient
   of a per-vertex scalar function has been replaced. Now `igl::grad` computes
   the linear operator (previous computed using `igl::gradMat`). The gradient
   values can still be recovered by multiplying the operator against the scalar
   field as a vector and reshaping to have gradients per row.
 - `MASSMATRIX_*` has become `MASSMATRIX_TYPE_*`
 - The function `igl::project_normals`, which cast a line for each vertex of
   mesh _A_ in the normal direction and found the closest intersection along
   these lines with mesh _B_, has been removed.
