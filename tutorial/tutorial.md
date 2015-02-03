title: libigl Tutorial
author: Daniele Panozzo and Alec Jacobson
date: 07 November 2015
css: style.css
html header:   <script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
<link rel="stylesheet" href="http://yandex.st/highlightjs/7.3/styles/default.min.css">
<script src="http://yandex.st/highlightjs/7.3/highlight.min.js"></script>
<script>hljs.initHighlightingOnLoad();</script>

# libigl tutorial notes

#### as presented by Daniele Panozzo and Alec Jacobson at SGP Graduate School 2015

![](images/libigl-logo.jpg)

Libigl is an open source C++ library for geometry processing research and
development.  Dropping the heavy data structures of tradition geometry
libraries, libigl is a simple header-only library of encapsulated functions.
This combines the rapid prototyping familiar to Matlab or Python programmers
with the performance and versatility of C++.  The tutorial is a self-contained,
hands-on introduction to libigl.  Via interactive, step-by-step examples, we
demonstrate how to accomplish common geometry processing tasks such as
computation of differential quantities and operators, real-time deformation,
parametrization, numerical optimization and remeshing. Each section of the
lecture notes links to a cross-platform example application.

# Table of contents

* [Chapter 1: Introduction to libigl](#100)
    * [Libigl design principles](#100b)
    * [101 Mesh representation](#101)
    * [102 Visualizing surfaces](#102)
    * [103 Interaction with keyboard and mouse](#103)
    * [104 Scalar field visualization](#104)
    * [105 Overlays](#105)
* [Chapter 2: Discrete Geometric Quantities and
  Operators](#chapter2:discretegeometricquantitiesandoperators)
    * [201 Normals](#normals)
        * [Per-face](#per-face)
        * [Per-vertex](#per-vertex)
        * [Per-corner](#per-corner)
    * [202 Gaussian Curvature](#gaussiancurvature)
    * [203 Curvature Directions](#curvaturedirections)
    * [204 Gradient](#gradient)
    * [204 Laplacian](#laplacian)
        * [Mass matrix](#massmatrix)
        * [Alternative construction of
          Laplacian](#alternativeconstructionoflaplacian)
* [Chapter 3: Matrices and Linear Algebra](#chapter3:matricesandlinearalgebra)
    * [301 Slice](#slice)
    * [302 Sort](#sort)
        * [Other Matlab-style functions](#othermatlab-stylefunctions)
    * [303 Laplace Equation](#laplaceequation)
        * [Quadratic energy minimization](#quadraticenergyminimization)
    * [304 Linear Equality Constraints](#linearequalityconstraints)
    * [305 Quadratic Programming](#quadraticprogramming)
* [Chapter 4: Shape Deformation](#chapter4:shapedeformation)
    * [401 Biharmonic Deformation](#biharmonicdeformation)
    * [402 Polyharmonic Deformation](#polyharmonicdeformation)
    * [403 Bounded Biharmonic Weights](#boundedbiharmonicweights)
    * [404 Dual Quaternion Skinning](#dualquaternionskinning)
    * [405 As-rigid-as-possible](#as-rigid-as-possible)
    * [406 Fast automatic skinning
      transformations](#fastautomaticskinningtransformations)
        * [ARAP with grouped edge-sets](#arapwithgroupededge-sets)
* [Chapter 5: Parametrization](#500)
    * [501 Harmonic parametrization](#501)
    * [502 Least-Square Conformal Maps](#502)
    * [503 As-Rigid-As-Possible](#503)
    * [504 N-Rotationally symmetric tangent fields](#504)
    * [505 Global, seamless integer-grid parametrization](#505)
    * [506 Anisotropic remeshing using frame fields](#506)
    * [507 N-PolyVector fields](#507)
    * [508 Conjugate vector fields](#508)
    * [509 Planarization](#509)
* [Chapter 6: External libraries](#600)
    * [601 State serialization](#601)
    * [602 Mixing Matlab code](#602)
        * [Saving a Matlab workspace](#savingamatlabworkspace)
        * [Dumping Eigen matrices to copy and paste into
          Matlab](#dumpingeigenmatricestocopyandpasteintomatlab)
    * [603 Calling libigl functions from Matlab](#603)
    * [604 Triangulation of closed polygons](#604)
    * [605 Tetrahedralization of closed surfaces](#605)
    * [606 Baking ambient occlusion](#606)
    * [607 Picking vertices and faces](#607)
    * [608 Locally Injective Maps](#608)
    * [609 Boolean Operations on Meshes](#609)
* [Chapter 7: Miscellaneous](#700)
    * [701 Mesh Statistics](#701)
    * [702 Generalized Winding Number](#702)
* [Chapter 8: Outlook for continuing development](#future)

# Chapter 1 [100]

We introduce libigl with a series of self-contained examples. The purpose of
each example is to showcase a feature of libigl while applying to a practical
problem in geometry processing. In this chapter, we will present the basic
concepts of libigl and introduce a simple mesh viewer that allows to
visualize a surface mesh and its attributes. All the tutorial examples are
cross-platform and can be compiled on MacOSX, Linux and Windows.

## libigl design principles [100b]

Before getting into the examples, we summarize the main design principles in
libigl:

1. **No complex data types.** We mostly use matrices and vectors. This greatly
  favors code reusability and forces the function authors to expose all the
  parameters used by the algorithm.  

2. **Minimal dependencies.** We use external libraries only when necessary and
  we wrap them in a small set of functions.

3. **Header-only.** It is straight forward to use our library since it is only
  one additional include directory in your project. (if you are worried about
  compilation speed, it is also possible to build the library as a [static
  library](../build/))

4. **Function encapsulation.** Every function (including its full
  implementation) is contained in a pair of .h/.cpp files with the same name of
  the function.


### Downloading libigl
libigl can be downloaded from our [github
repository](https://github.com/libigl/libigl) or cloned with git:

```bash
git clone https://github.com/libigl/libigl.git
```

The core libigl functionality only depends on the C++ Standard Library and
Eigen.

The examples in this tutorial depend on [glfw](http://www.glfw.org),
[glew](http://glew.sourceforge.net) and [AntTweakBar](http://anttweakbar.sourceforge.net/doc/).
The source code of each library is bundled with libigl
and they can be compiled all at once using:

```bash
sh compile_dependencies_macosx.sh (MACOSX)
sh compile_dependencies_linux.sh (LINUX)
```

For windows, precompiled binaries are provided (Visual Studio 2014 64bit).

To build all the examples in the tutorial, you can use the CMakeLists.txt in
the tutorial folder:

```bash
cd tutorial
mkdir build
cd build
cmake -DCMAKE_BUILD_TYPE=Release ../
make
```

The examples can also be built independently using the CMakeLists.txt
inside each example folder.

A few examples in Chapter 5 requires the [CoMiSo
solver](http://www.graphics.rwth-aachen.de/software/comiso) which has to be
downloaded and compiled separately.

## Mesh representation [101]

libigl uses the [Eigen](http://eigen.tuxfamily.org/) library to encode vector
and matrices. We suggest that you keep the
[dense](http://eigen.tuxfamily.org/dox/group__QuickRefPage.html) and
[sparse](http://eigen.tuxfamily.org/dox/group__SparseQuickRefPage.html) quick
reference guides at hand while you read the examples in this tutorial.

A triangular mesh is encoded as a pair of matrices:

```cpp
Eigen::MatrixXd V;
Eigen::MatrixXi F;
```

`V` is a #N by 3 matrix which stores the coordinates of the vertices. Each
row stores the coordinate of a vertex, with its x,y and z coordinates in the first,
second and third column, respectively. The matrix `F` stores the triangle
connectivity: each line of `F` denotes a triangle whose 3 vertices are
represented as indices pointing to rows of `V`.

![A simple mesh made of 2 triangles and 4 vertices.](images/VF.png)

Note that the order of the vertex indices in `F` determines the orientation of
the triangles and it should thus be consistent for the entire surface.
This simple representation has many advantages:

1. it is memory efficient and cache friendly
2. the use of indices instead of pointers greatly simplifies debugging
3. the data can be trivially copied and serialized

libigl provides input [output] functions to read [write] many common mesh formats.
The IO functions are contained in the files read\*.h and write\*.h. As a general
rule each libigl function is contained in a pair of .h/.cpp files with the same name.
By default, the .h files include the corresponding cpp files, making the library header-only.

Reading a mesh from a file requires a single libigl function call:

```cpp
igl::readOFF("../shared/cube.off", V, F);
```

The function reads the mesh cube.off and it fills the provided `V` and `F` matrices.
Similarly, a mesh can be written in an OBJ file using:

```cpp
igl::writeOBJ("cube.obj",V,F);
```

[Example 101](101_FileIO/main.cpp) contains a simple mesh
converter from OFF to OBJ format.

## Visualizing surfaces [102]

Libigl provides an glfw-based OpenGL 3.2 viewer to visualize surfaces, their
properties and additional debugging informations.

The following code ([Example 102](102_DrawMesh/main.cpp)) is a basic skeleton
for all the examples that will be used in the tutorial.
It is a standalone application that loads a mesh and uses the viewer to
render it.

```cpp
#include <igl/readOFF.h>
#include <igl/viewer/Viewer.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF("../shared/bunny.off", V, F);

  // Plot the mesh
  igl::Viewer viewer;
  viewer.data.set_mesh(V, F);
  viewer.launch();
}
```

The function `set_mesh` copies the mesh into the viewer.
`Viewer.launch()`  creates a window, an OpenGL context and it starts the draw loop.
Additional properties can be plotted on the mesh (as we will see later),
and it is possible to extend the viewer with standard OpenGL code.
Please see the documentation in
[Viewer.h](../include/igl/Viewer/Viewer.h) for more details.

![([Example 102](102_DrawMesh/main.cpp)) loads and draws a
mesh.](images/102_DrawMesh.png)

## Interaction with keyboard and mouse [103]

Keyboard and mouse events triggers callbacks that can be registered in the
viewer. The viewer supports the following callbacks:

```cpp
bool (*callback_pre_draw)(Viewer& viewer);
bool (*callback_post_draw)(Viewer& viewer);
bool (*callback_mouse_down)(Viewer& viewer, int button, int modifier);
bool (*callback_mouse_up)(Viewer& viewer, int button, int modifier);
bool (*callback_mouse_move)(Viewer& viewer, int mouse_x, int mouse_y);
bool (*callback_mouse_scroll)(Viewer& viewer, float delta_y);
bool (*callback_key_down)(Viewer& viewer, unsigned char key, int modifiers);
bool (*callback_key_up)(Viewer& viewer, unsigned char key, int modifiers);
```

A keyboard callback can be used to visualize multiple meshes or different
stages of an algorithm, as demonstrated in [Example 103](103_Events/main.cpp), where
the keyboard callback changes the visualized mesh depending on the key pressed:

```cpp
bool key_down(igl::Viewer& viewer, unsigned char key, int modifier)
{
  if (key == '1')
  {
    viewer.data.clear();
    viewer.data.set_mesh(V1, F1);
    viewer.core.align_camera_center(V1,F1);
  }
  else if (key == '2')
  {
    viewer.data.clear();
    viewer.data.set_mesh(V2, F2);
    viewer.core.align_camera_center(V2,F2);
  }
  return false;
}
```

The callback is registered in the viewer as follows:

```cpp
viewer.callback_key_down = &key_down;
```

Note that the mesh is cleared before using set_mesh. This has to be called
every time the number of vertices or faces of the plotted mesh changes. Every
callback returns a boolean value that tells the viewer if the event has been
handled by the plugin, or if the viewer should process it normally. This is
useful, for example, to disable the default mouse event handling if you want to
control the camera directly in your code.

The viewer can be extended using plugins, which are classes that implements all
the viewer's callbacks. See the
[Viewer_plugin](../include/igl/viewer/ViewerPlugin.h) for more details.

## Scalar field visualization [104]

Colors and normals can be associated to faces or vertices using the
set_colors function:

```cpp
viewer.data.set_colors(C);
```

`C` is a #C by 3 matrix with one RGB color per row. `C` must have as many
rows as the number of faces **or** the number of vertices of the mesh.
Depending on the size of `C`, the viewer applies the color to the faces or to
the vertices.

Colors can be used to visualize a scalar function defined on a surface.  The
scalar function is converted to colors using a color transfer function, which
maps a scalar value between 0 and 1 to a color. A simple example of a scalar
field defined on a surface is the z coordinate of each point, which can be
extract from our mesh representation by taking the last column of `V`
([Example 104](104_Colors/main.cpp)). The function `igl::jet` can be used to
convert it to colors:

```cpp
Eigen::VectorXd Z = V.col(2);
igl::jet(Z,true,C);
```

The first row extracts the third column from `V` (the z coordinate of each
vertex) and the second calls a libigl functions that converts a scalar field to colors. The second parameter of jet normalizes the scalar field to lie between 0 and 1 before applying the transfer function.

![([Example 104](104_Colors/main.cpp)) igl::jet converts a scalar field to a
color field.](images/104_Colors.png)

`igl::jet` is an example of a standard function in libigl: it takes simple
types and can be easily reused for many different tasks.  Not committing to
heavy data structures types favors simplicity, ease of use and reusability.

## Overlays [105]

In addition to plotting the surface, the viewer supports the visualization of points, lines and text labels: these overlays can be very helful while developing geometric processing algorithms to plot debug informations.

```cpp
viewer.data.add_points(P,Eigen::RowVector3d(r,g,b));
```

Draws a point of color r,g,b for each row of P. The point is placed at the coordinates specified in each row of P, which is a #P by 3 matrix.

```cpp
viewer.data.add_edges(P1,P2,Eigen::RowVector3d(r,g,b);
```

Draws a line of color r,g,b for each row of P1 and P2, which connects the 3D point in to the point in P2. Both P1 and P2 are of size #P by 3.

```cpp
viewer.data.add_label(p,str);
```

Draws a label containing the string str at the position p, which is a vector of length 3.

These functions are demonstrate in [Example 105](105_Overlays/main.cpp) where
the bounding box of a mesh is plotted using lines and points.
Using matrices to encode the mesh and its attributes allows to write short and
efficient code for many operations, avoiding to write for loops. For example,
the bounding box of a mesh can be found by taking the colwise maximum and minimum of `V`:

```cpp
Eigen::Vector3d m = V.colwise().minCoeff();
Eigen::Vector3d M = V.colwise().maxCoeff();
```

![([Example 105](105_Overlays/main.cpp)) The bounding box of a mesh is shown
using overlays.](images/105_Overlays.png)

# Chapter 2: Discrete Geometric Quantities and Operators
This chapter illustrates a few discrete quantities that libigl can compute on a
mesh and the libigl functions that construct popular discrete differential
geometry operators. It also provides an introduction to basic drawing and coloring routines of our viewer.

## Normals
Surface normals are a basic quantity necessary for rendering a surface. There
are a variety of ways to compute and store normals on a triangle mesh. [Example 201](201_Normals/main.cpp) demonstrates how to compute and visualize normals with libigl.

### Per-face
Normals are well defined on each triangle of a mesh as the vector orthogonal to
triangle's plane. These piecewise-constant normals produce piecewise-flat
renderings: the surface appears non-smooth and reveals its underlying
discretization.

### Per-vertex
Normals can be computed and stored on vertices, and interpolated in the interior of the triangles to produce smooth renderings ([Phong shading](http://en.wikipedia.org/wiki/Phong_shading)).
Most techniques for computing per-vertex normals take an average of incident face normals. The main difference between these techniques is their weighting scheme: Uniform
weighting is heavily biased by the discretization choice, whereas area-based
or angle-based weighting is more forgiving.

The typical half-edge style computation of area-based weights has this structure:

```cpp
N.setZero(V.rows(),3);
for(int i : vertices)
{
  for(face : incident_faces(i))
  {
    N.row(i) += face.area * face.normal;
  }
}
N.rowwise().normalize();
```

At first glance, it might seem inefficient to loop over incident faces---and thus constructing the per-vertex normals--- without using an half-edge data structure. However, per-vertex normals may be _throwing_ each face normal to
running sums on its corner vertices:

```cpp
N.setZero(V.rows(),3);
for(int f = 0; f < F.rows();f++)
{
  for(int c = 0; c < 3;c++)
  {
    N.row(F(f,c)) += area(f) * face_normal.row(f);
  }
}
N.rowwise().normalize();
```

### Per-corner

Storing normals per-corner is an efficient and convenient way of supporting both
smooth and sharp (e.g. creases and corners) rendering. This format is common to
OpenGL and the .obj mesh file format. Often such normals are tuned by the mesh
designer, but creases and corners can also be computed automatically. Libigl
implements a simple scheme which computes corner normals as averages of
normals of faces incident on the corresponding vertex which do not deviate by more than a specified dihedral angle (e.g. 20°).

![The `Normals` example computes per-face (left), per-vertex (middle) and
per-corner (right) normals](images/fandisk-normals.jpg)

## Gaussian curvature
Gaussian curvature on a continuous surface is defined as the product of the
principal curvatures:

 $k_G = k_1 k_2.$

As an _intrinsic_ measure, it depends on the metric and
not the surface's embedding.

Intuitively, Gaussian curvature tells how locally spherical or _elliptic_ the
surface is ( $k_G>0$ ), how locally saddle-shaped or _hyperbolic_ the surface
is ( $k_G<0$ ), or how locally cylindrical or _parabolic_ ( $k_G=0$ ) the
surface is.

In the discrete setting, one definition for a "discrete Gaussian curvature"
on a triangle mesh is via a vertex's _angular deficit_:

 $k_G(v_i) = 2π - \sum\limits_{j\in N(i)}θ_{ij},$

where $N(i)$ are the triangles incident on vertex $i$ and $θ_{ij}$ is the angle
at vertex $i$ in triangle $j$ [][#meyer_2003].

Just like the continuous analog, our discrete Gaussian curvature reveals
elliptic, hyperbolic and parabolic vertices on the domain, as demonstrated in [Example 202](202GaussianCurvature/main.cpp).

![The `GaussianCurvature` example computes discrete Gaussian curvature and
visualizes it in pseudocolor.](images/bumpy-gaussian-curvature.jpg)

## Curvature directions
The two principal curvatures $(k_1,k_2)$ at a point on a surface measure how
much the surface bends in different directions. The directions of maximum and
minimum (signed) bending are called principal directions and are always
orthogonal.

Mean curvature is defined as the average of principal curvatures:

 $H = \frac{1}{2}(k_1 + k_2).$

One way to extract mean curvature is by examining the Laplace-Beltrami operator
applied to the surface positions. The result is a so-called mean-curvature
normal:

  $-\Delta \mathbf{x} = H \mathbf{n}.$

It is easy to compute this on a discrete triangle mesh in libigl using the
cotangent Laplace-Beltrami operator [][#meyer_2003].

```cpp
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
...
MatrixXd HN;
SparseMatrix<double> L,M,Minv;
igl::cotmatrix(V,F,L);
igl::massmatrix(V,F,igl::MASSMATRIX_TYPE_VORONOI,M);
igl::invert_diag(M,Minv);
HN = -Minv*(L*V);
H = HN.rowwise().norm(); //up to sign
```

Combined with the angle defect definition of discrete Gaussian curvature, one
can define principal curvatures and use least squares fitting to find
directions [][#meyer_2003].

Alternatively, a robust method for determining principal curvatures is via
quadric fitting [][#panozzo_2010]. In the neighborhood around every vertex, a
best-fit quadric is found and principal curvature values and directions are
analytically computed on this quadric ([Example
203](203_curvatureDirections/main.cpp)).

![The `CurvatureDirections` example computes principal curvatures via quadric
fitting and visualizes mean curvature in pseudocolor and principal directions
with a cross field.](images/fertility-principal-curvature.jpg)

## Gradient
Scalar functions on a surface can be discretized as a piecewise linear function
with values defined at each mesh vertex:

 $f(\mathbf{x}) \approx \sum\limits_{i=1}^n \phi_i(\mathbf{x})\, f_i,$

where $\phi_i$ is a piecewise linear hat function defined by the mesh so that
for each triangle $\phi_i$ is _the_ linear function which is one only at
vertex $i$ and zero at the other corners.

![Hat function $\phi_i$ is one at vertex $i$, zero at all other vertices, and
linear on incident triangles.](images/hat-function.jpg)

Thus gradients of such piecewise linear functions are simply sums of gradients
of the hat functions:

 $\nabla f(\mathbf{x}) \approx
 \nabla \sum\limits_{i=1}^n \phi_i(\mathbf{x})\, f_i =
 \sum\limits_{i=1}^n \nabla \phi_i(\mathbf{x})\, f_i.$

This reveals that the gradient is a linear function of the vector of $f_i$
values. Because the $\phi_i$ are linear in each triangle, their gradients are
_constant_ in each triangle. Thus our discrete gradient operator can be written
as a matrix multiplication taking vertex values to triangle values:

 $\nabla f \approx \mathbf{G}\,\mathbf{f},$

where $\mathbf{f}$ is $n\times 1$ and $\mathbf{G}$ is an $md\times n$ sparse
matrix. This matrix $\mathbf{G}$ can be derived geometrically, e.g.
[ch. 2][#jacobson_thesis_2013].
Libigl's `grad` function computes $\mathbf{G}$ for
triangle and tetrahedral meshes ([Example 204](204_Gradient/main.cpp)):

![The `Gradient` example computes gradients of an input function on a mesh and
visualizes the vector field.](images/cheburashka-gradient.jpg)

## Laplacian

The discrete Laplacian is an essential geometry processing tool. Many
interpretations and flavors of the Laplace and Laplace-Beltrami operator exist.

In open Euclidean space, the _Laplace_ operator is the usual divergence of
gradient (or equivalently the Laplacian of a function is the trace of its
Hessian):

 $\Delta f =
 \frac{\partial^2 f}{\partial x^2} +
 \frac{\partial^2 f}{\partial y^2} +
 \frac{\partial^2 f}{\partial z^2}.$

The _Laplace-Beltrami_ operator generalizes this to surfaces.

When considering piecewise-linear functions on a triangle mesh, a discrete
Laplacian may be derived in a variety of ways. The most popular in geometry
processing is the so-called ``cotangent Laplacian'' $\mathbf{L}$, arising
simultaneously from FEM, DEC and applying divergence theorem to vertex
one-rings. As a linear operator taking vertex values to vertex values, the
Laplacian $\mathbf{L}$ is a $n\times n$ matrix with elements:

$L_{ij} = \begin{cases}j \in N(i) &\cot \alpha_{ij} + \cot \beta_{ij},\\
j \notin N(i) & 0,\\
i = j & -\sum\limits_{k\neq i} L_{ik},
\end{cases}$

where $N(i)$ are the vertices adjacent to (neighboring) vertex $i$, and
$\alpha_{ij},\beta_{ij}$ are the angles opposite to edge ${ij}$.
This formula leads to a typical half-edge style implementation for
constructing $\mathbf{L}$:

```cpp
for(int i : vertices)
{
  for(int j : one_ring(i))
  {
    for(int k : triangle_on_edge(i,j))
    {
      L(i,j) = cot(angle(i,j,k));
      L(i,i) -= cot(angle(i,j,k));
    }
  }
}
```

Similarly as before, it may seem to loop over one-rings without having an half-edge data structure. However, this is not the case, since the Laplacian may be built by summing together contributions for each triangle, much in spirit with its FEM discretization
of the Dirichlet energy (sum of squared gradients):

```cpp
for(triangle t : triangles)
{
  for(edge i,j : t)
  {
    L(i,j) += cot(angle(i,j,k));
    L(j,i) += cot(angle(i,j,k));
    L(i,i) -= cot(angle(i,j,k));
    L(j,j) -= cot(angle(i,j,k));
  }
}
```

Libigl implements discrete "cotangent" Laplacians for triangles meshes and
tetrahedral meshes, building both with fast geometric rules rather than "by the
book" FEM construction which involves many (small) matrix inversions, cf.
[#sharf_2007][].

The operator applied to mesh vertex positions amounts to smoothing by _flowing_
the surface along the mean curvature normal direction ([Example 205](205_Laplacian/main.cpp)). Note that this is equivalent to minimizing surface area.

![The `Laplacian` example computes conformalized mean curvature flow using the
cotangent Laplacian [#kazhdan_2012][].](images/cow-curvature-flow.jpg)

### Mass matrix
The mass matrix $\mathbf{M}$ is another $n \times n$ matrix which takes vertex
values to vertex values. From an FEM point of view, it is a discretization of
the inner-product: it accounts for the area around each vertex. Consequently,
$\mathbf{M}$ is often a diagonal matrix, such that $M_{ii}$ is the barycentric
or voronoi area around vertex $i$ in the mesh [#meyer_2003][]. The inverse of
this matrix is also very useful as it transforms integrated quantities into
point-wise quantities, e.g.:

 $\Delta f \approx \mathbf{M}^{-1} \mathbf{L} \mathbf{f}.$

In general, when encountering squared quantities integrated over the surface,
the mass matrix will be used as the discretization of the inner product when
sampling function values at vertices:

 $\int_S x\, y\ dA \approx \mathbf{x}^T\mathbf{M}\,\mathbf{y}.$

An alternative mass matrix $\mathbf{T}$ is a $md \times md$ matrix which takes
triangle vector values to triangle vector values. This matrix represents an
inner-product accounting for the area associated with each triangle (i.e. the
triangles true area).

### Alternative construction of Laplacian

An alternative construction of the discrete cotangent Laplacian is by
"squaring" the discrete gradient operator. This may be derived by applying
Green's identity (ignoring boundary conditions for the moment):

  $\int_S \|\nabla f\|^2 dA = \int_S f \Delta f dA$

Or in matrix form which is immediately translatable to code:

  $\mathbf{f}^T \mathbf{G}^T \mathbf{T} \mathbf{G} \mathbf{f} =
  \mathbf{f}^T \mathbf{M} \mathbf{M}^{-1} \mathbf{L} \mathbf{f} =
  \mathbf{f}^T \mathbf{L} \mathbf{f}.$

So we have that $\mathbf{L} = \mathbf{G}^T \mathbf{T} \mathbf{G}$. This also
hints that we may consider $\mathbf{G}^T$ as a discrete _divergence_ operator,
since the Laplacian is the divergence of the gradient. Naturally, $\mathbf{G}^T$ is
a $n \times md$ sparse matrix which takes vector values stored at triangle faces
to scalar divergence values at vertices.

# Chapter 3: Matrices and linear algebra
Libigl relies heavily on the Eigen library for dense and sparse linear algebra
routines. Besides geometry processing routines, libigl has linear algebra
routines which bootstrap Eigen and make it feel even more similar to a high-level
algebra library such as Matlab.

## Slice
A very familiar and powerful routine in Matlab is array slicing. This allows
reading from or writing to a possibly non-contiguous sub-matrix. Let's consider
the Matlab code:

```matlab
B = A(R,C);
```

If `A` is a $m \times n$ matrix and `R` is a $j$-long list of row-indices
(between 1 and $m$) and `C` is a $k$-long list of column-indices, then as a
result `B` will be a $j \times k$ matrix drawing elements from `A` according to
`R` and `C`. In libigl, the same functionality is provided by the `slice`
function ([Example 301](301_Slice/main.cpp)):

```cpp
VectorXi R,C;
MatrixXd A,B;
...
igl::slice(A,R,C,B);
```

Note that `A` and `B` could also be sparse matrices.

Similarly, consider the Matlab code:

```matlab
A(R,C) = B;
```

Now, the selection is on the left-hand side so the $j \times k$ matrix  `B` is
being _written into_ the submatrix of `A` determined by `R` and `C`. This
functionality is provided in libigl using `slice_into`:

```cpp
igl::slice_into(B,R,C,A);
```

![The example `Slice` shows how to use `igl::slice` to change the colors for
triangles on a mesh.](images/decimated-knight-slice-color.jpg)

## Sort

Matlab and other higher-level languages make it very easy to extract indices of
sorting and comparison routines. For example in Matlab, one can write:

```matlab
[Y,I] = sort(X,1,'ascend');
```

so if `X` is a $m \times n$ matrix then `Y` will also be an $m \times n$ matrix
with entries sorted along dimension `1` in `'ascend'`ing order. The second
output `I` is a $m \times n$ matrix of indices such that `Y(i,j) =
X(I(i,j),j);`. That is, `I` reveals how `X` is sorted into `Y`.

This same functionality is supported in libigl:

```cpp
igl::sort(X,1,true,Y,I);
```

Similarly, sorting entire rows can be accomplished in Matlab using:

```matlab
[Y,I] = sortrows(X,'ascend');
```

where now `I` is a $m$ vector of indices such that `Y = X(I,:)`.

In libigl, this is supported with

```cpp
igl::sortrows(X,true,Y,I);
```
where again `I` reveals the index of sort so that it can be reproduced with
`igl::slice(X,I,1,Y)`.

Analogous functions are available in libigl for: `max`, `min`, and `unique`.

![The example `Sort` shows how to use `igl::sortrows` to
pseudocolor triangles according to their barycenters' sorted
order ([Example 302](302_Sort/main.cpp)).](images/decimated-knight-sort-color.jpg)


### Other Matlab-style functions
Libigl implements a variety of other routines with the same api and
functionality as common Matlab functions.

| Name                     | Description                                                                         |
| :----------------------- | :---------------------------------------------------------------------------------- |
| `igl::any_of`            | Whether any elements are non-zero (true)                                            |
| `igl::cat`               | Concatenate two matrices (especially useful for dealing with Eigen sparse matrices) |
| `igl::ceil`              | Round entries up to nearest integer |
| `igl::cumsum`            | Cumulative sum of matrix elements |
| `igl::colon`             | Act like Matlab's `:`, similar to Eigen's `LinSpaced` |
| `igl::cross`             | Cross product per-row |
| `igl::dot`               | dot product per-row |
| `igl::find`              | Find subscripts of non-zero entries |
| `igl::floot`             | Round entries down to nearest integer |
| `igl::histc`             | Counting occurrences for building a histogram |
| `igl::hsv_to_rgb`        | Convert HSV colors to RGB (cf. Matlab's `hsv2rgb`) |
| `igl::intersect`         | Set intersection of matrix elements. |
| `igl::jet`               | Quantized colors along the rainbow. |
| `igl::kronecker_product` | Compare to Matlab's `kronprod` |
| `igl::median`            | Compute the median per column |
| `igl::mode`              | Compute the mode per column |
| `igl::orth`              | Orthogonalization of a basis |
| `igl::setdiff`           | Set difference of matrix elements |
| `igl::speye`             | Identity as sparse matrix |

## Laplace equation
A common linear system in geometry processing is the Laplace equation:

 $∆z = 0$

subject to some boundary conditions, for example Dirichlet boundary conditions
(fixed value):

 $\left.z\right|_{\partial{S}} = z_{bc}$

In the discrete setting, the linear system can be written as:

 $\mathbf{L} \mathbf{z} = \mathbf{0}$

where $\mathbf{L}$ is the $n \times n$ discrete Laplacian and $\mathbf{z}$ is a
vector of per-vertex values. Most of $\mathbf{z}$ correspond to interior
vertices and are unknown, but some of $\mathbf{z}$ represent values at boundary
vertices. Their values are known so we may move their corresponding terms to
the right-hand side.

Conceptually, this is very easy if we have sorted $\mathbf{z}$ so that interior
vertices come first and then boundary vertices:

 $$\left(\begin{array}{cc}
 \mathbf{L}_{in,in} & \mathbf{L}_{in,b}\\
 \mathbf{L}_{b,in} & \mathbf{L}_{b,b}\end{array}\right)
 \left(\begin{array}{c}
 \mathbf{z}_{in}\\
 \mathbf{z}_{b}\end{array}\right) =
 \left(\begin{array}{c}
 \mathbf{0}_{in}\\
 \mathbf{z}_{bc}\end{array}\right)$$

The bottom block of equations is no longer meaningful so we'll only consider
the top block:

 $$\left(\begin{array}{cc}
 \mathbf{L}_{in,in} & \mathbf{L}_{in,b}\end{array}\right)
 \left(\begin{array}{c}
 \mathbf{z}_{in}\\
 \mathbf{z}_{b}\end{array}\right) =
 \mathbf{0}_{in}$$

We can move the known values to the right-hand side:

 $$\mathbf{L}_{in,in}
 \mathbf{z}_{in} = -
 \mathbf{L}_{in,b}
 \mathbf{z}_{b}$$

Finally we can solve this equation for the unknown values at interior vertices
$\mathbf{z}_{in}$.

However, our vertices will often not be sorted in this way. One option would be to sort `V`,
then proceed as above and then _unsort_ the solution `Z` to match `V`. However,
this solution is not very general.

With array slicing no explicit sort is needed. Instead we can _slice-out_
submatrix blocks ($\mathbf{L}_{in,in}$, $\mathbf{L}_{in,b}$, etc.) and follow
the linear algebra above directly. Then we can slice the solution _into_ the
rows of `Z` corresponding to the interior vertices ([Example 303](303_LaplaceEquation/main.cpp)).

![The `LaplaceEquation` example solves a Laplace equation with Dirichlet
boundary conditions.](images/camelhead-laplace-equation.jpg)

### Quadratic energy minimization

The same Laplace equation may be equivalently derived by minimizing Dirichlet
energy subject to the same boundary conditions:

 $\mathop{\text{minimize }}_z \frac{1}{2}\int\limits_S \|\nabla z\|^2 dA$

On our discrete mesh, recall that this becomes

 $\mathop{\text{minimize }}_\mathbf{z}  \frac{1}{2}\mathbf{z}^T \mathbf{G}^T \mathbf{D}
 \mathbf{G} \mathbf{z} \rightarrow \mathop{\text{minimize }}_\mathbf{z} \mathbf{z}^T \mathbf{L} \mathbf{z}$

The general problem of minimizing some energy over a mesh subject to fixed
value boundary conditions is so wide spread that libigl has a dedicated api for
solving such systems.

Let us consider a general quadratic minimization problem subject to different
common constraints:

 $$\mathop{\text{minimize }}_\mathbf{z}  \frac{1}{2}\mathbf{z}^T \mathbf{Q} \mathbf{z} +
 \mathbf{z}^T \mathbf{B} + \text{constant},$$

 subject to

 $$\mathbf{z}_b = \mathbf{z}_{bc} \text{ and } \mathbf{A}_{eq} \mathbf{z} =
 \mathbf{B}_{eq},$$

where

  - $\mathbf{Q}$ is a (usually sparse) $n \times n$ positive semi-definite
    matrix of quadratic coefficients (Hessian),
  - $\mathbf{B}$ is a $n \times 1$ vector of linear coefficients,
  - $\mathbf{z}_b$ is a $|b| \times 1$ portion of
$\mathbf{z}$ corresponding to boundary or _fixed_ vertices,
  - $\mathbf{z}_{bc}$ is a $|b| \times 1$ vector of known values corresponding to
    $\mathbf{z}_b$,
  - $\mathbf{A}_{eq}$ is a (usually sparse) $m \times n$ matrix of linear
    equality constraint coefficients (one row per constraint), and
  - $\mathbf{B}_{eq}$ is a $m \times 1$ vector of linear equality constraint
    right-hand side values.

This specification is overly general as we could write $\mathbf{z}_b =
\mathbf{z}_{bc}$ as rows of $\mathbf{A}_{eq} \mathbf{z} =
\mathbf{B}_{eq}$, but these fixed value constraints appear so often that they
merit a dedicated place in the API.

In libigl, solving such quadratic optimization problems is split into two
routines: precomputation and solve. Precomputation only depends on the
quadratic coefficients, known value indices and linear constraint coefficients:

```cpp
igl::min_quad_with_fixed_data mqwf;
igl::min_quad_with_fixed_precompute(Q,b,Aeq,true,mqwf);
```

The output is a struct `mqwf` which contains the system matrix factorization
and is used during solving with arbitrary linear terms, known values, and
constraint in the right-hand sides:

```cpp
igl::min_quad_with_fixed_solve(mqwf,B,bc,Beq,Z);
```

The output `Z` is a $n \times 1$ vector of solutions with fixed values
correctly placed to match the mesh vertices `V`.

## Linear equality constraints
We saw above that `min_quad_with_fixed_*` in libigl provides a compact way to
solve general quadratic programs. Let's consider another example, this time
with active linear equality constraints. Specifically let's solve the
`bi-Laplace equation` or equivalently minimize the Laplace energy:

 $$\Delta^2 z = 0 \leftrightarrow \mathop{\text{minimize }}\limits_z \frac{1}{2}
 \int\limits_S (\Delta z)^2 dA$$

subject to fixed value constraints and a linear equality constraint:

 $z_{a} = 1, z_{b} = -1$ and $z_{c} = z_{d}$.

Notice that we can rewrite the last constraint in the familiar form from above:

 $z_{c} - z_{d} = 0.$

Now we can assembly `Aeq` as a $1 \times n$ sparse matrix with a coefficient
$1$ in the column corresponding to vertex $c$ and a $-1$ at $d$. The right-hand
side `Beq` is simply zero.

Internally, `min_quad_with_fixed_*` solves using the Lagrange Multiplier
method. This method adds additional variables for each linear constraint (in
general a $m \times 1$ vector of variables $\lambda$) and then solves the
saddle problem:

 $$\mathop{\text{find saddle }}_{\mathbf{z},\lambda}\, \frac{1}{2}\mathbf{z}^T \mathbf{Q} \mathbf{z} +
  \mathbf{z}^T \mathbf{B} + \text{constant} + \lambda^T\left(\mathbf{A}_{eq}
 \mathbf{z} - \mathbf{B}_{eq}\right)$$

This can be rewritten in a more familiar form by stacking $\mathbf{z}$ and
$\lambda$ into one $(m+n) \times 1$ vector of unknowns:

 $$\mathop{\text{find saddle }}_{\mathbf{z},\lambda}\,
 \frac{1}{2}
 \left(
  \mathbf{z}^T
  \lambda^T
 \right)
 \left(
  \begin{array}{cc}
  \mathbf{Q}      & \mathbf{A}_{eq}^T\\
  \mathbf{A}_{eq} & 0
  \end{array}
 \right)
 \left(
  \begin{array}{c}
  \mathbf{z}\\
  \lambda
  \end{array}
 \right) +
 \left(
  \mathbf{z}^T
  \lambda^T
 \right)
 \left(
  \begin{array}{c}
  \mathbf{B}\\
  -\mathbf{B}_{eq}
  \end{array}
  \right)
  + \text{constant}$$

Differentiating with respect to $\left( \mathbf{z}^T \lambda^T \right)$ reveals
a linear system and we can solve for $\mathbf{z}$ and $\lambda$. The only
difference from the straight quadratic _minimization_ system, is that this
saddle problem system will not be positive definite. Thus, we must use a
different factorization technique (LDLT rather than LLT): libigl's
`min_quad_with_fixed_precompute` automatically chooses the correct solver in
the presence of linear equality constraints ([Example 304](304_LinearEqualityConstraints/main.cpp)).

![The example `LinearEqualityConstraints` first solves with just fixed value
constraints (left: 1 and -1 on the left hand and foot respectively), then
solves with an additional linear equality constraint (right: points on right
hand and foot constrained to be equal).](images/cheburashka-biharmonic-leq.jpg)

## Quadratic programming

We can generalize the quadratic optimization in the previous section even more
by allowing inequality constraints. Specifically box constraints (lower and
upper bounds):

 $\mathbf{l} \le \mathbf{z} \le \mathbf{u},$

where $\mathbf{l},\mathbf{u}$ are $n \times 1$ vectors of lower and upper
bounds
and general linear inequality constraints:

 $\mathbf{A}_{ieq} \mathbf{z} \le \mathbf{B}_{ieq},$

where $\mathbf{A}_{ieq}$ is a $k \times n$ matrix of linear coefficients and
$\mathbf{B}_{ieq}$ is a $k \times 1$ matrix of constraint right-hand sides.

Again, we are overly general as the box constraints could be written as
rows of the linear inequality constraints, but bounds appear frequently enough
to merit a dedicated api.

Libigl implements its own active set routine for solving _quadratric programs_
(QPs). This algorithm works by iteratively "activating" violated inequality
constraints by enforcing them as equalities and "deactivating" constraints
which are no longer needed.

After deciding which constraints are active at each iteration, the problem
reduces to a quadratic minimization subject to linear _equality_ constraints,
and the method from the previous section is invoked. This is repeated until convergence.

Currently the implementation is efficient for box constraints and sparse
non-overlapping linear inequality constraints.

Unlike alternative interior-point methods, the active set method benefits from
a warm-start (initial guess for the solution vector $\mathbf{z}$).

```cpp
igl::active_set_params as;
// Z is optional initial guess and output
igl::active_set(Q,B,b,bc,Aeq,Beq,Aieq,Bieq,lx,ux,as,Z);
```

![ [Example 305](305_QuadraticProgramming/main.cpp) uses an active set solver to optimize
discrete biharmonic kernels [#rustamov_2011][] at multiple scales
.](images/cheburashka-multiscale-biharmonic-kernels.jpg)

# Chapter 4: Shape deformation
Modern mesh-based shape deformation methods satisfy user deformation
constraints at handles (selected vertices or regions on the mesh) and propagate
these handle deformations to the rest of shape _smoothly_ and _without removing
or distorting details_. Libigl provides implementations of a variety of
state-of-the-art deformation techniques, ranging from quadratic mesh-based
energy minimizers, to skinning methods, to non-linear elasticity-inspired
techniques.

## Biharmonic deformation
The period of research between 2000 and 2010 produced a collection of
techniques that cast the problem of handle-based shape deformation as a
quadratic energy minimization problem or equivalently the solution to a linear
partial differential equation.

There are many flavors of these techniques, but a prototypical subset are those
that consider solutions to the bi-Laplace equation, that is a biharmonic
function [#botsch_2004][]. This fourth-order PDE provides sufficient
flexibility in boundary conditions to ensure $C^1$ continuity at handle
constraints (in the limit under refinement) [#jacobson_mixed_2010][].

### Biharmonic surfaces
Let us first begin our discussion of biharmonic _deformation_, by considering
biharmonic _surfaces_. We will casually define biharmonic surfaces as surface
whose _position functions_ are biharmonic with respect to some initial
parameterization:

 $\Delta^2 \mathbf{x}' = 0$

and subject to some handle constraints, conceptualized as "boundary
conditions":

 $\mathbf{x}'_{b} = \mathbf{x}_{bc}.$

where $\mathbf{x}'$ is the unknown 3D position of a point on the surface. So we
are asking that the bi-Laplacian of each of spatial coordinate function to be
zero.

In libigl, one can solve a biharmonic problem with `igl::harmonic`
and setting $k=2$ (_bi_-harmonic):

```cpp
// U_bc contains deformation of boundary vertices b
igl::harmonic(V,F,b,U_bc,2,U);
```

This produces a smooth surface that interpolates the handle constraints, but all
original details on the surface will be _smoothed away_. Most obviously, if the
original surface is not already biharmonic, then giving all handles the
identity deformation (keeping them at their rest positions) will **not**
reproduce the original surface. Rather, the result will be the biharmonic
surface that does interpolate those handle positions.

Thus, we may conclude that this is not an intuitive technique for shape
deformation.

### Biharmonic deformation fields
Now we know that one useful property for a deformation technique is "rest pose
reproduction": applying no deformation to the handles should apply no
deformation to the shape.

To guarantee this by construction we can work with _deformation fields_ (ie.
displacements)
$\mathbf{d}$ rather
than directly with positions $\mathbf{x}$. Then the deformed positions can be
recovered as

 $\mathbf{x}' = \mathbf{x}+\mathbf{d}.$

A smooth deformation field $\mathbf{d}$ which interpolates the deformation
fields of the handle constraints will impose a smooth deformed shape
$\mathbf{x}'$. Naturally, we consider _biharmonic deformation fields_:

 $\Delta^2 \mathbf{d} = 0$

subject to the same handle constraints, but rewritten in terms of their implied
deformation field at the boundary (handles):

 $\mathbf{d}_b = \mathbf{x}_{bc} - \mathbf{x}_b.$

Again we can use `igl::harmonic` with $k=2$, but this time solve for the
deformation field and then recover the deformed positions:

```cpp
// U_bc contains deformation of boundary vertices b
D_bc = U_bc - igl::slice(V,b,1);
igl::harmonic(V,F,b,D_bc,2,D);
U = V+D;
```

![The [BiharmonicDeformation](401_BiharmonicDeformation/main.cpp) example deforms a statue's head as a _biharmonic
surface_ (top) and using a _biharmonic displacements_
(bottom).](images/max-biharmonic.jpg)

#### Relationship to "differential coordinates" and Laplacian surface editing
Biharmonic functions (whether positions or displacements) are solutions to the
bi-Laplace equation, but also minimizers of the "Laplacian energy". For
example, for displacements $\mathbf{d}$, the energy reads

 $\int\limits_S \|\Delta \mathbf{d}\|^2 dA,$

where we define $\Delta \mathbf{d}$ to simply apply the Laplacian
coordinate-wise.

By linearity of the Laplace(-Beltrami) operator we can reexpress this energy in
terms of the original positions $\mathbf{x}$ and the unknown positions
$\mathbf{x}' = \mathbf{x} - \mathbf{d}$:

 $\int\limits_S \|\Delta (\mathbf{x}' - \mathbf{x})\|^2 dA = \int\limits_S
 \|\Delta \mathbf{x}' - \Delta \mathbf{x})\|^2 dA.$

In the early work of Sorkine et al., the quantities $\Delta \mathbf{x}'$ and
$\Delta \mathbf{x}$ were dubbed "differential coordinates" [#sorkine_2004][].
Their deformations (without linearized rotations) is thus equivalent to
biharmonic deformation fields.

## Polyharmonic deformation
We can generalize biharmonic deformation by considering different powers of
the Laplacian, resulting in a series of PDEs of the form:

 $\Delta^k \mathbf{d} = 0.$

with $k\in{1,2,3,\dots}$. The choice of $k$ determines the level of continuity
at the handles. In particular, $k=1$ implies $C^0$ at the boundary, $k=2$
implies $C^1$, $k=3$ implies $C^2$ and in general $k$ implies $C^{k-1}$.

```cpp
int k = 2;// or 1,3,4,...
igl::harmonic(V,F,b,bc,k,Z);
```

![The [PolyharmonicDeformation](402_PolyharmonicDeformation/main.cpp) example deforms a flat domain (left) into a bump as a
solution to various $k$-harmonic PDEs.](images/bump-k-harmonic.jpg)

## Bounded biharmonic weights
In computer animation, shape deformation is often referred to as "skinning".
Constraints are posed as relative rotations of internal rigid "bones" inside a
character. The deformation method, or skinning method, determines how the
surface of the character (i.e. its skin) should move as a function of the bone
rotations.

The most popular technique is linear blend skinning. Each point on the shape
computes its new location as a linear combination of bone transformations:

 $\mathbf{x}' = \sum\limits_{i = 1}^m w_i(\mathbf{x}) \mathbf{T}_i
 \left(\begin{array}{c}\mathbf{x}_i\\1\end{array}\right),$

where $w_i(\mathbf{x})$ is the scalar _weight function_ of the ith bone evaluated at
$\mathbf{x}$ and $\mathbf{T}_i$ is the bone transformation as a $4 \times 3$
matrix.

This formula is embarassingly parallel (computation at one point does not
depend on shared data need by computation at another point). It is often
implemented as a vertex shader. The weights and rest positions for each vertex
are sent as vertex shader _attributes_ and bone transformations are sent as
_uniforms_. Then vertices are transformed within the vertex shader, just in
time for rendering.

As the skinning formula is linear (hence its name), we can write it as matrix
multiplication:

 $\mathbf{X}' = \mathbf{M} \mathbf{T},$

where $\mathbf{X}'$ is $n \times 3$ stack of deformed positions as row
vectors, $\mathbf{M}$ is a $n \times m\cdot dim$ matrix containing weights and
rest positions and $\mathbf{T}$ is a $m\cdot (dim+1) \times dim$ stack of
transposed bone transformations.

Traditionally, the weight functions $w_j$ are painted manually by skilled
rigging professionals. Modern techniques now exist to compute weight functions
automatically given the shape and a description of the skeleton (or in general
any handle structure such as a cage, collection of points, selected regions,
etc.).

Bounded biharmonic weights are one such technique that casts weight computation
as a constrained optimization problem [#jacobson_2011][]. The weights enforce
smoothness by minimizing the familiar Laplacian energy:

 $\sum\limits_{i = 1}^m \int_S (\Delta w_i)^2 dA$

subject to constraints which enforce interpolation of handle constraints:

 $w_i(\mathbf{x}) = \begin{cases} 1 & \text{ if } \mathbf{x} \in H_i\\ 0 &
 \text{ otherwise } \end{cases},$

where $H_i$ is the ith handle, and constraints which enforce non-negativity,
parition of unity and encourage sparsity:

 $0\le w_i \le 1$ and $\sum\limits_{i=1}^m w_i = 1.$

This is a quadratic programming problem and libigl solves it using its active
set solver or by calling out to [Mosek](http://www.mosek.com).

![The example [BoundedBiharmonicWeights](403_BoundedBiharmonicWeights/main.cpp) computes weights for a tetrahedral
mesh given a skeleton (top) and then animates a linear blend skinning
deformation (bottom).](images/hand-bbw.jpg)

## Dual quaternion skinning
Even with high quality weights, linear blend skinning is limited. In
particular, it suffers from known artifacts stemming from blending rotations as
as matrices: a weight combination of rotation matrices is not necessarily a
rotation. Consider an equal blend between rotating by $-\pi/2$ and by $\pi/2$
about the $z$-axis. Intuitively one might expect to get the identity matrix,
but instead the blend is a degenerate matrix scaling the $x$ and $y$
coordinates by zero:

 $0.5\left(\begin{array}{ccc}0&-1&0\\1&0&0\\0&0&1\end{array}\right)+
 0.5\left(\begin{array}{ccc}0&1&0\\-1&0&0\\0&0&1\end{array}\right)=
 \left(\begin{array}{ccc}0&0&0\\0&0&0\\0&0&1\end{array}\right)$

In practice, this means the shape shrinks and collapses in regions where bone
weights overlap: near joints.

Dual quaternion skinning presents a solution [#kavan_2008]. This method
represents rigid transformations as a pair of unit quaternions,
$\hat{\mathbf{q}}$. The linear blend skinning formula is replaced with a
linear blend of dual quaternions:

 $\mathbf{x}' =
 \cfrac{\sum\limits_{i=1}^m w_i(\mathbf{x})\hat{\mathbf{q}_i}}
 {\left\|\sum\limits_{i=1}^m w_i(\mathbf{x})\hat{\mathbf{q}_i}\right\|}
 \mathbf{x},$

where $\hat{\mathbf{q}_i}$ is the dual quaternion representation of the rigid
transformation of bone $i$. The normalization forces the result of the linear
blending to again be a unit dual quaternion and thus also a rigid
transformation.

Like linear blend skinning, dual quaternion skinning is best performed in the
vertex shader. The only difference being that bone transformations are sent as
dual quaternions rather than affine transformation matrices.  Libigl supports
CPU-side dual quaternion skinning with the `igl::dqs` function, which takes a
more traditional representation of rigid transformations as input and
internally converts to the dual quaternion representation before blending:

```cpp
// vQ is a list of rotations as quaternions
// vT is a list of translations
igl::dqs(V,W,vQ,vT,U);
```

![The example [DualQuaternionSkinning](404_DualQuaternionSkinning/main.cpp) compares linear blend skinning (top) to dual
quaternion skinning (bottom), highlighting LBS's candy wrapper effect (middle)
and joint collapse (right).](images/arm-dqs.jpg)

## As-rigid-as-possible

Skinning and other linear methods for deformation are inherently limited.
Difficult arises especially when large rotations are imposed by the handle
constraints.

In the context of energy-minimization approaches, the problem stems from
comparing positions (our displacements) in the coordinate frame of the
undeformed shape. These quadratic energies are at best invariant to global
rotations of the entire shape, but not smoothly varying local rotations. Thus
linear techniques will not produce non-trivial bending and twisting.

Furthermore, when considering solid shapes (e.g. discretized with tetrahedral
meshes) linear methods struggle to maintain local volume, and they often suffer from
shrinking and bulging artifacts.

Non-linear deformation techniques present a solution to these problems.
They work by comparing the deformation of a mesh
vertex to its rest position _rotated_ to a new coordinate frame which best
matches the deformation. The non-linearity stems from the mutual dependence of
the deformation and the best-fit rotation. These techniques are often labeled
"as-rigid-as-possible" as they penalize the sum of all local deformations'
deviations from rotations.

To arrive at such an energy, let's consider a simple per-triangle energy:

 $E_\text{linear}(\mathbf{X}') = \sum\limits_{t \in T} a_t \sum\limits_{\{i,j\}
 \in t} w_{ij} \left\|
 \left(\mathbf{x}'_i - \mathbf{x}'_j\right) -
 \left(\mathbf{x}_i - \mathbf{x}_j\right)\right\|^2$

where $\mathbf{X}'$ are the mesh's unknown deformed vertex positions, $t$ is a
triangle in a list of triangles $T$, $a_t$ is the area of triangle $t$ and
$\{i,j\}$ is an edge in triangle $t$. Thus, this energy measures the norm of
change between an edge vector in the original mesh $\left(\mathbf{x}_i -
\mathbf{x}_j\right)$ and the unknown mesh $\left(\mathbf{x}'_i -
\mathbf{x}'_j\right)$.

This energy is **not** rotation invariant. If we rotate the mesh by 90 degrees
the change in edge vectors not aligned with the axis of rotation will be large,
despite the overall deformation being perfectly rigid.

So, the "as-rigid-as-possible" solution is to append auxiliary variables
$\mathbf{R}_t$
for each triangle $t$ which are constrained to be rotations. Then the energy is
rewritten, this time comparing deformed edge vectors to their rotated rest
counterparts:


 $E_\text{arap}(\mathbf{X}',\{\mathbf{R}_1,\dots,\mathbf{R}_{|T|}\}) = \sum\limits_{t \in T} a_t \sum\limits_{\{i,j\}
 \in t} w_{ij} \left\|
 \left(\mathbf{x}'_i - \mathbf{x}'_j\right)-
 \mathbf{R}_t\left(\mathbf{x}_i - \mathbf{x}_j\right)\right\|^2.$

The separation into the primary vertex position variables $\mathbf{X}'$ and the
rotations $\{\mathbf{R}_1,\dots,\mathbf{R}_{|T|}\}$ lead to strategy for
optimization, too. If the rotations $\{\mathbf{R}_1,\dots,\mathbf{R}_{|T|}\}$
are held fixed then the energy is quadratic in the remaining variables
$\mathbf{X}'$ and can be optimized by solving a (sparse) global linear system.
Alternatively, if $\mathbf{X}'$ are held fixed then each rotation is the
solution to a localized _Procrustes_ problem (found via $3 \times 3$ SVD or
polar decompostion). These two steps---local and global---each weakly decrease
the energy, thus we may safely iterate them until convergence.

The different flavors of "as-rigid-as-possible" depend on the dimension and
codimension of the domain and the edge-sets $T$. The proposed surface
manipulation technique by Sorkine and Alexa [#sorkine_2007][], considers $T$ to
be the set of sets of edges emanating from each vertex (spokes). Later, Chao et
al.  derived the relationship between "as-rigid-as-possible" mesh energies and
co-rotational elasticity considering 0-codimension elements as edge-sets:
triangles in 2D and tetrahedra in 3D [#chao_2010][]. They also showed how
Sorkine and Alexa's edge-sets are not a discretization of a continuous energy,
proposing instead edge-sets for surfaces containing all edges of elements
incident on a vertex (spokes and rims). They show that this amounts to
measuring bending, albeit in a discretization-dependent way.

Libigl, supports these common flavors. Selecting one is a matter of setting the
energy type before the precompuation phase:

```cpp
igl::ARAPData data;
arap_data.energy = igl::ARAP_ENERGY_TYPE_SPOKES;
//arap_data.energy = igl::ARAP_ENERGY_TYPE_SPOKES_AND_RIMS;
//arap_data.energy = igl::ARAP_ENERGY_TYPE_ELEMENTS; //triangles or tets
igl::arap_precomputation(V,F,dim,b,data);
```

Just like `igl::min_quad_with_fixed_*`, this precomputation phase only depends
on the mesh, fixed vertex indices `b` and the energy parameters. To solve with
certain constraints on the positions of vertices in `b`, we may call:

```cpp
igl::arap_solve(bc,data,U);
```

which uses `U` as an initial guess and then computes the solution into it.

Libigl's implementation of as-rigid-as-possible deformation takes advantage of
the highly optimized singular value decomposition code from McAdams et al.
[#mcadams_2011][] which leverages SSE intrinsics.

![The example [AsRigidAsPossible](405_AsRigidAsPossible/main.cpp) deforms a surface as if it were made of an
elastic material](images/decimated-knight-arap.jpg)

The concept of local rigidity will be revisited shortly in the context of
surface parameterization.

## Fast automatic skinning transformations

Non-linear optimization is, unsurprisingly, slower than its linear cousins. In
the case of the as-rigid-as-possible optimization, the bottleneck is typically
the large number of polar decompositions necessary to recover best fit
rotations for each edge-set (i.e. for each triangle, tetrahedron, or vertex
cell). Even if this code is optimized, the number of primary degrees of freedom
is tied to the discretization level, despite the deformations' low frequency
behavior.

This invites two routes toward fast non-linear optimization. First, is it
necessary (or even advantageous) to find so many best-fit rotations? Second,
can we reduce the degrees of freedom to better reflect the frequency of the
desired deformations.

Taken in turn, these optimizations culminate in a method which optimizes over
the space of linear blend skinning deformations spanned by high-quality weights
(i.e. manually painted ones or bounded biharmonic weights). This space is a
low-dimensional subspace of all possible mesh deformations, captured by writing
linear blend skinning in matrix form:

 $\mathbf{X}' = \mathbf{M}\mathbf{T}$

where the mesh vertex positions in the $n \times 3$ matrix $\mathbf{X}'$ are
replaced by a linear combination of a small number of degrees of freedom in the
$(3+1)m \times 3$ stack of transposed "handle" transformations. Swapping in
$\mathbf{M}\mathbf{T}$ for $\mathbf{X}'$ in the ARAP energies above immediately
sees performance gains during the global solve step as $m << n$.

The complexity of the local step---fitting rotations---is still bound
to the original mesh discretization. However, if the skinning is well behaved,
we can make the assumption that places on the shape with similar skinning
weights will deform similarly and thus imply similar best-fit rotations.
Therefore, we cluster edge-sets according to their representation in
_weight-space_: where a vertex $\mathbf{x}$ takes the coordinates
$[w_1(\mathbf{x}),w_2(\mathbf{x}),\dots,w_m(\mathbf{x})]$. The number of
clustered edge-sets show diminishing returns on the deformation quality so we
may choose a small number of clusters, proportional to the number of skinning
weight functions (rather than the number of discrete mesh vertices).

This proposed deformation model [#jacobson_2012][], can simultaneously be seen as a
fast, subspace optimization for ARAP and as an automatic method for finding
_the best_ skinning transformation degrees of freedom.

A variety of user interfaces are supported via linear equality constraints on
the skinning transformations associated with handles. To fix a transformation
entirely we simply add the constraint:

 $\left(\begin{array}{cccc}
 1 & 0 & 0 & 0\\
 0 & 1 & 0 & 0\\
 0 & 0 & 1 & 0\\
 0 & 0 & 0 & 1\end{array}\right)
 \mathbf{T}_i^T = \hat{\mathbf{T}}_i^T,$

where $\hat{\mathbf{T}}_i^T$ is the $(3+1) \times 3$ transposed fixed
transformation for handle $i$.

To fix only the origin of a handle, we add a constraint requiring the
transformation to interpolate a point in space (typically the centroid of all
points with $w_i = 1$:

 $\mathbf{c}'^T\mathbf{T}_i^T = \mathbf{c}^T,$

where $\mathbf{c}^T$ is the $1 \times (3+1)$ position of the point at rest in
transposed homogeneous coordinates, and $\mathbf{c}'^T$ the point given by the
user.

We can similarly fix just the linear part of the transformation at a handle,
freeing the translation component (producing a "chickenhead" effect):

 $\left(\begin{array}{cccc}
 1&0&0&0\\
 0&1&0&0\\
 0&0&1&0\end{array}\right)
 \mathbf{T}_i^T = \hat{\mathbf{L}}_i^T,$

where $\hat{\mathbf{L}}_i^T$ is the fixed $3 \times 3$ linear part of the
transformation at handle $i$.

And lastly we can allow the user to entirely _free_ the transformation's
degrees of freedom, delegating the optimization to find the best possible
values for all elements. To do this, we simply abstain from adding a
corresponding constraint.

### ARAP with grouped edge-sets

Being a subspace method, an immediate disadvantage is the reduced degrees of
freedom. This brings performance, but in some situations limits behavior too
much. In such cases one can use the skinning subspace to build an effective
clustering of rotation edge-sets for a traditional ARAP optimization: forgoing
the subspace substitution. This has an two-fold effect. The cost of the
rotation fitting, local step drastically reduces, and the deformations are
"regularized" according the clusters. From a high level point of view, if the clusters
are derived from skinning weights, then they will discourage bending,
especially along isolines of the weight functions.

In this light, we can think of the "spokes+rims" style surface ARAP as a (slight and
redundant) clustering of the per-triangle edge-sets.

![The example [FastAutomaticSkinningTransformations](406_FastAutomaticSkinningTransformations/main.cpp) compares a full (slow)
ARAP deformation on a detailed shape (left of middle), to ARAP with grouped
rotation edge sets (right of middle), to the very fast subpsace method
(right).](images/armadillo-fast.jpg)

# Chapter 5: Parametrization [500]

In computer graphics, we denote as surface parametrization a map from the
surface to \\(\mathbf{R}^2\\). It is usually encoded by a new set of 2D
coordinates for each vertex of the mesh (and possibly also by a new set of
faces in one to one correspondence with the faces of the original surface).
Note that
this definition is the *inverse* of the classical differential geometry
definition.

A parametrization has many applications, ranging from texture mapping to
surface remeshing. Many algorithms have been proposed, and they can be broadly
divided in four families:

1. **Single patch, fixed boundary**: these algorithm can parametrize a
disk-like part of the surface given fixed 2D positions for its boundary. These
algorithms are efficient and simple, but they usually produce high-distortion maps due to the fixed boundary.

2. **Single patch, free boundary:** these algorithms let the boundary
deform freely, greatly reducing the map distortion. Care should be taken to
prevent the border to self-intersect.

3. **Global parametrization**: these algorithms work on meshes with arbitrary
genus. They initially cut the mesh in multiple patches that can be separately parametrized. The generated maps are discontinuous on the cuts (often referred as *seams*).

4. **Global seamless parametrization**: these are global parametrization algorithm that hides the seams, making the parametrization "continuous", under specific assumptions that we will discuss later.

## Harmonic parametrization [501]

Harmonic parametrization [#eck_2005][] is a single patch, fixed boundary parametrization
algorithm that computes the 2D coordinates of the flattened mesh as two
harmonic functions.

The algorithm is divided in 3 steps:

1. Detect of the boundary vertices
2. Map the boundary vertices to a circle
3. Compute two harmonic functions (one for u and one for the v coordinate). The harmonic functions use the fixed vertices on the circle as boundary constraints.

The algorithm can be coded using libigl as follows:

```cpp
Eigen::VectorXi bnd;
igl::boundary_loop(V,F,bnd);

Eigen::MatrixXd bnd_uv;
igl::map_vertices_to_circle(V,bnd,bnd_uv);

igl::harmonic(V,F,bnd,bnd_uv,1,V_uv);
```

where `bnd` contains the indices of the boundary vertices, bnd_uv their position on the UV plane, and "1" denotes that we want to compute an harmonic function (2 will be for biharmonic, 3 for triharmonic, etc.). Note that each of the three
functions is designed to be reusable in other parametrization algorithms.

A UV parametrization can be visualized in the viewer with:

```cpp
viewer.data.set_uv(V_uv);
```

The UV coordinates are then used to apply a procedural checkerboard texture to the
mesh ([Example 501](501_HarmonicParam/main.cpp)).

![([Example 501](501_HarmonicParam/main.cpp)) Harmonic parametrization. (left)
mesh with texture, (right) UV parametrization with
texture](images/501_HarmonicParam.png)

## Least squares conformal maps [502]

Least squares conformal maps parametrization [#levy_2002][] minimizes the
conformal (angular) distortion of the parametrization. Differently from
harmonic parametrization, it does not need to have a fixed boundary.

LSCM minimizes the following energy:

\\[ E_{LSCM}(\mathbf{u},\mathbf{v}) = \int_X \frac{1}{2}| \nabla \mathbf{u}^{\perp} - \nabla \mathbf{v} |^2 dA \\]

which can be rewritten in matrix form as [#mullen_2008][]:

\\[ E_{LSCM}(\mathbf{u},\mathbf{v}) = \frac{1}{2} [\mathbf{u},\mathbf{v}]^t (L_c - 2A) [\mathbf{u},\mathbf{v}] \\]

where $L_c$ is the cotangent Laplacian matrix and A is a matrix such that \\(
[\mathbf{u},\mathbf{v}]^t A  [\mathbf{u},\mathbf{v}]$ is equal to the
[vector area](http://en.wikipedia.org/wiki/Vector_area) of the mesh.

Using libigl, this matrix energy can be written in a few lines of codes. The
cotangent matrix can be computed using `igl::cotmatrix`:

```cpp
SparseMatrix<double> L;
igl::cotmatrix(V,F,L);
```

Note that we want to apply the Laplacian matrix to the u and v coordinates at
the same time, thus we need to extend it taking the left
Kronecker product with a 2x2 identity matrix:

```cpp
SparseMatrix<double> L_flat;
igl::repdiag(L,2,L_flat);
```

The area matrix is computed with `igl::vector_area_matrix`:

```cpp
SparseMatrix<double> A;
igl::vector_area_matrix(F,A);
```

The final energy matrix is the sum of these two matrices. Note that in this
case we do not need to fix the boundary. To remove the null space of the energy and make the minimum unique, it is sufficient to fix two arbitrary
vertices to two arbitrary positions. The full source code is provided in [Example 502](502_LSCMParam/main.cpp).


![([Example 502](502_LSCMParam/main.cpp)) LSCM parametrization. (left) mesh
with texture, (right) UV parametrization](images/502_LSCMParam.png)

## As-rigid-as-possible parametrization [503]

As-rigid-as-possible parametrization [#liu_2008][] is a powerful single-patch,
non-linear algorithm to compute a parametrization that strives to preserve
distances (and thus angles). The idea is very similar to ARAP surface
deformation: each triangle is mapped to the plane trying to preserve its
original shape, up to a rigid rotation.

The algorithm can be implemented reusing the functions discussed in the
deformation chapter: `igl::arap_precomputation` and `igl::arap_solve`. The only
difference is that the optimization has to be done in 2D instead of 3D and that
we need to compute a starting point. While for 3D deformation the optimization
is bootstrapped with the original mesh, this is not the case for ARAP
parametrization since the starting point must be a 2D mesh. In [Example
503](503_ARAPParam/main.cpp), we initialize the optimization with harmonic
parametrization. Similarly to LSCM, the boundary is free to deform to minimize
the distortion.

![([Example 503](502_ARAPParam/main.cpp)) As-Rigid-As-Possible parametrization.
(left) mesh with texture, (right) UV parametrization with
texture](images/503_ARAPParam.png)

## N-rotationally symmetric tangent fields [504]

The design of tangent fields is a basic tool used to design guidance fields for
uniform quadrilateral and hexahedral remeshing. Libigl contains an
implementation of all the state-of-the-art algorithms to design N-RoSy fields
and their generalizations.

In libigl, tangent unit-length vector fields are piece-wise constant on the
faces of a triangle mesh, and they are described by one or more vectors per-face. The function

```cpp
igl::nrosy(V,F,b,bc,b_soft,b_soft_weight,bc_soft,N,0.5,
           output_field,output_singularities);
```

creates a smooth unit-length vector field (N=1) starting from a sparse set of
constrained faces, whose indices are listed in b and their constrained value is
specified in bc. The functions supports soft_constraints (b_soft,
b_soft_weight, bc_soft), and returns the interpolated field for each face of
the triangle mesh (output_field), plus the singularities of the field
(output_singularities).

![Design of a unit-length vector field](images/504_vector_field.png)

The singularities are vertices where the field vanishes (highlighted in red in
the figure above). `igl::nrosy` can also generate N-RoSy fields [#levy_2008][],
which are a generalization of vector fields where in every face the vector is
defined up to a constant rotation of $2\pi / N$. As can be observed in
the following figure, the singularities of the fields generated with different
N are of different types and they appear in different positions.

![Design of a 2-,4- and 9-RoSy field](images/504_nrosy_field.png)

We demonstrate how to call and plot N-RoSy fields in [Example
504](504_NRosyDesign/main.cpp), where the degree of the field can be change
pressing the number keys. `igl::nrosy` implements the algorithm proposed in
[#bommes_2009][]. N-RoSy fields can also be interpolated with the algorithm
proposed in [#knoppel_2013][], see Section [507] for more details
([igl::n_polyvector](../include/igl/n_polyvector.h)).

### Global, seamless integer-grid parametrization [505]

The previous parametrization methods were focusing on creating parametrizations
of surface patches aimed at texture mapping or baking of other surface
properties such as normals and high-frequency details. Global, seamless
parametrization aims at parametrizing complex shapes with a parametrization
that is aligned with a given set of directions for the purpose of surface
remeshing. In libigl, we provide a reference  implementation of the pipeline
proposed in the mixed integer quadrangulation paper [#bommes_2009][].

The first step involves the design of a 4-RoSy field (sometimes called *cross*
field) that describes the alignment of the edges of the desired quadrilateral
remeshing. The field constraints are usually manually specified or extracted
from the principal curvature directions. In [[Example
506](506_FrameField/main.cpp)], we simply fix one face in a random direction.

![Initial cross field prescribing the edge alignment.](images/505_MIQ_1.png)

### Combing and cutting

Given the cross field, we now want to cut the surface so that it becomes
homeomorphic to a disk. While this could be done directly on the cross-field, we
opt to perform this operation on its bisector field (a copy of the field
rotated by 45 degrees) since it is more stable and generic. Working on the
bisectors allow us to take as input generalized, non-orthogonal and non-unit
length cross fields.

We thus rotate the field,

![Bisector field.](images/505_MIQ_2.png)

and we remove the rotation ambiguity by assigning to each face a u and a v
direction. The assignment is done with a breadth-first search starting from a
random face.

![Combed bisector field.](images/505_MIQ_3.png)

You can imagine this process as combing an hairy surface: you will be able to
comb part of it, but at some point you will not be able to consistently comb
the entire surface ([Hairy ball
theorem](http://en.wikipedia.org/wiki/Hairy_ball_theorem)). The discontinuities
in the combing define the cut graph:

![Cut graph.](images/505_MIQ_4.png)

Finally, we rotate the combed field by 45 degrees to undo the initial degrees
rotation:

![Combed cross field.](images/505_MIQ_5.png)

The combed cross field can be seen as the ideal Jacobian of the parametrization
that will be computed in the next section.

### Poisson parametrization

The mesh is cut along the seams and a parametrization is computed trying to
find two scalar functions whose gradient matches the combed cross field
directions. This is a classical Poisson problem, that is solved minimizing the
following quadratic energy:

\\[ E(\mathbf{u},\mathbf{v}) = |\nabla \mathbf{u} - X_u|^2 + |\nabla \mathbf{v} - X_v|^2 \\]

where $X_u$ and $X_u$ denotes the combed cross field. Solving this
problem generates a parametrization whose u and v isolines are aligned with the
input cross field.

![Poisson parametrization.](images/505_MIQ_8.png)

We hide the seams by adding integer constraints to the Poisson problem
that align the isolines on both sides of each seam [#bommes_2009].

![Seamless Poisson parametrization.](images/505_MIQ_7.png)

Note that this parametrization can only be used for remeshing purposes, since
it contains many overlaps.

![Seamless Poisson parametrization (in 2D).](images/505_MIQ_6.png)

A quad mesh can be extracted from this parametrization using
[libQEx](https://github.com/hcebke/libQEx) (not included in libigl).
The full pipeline is implemented in [Example 505](505_MIQ/main.cpp).

## Anisotropic remeshing [506]

Anisotropic and non-uniform quad remeshing is important to concentrate the
elements in the regions with more details. It is possible to extend the MIQ
quad meshing framework to generate anisotropic quad meshes using a mesh
deformation approach [#panozzo_2014][].

The input of the anisotropic remeshing algorithm is a sparse set of constraints
that define the shape and scale of the desired quads. This can be encoded as a
frame field, which is a pair of non-orthogonal and non-unit length vectors. The
frame field can be interpolated by decomposing it in a 4-RoSy field and a
unique affine transformation. The two parts can then be interpolated
separately, using `igl::nrosy` for the cross field, and an harmonic interpolant
for the affine part.

![Interpolation of a frame field. Colors on the vectors denote the desired
scale. The red faces contains the frame field
constraints.](images/506_FrameField_1.png)

After the interpolation, the surface is warped to transform each frame into an
orthogonal and unit length cross (i.e. removing the scaling and skewness from
the frame). This deformation defines a new embedding (and a new metric) for the
surface.

![The surface is deformed to transform the frame field in a cross
field.](images/506_FrameField_2.png)

The deformed surface can the be isotropically remeshed using the MIQ algorithm
that has been presented in the previous section.

![The deformed surface is isotropically remeshed.](images/506_FrameField_3.png)

The UV coordinates of the deformed surface can then be used to transport the
parametrization to the original surface, where the isolines will trace a quad
mesh whose elements are similar to the shape prescribed in the input frame
field.

![The global parametrization is lifted to the original surface to create the
anisotropic quad meshing.](images/506_FrameField_4.png)

Our implementation ([Example 506](506_FrameField/main.cpp)) uses MIQ to
generate the UV parametrization, but other algorithms could be applied: the
only desiderata is that the generated quad mesh should be as isotropic as
possible.

## N-PolyVector fields [507]

N-RoSy vector fields can be further generalized to represent arbitrary
vector-sets, with arbitrary angles between them and with arbitrary lengths
[#diamanti_2014][].  This generalization is called  N-PolyVector field, and
libigl provides the function `igl::n_polyvector` to design them starting from a
sparse set of constraints ([Example 507](507_PolyVectorField/main.cpp)).

![Interpolation of a 6-PolyVector field (right) and a 12-PolyVector field from a sparse set of random constraints.](images/507_PolyVectorField.png)

The core idea is to represent the vector set as the roots of a complex
polynomial: The polynomial coefficients are then harmonically interpolated
leading to polynomials whose roots smoothly vary over the surface.

Globally optimal direction fields [#knoppel_2013][] are a special case of
PolyVector fields. If the constraints are taken from an N-RoSy field,
`igl::n_polyvector` generates a field that is equivalent, after normalization,
to a globally optimal direction field.

## Conjugate vector fields [508]

Two tangent vectors lying on a face of a triangle mesh are conjugate if

\\[ k_1 (u^T d_1)(v^T d_1) + k_2(u^T d_2)(v^T d_2) = 0. \\]

This condition is very important in architectural geometry: The faces of an
infinitely dense quad mesh whose edges are aligned with a conjugate field are
planar. Thus, a quad mesh whose edges follow a conjugate field  are easier to
planarize [#liu_2011].

Finding a conjugate vector field that satisfies given directional constraints
is a standard problem in architectural geometry, which can be tackled by
deforming a Poly-Vector field to the closest conjugate field.

This algorithm [#diamanti_2014] alternates a global step, which enforces
smoothness, with a local step, that projects the field on every face to the
closest conjugate field ([Example 508](508_ConjugateField/main.cpp)).

![A smooth 4-PolyVector field (left) is deformed to become a conjugate field
(right).](images/508_ConjugateField.png)

## Planarization [509]

A quad mesh can be transformed in a planar quad mesh with Shape-Up
[#bouaziz_2012], a local/global approach that uses the global step to enforce
surface continuity and the local step to enforce planarity.

[Example 509](509_Planarization/main.cpp) planarizes a quad mesh until it
satisfies a user-given planarity threshold.

![A non-planar quad mesh (left) is planarized using the libigl function
igl::palanarize (right). The colors represent the planarity of the
quads.](images/509_Planarization.png)

# Chapter 6: External libraries [600]

An additional positive side effect of using matrices as basic types is that it
is easy to exchange data between libigl and other softwares and libraries.

## State serialization [601]

Geometry processing applications often require a considerable amount of
computational time and/or manual input. Serializing the state of the application is a simple strategy to greatly increase the development efficiency. It allows to quickly start debugging just
before the crash happens, avoiding to wait for the precomputation to take place
every time and it also makes your experiments reproducible, allowing to quickly test algorithms variants on the same input data.

Serialization is often not considered in geometry processing due
to the extreme difficulty in serializing pointer-based data structured, such as
an half-edge data structure ([OpenMesh](http://openmesh.org), [CGAL](http://www.cgal.org)), or a pointer based indexed structure ([VCG](http://vcg.isti.cnr.it/~cignoni/newvcglib/html/)).

In libigl, serialization is much simpler, since the majority of the functions use basic types, and pointers are used in very rare cases (usually to interface
with external libraries). Libigl bundles a simple and self-contained XML serialization framework, that drastically reduces the overhead required to add
serialization to your applications.

Assume that the state of your application is a mesh and a set of
integer ids:

```cpp
class State : public igl::XMLSerialization
{
public:
  State() : XMLSerialization("dummy") {}

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  std::vector<int> ids;

  void InitSerialization()
  {
    xmlSerializer->Add(V  , "V");
    xmlSerializer->Add(F  , "F");
    xmlSerializer->Add(ids, "ids");
  }
};
```

Any class can be made serializable by inheriting from ``igl::XMLSerialization` and trivially implementing the `InitSerialization` method. The library can serialize all the basic `stl` types, all `Eigen` types and any class inheriting
from `igl::XMLSerialization`.

The state can be saved into an xml file with:

```cpp
igl::XMLSerializer serializer_save("601_Serialization");
serializer_save.Add(state,"State");
serializer_save.Save("temp.xml",true);
```

This code generates the following xml file (assuming `V` and `F` contains a simple mesh with two triangles, and `ids` contains the numbers 6 and 7):

```xml
<:::601_Serialization>
    <State>
        <V rows="4" cols="3" matrix="
0,0,0,
1,0,0,
1,1,1,
2,1,0"/>
        <F rows="2" cols="3" matrix="
0,1,2,
1,3,2"/>
        <ids size="2" vector_int="
6,7"/>
    </State>
</:::601_Serialization>
```

The xml file can be loaded in a similar way:

```cpp
State loaded_state;
igl::XMLSerializer serializer_load("601_Serialization");
serializer_load.Add(loaded_state,"State");
serializer_load.Load("temp.xml");
```

The serialization framework can also be used as a convenient interface to
provide parameters to command line applications, since the xml files can be
directly edited with a standard text editor.

The code snippets above are extracted from [Example
601](601_Serialization/main.cpp). We strongly suggest that you make the entire
state of your application always serializable since it will save you a lot of
troubles when you will be preparing figures for a scientific report. It is very
common to have to do small changes to figures, and being able to serialize the
entire state just before you take screenshots will save you many painful hours
before a submission deadline.

## Mixing Matlab code [602]

Libigl can be interfaced with Matlab to offload numerically heavy computation
to a Matlab script. The major advantage of this approach is that you will be
able to develop efficient and complex user-interfaces in C++, while exploting
the syntax and fast protototyping features of matlab. In particular, the use of
an external Matlab script in a libigl application allows to change the Matlab
code while the C++ application is running, greatly increasing coding
efficiency.

We demonstrate how to integrate Matlab in a libigl application in [Example
602](602_Matlab/main.cpp). The example uses Matlab to compute the
Eigenfunctions of the discrete Laplacian operator, relying on libigl for mesh
IO, visualization and for computing the Laplacian operator.

Libigl can connect to an existing instance of Matlab (or launching a new one on
Linux/MacOSX) using:

```cpp
igl::mlinit(&engine);
```

The cotangent Laplacian is computed using igl::cotmatrix and uploaded to the
Matlab workspace:

```cpp
igl::cotmatrix(V,F,L);
igl::mlsetmatrix(&engine,"L",L);
```

It is now possible to use any Matlab function on the data. For example, we can
see the sparsity pattern of L using spy:

```cpp
igl::mleval(&engine,"spy(L)");
```

![The Matlab spy function is called from a libigl-based
application.](images/602_Matlab_1.png)

The results of Matlab computations can be returned back to the C++ application

```cpp
igl::mleval(&engine,"[EV,~] = eigs(-L,10,'sm')");
igl::mlgetmatrix(&engine,"EV",EV);
```

and plotted using the libigl viewer.

![4 Eigenfunctions of the Laplacian plotted in the libigl
viewer.](images/602_Matlab_2.png)


### Saving a Matlab workspace
To aid debugging, libigl also supplies functions to write Matlab `.mat`
"Workspaces". This C++ snippet saves a mesh and it's sparse Laplacian matrix to
a file:

```cpp
igl::readOFF("../shared/fertility.off", V, F);
igl::cotmatrix(V,F,L);
igl::MatlabWorkspace mw;
mw.save(V,"V");
mw.save_index(F,"F");
mw.save(L,"L");
mw.write("fertility.mat");
```

Then this workspace can be loaded into a Matlab IDE:

```matlab
load fertility.mat
```

The `igl::MatlabWorkspace` depends on Matlab libraries to compile and run,
but---in contrast to the engine routines above---will avoid launching a Matlab
instance upon execution.

### Dumping Eigen matrices to copy and paste into Matlab
Eigen supplies a sophisticated API for printing its matrix types to the screen.
Libigl has wrapped up a particularly useful formatting which makes it simple to
copy standard output from a C++ program into a Matlab IDE. The code:

```cpp
igl::readOFF("../shared/2triangles.off", V, F);
igl::cotmatrix(V,F,L);
std::cout<<igl::matlab_format(V,"V")<<std::endl;
std::cout<<igl::matlab_format((F.array()+1).eval(),"F")<<std::endl;
std::cout<<igl::matlab_format(L,"L")<<std::endl;
```

produces the output:

```matlab
V = [
  0 0 0
  1 0 0
  1 1 1
  2 1 0
];
F = [
  1 2 3
  2 4 3
];
LIJV = [
1  1    -0.7071067811865476
2  1     0.7071067811865475
3  1  1.570092458683775e-16
1  2     0.7071067811865475
2  2     -1.638010440969447
3  2     0.6422285251880865
4  2     0.2886751345948129
1  3  1.570092458683775e-16
2  3     0.6422285251880865
3  3    -0.9309036597828995
4  3     0.2886751345948129
2  4     0.2886751345948129
3  4     0.2886751345948129
4  4    -0.5773502691896258
];
L = sparse(LIJV(:,1),LIJV(:,2),LIJV(:,3));
```

which is easily copied and pasted into Matlab for debugging, etc.

## Calling libigl functions from Matlab [603]

It is also possible to call libigl functions from matlab, compiling them as MEX
functions. This can be used to offload to C++ code the computationally
intensive parts of a Matlab application.

We provide a wrapper for `igl::readOBJ` in [Example 603](603_MEX/compileMEX.m).
We plan to provide wrappers for all our functions in the future, if you are
interested in this feature (or if you want to help implementing it) please let
us know.

## Triangulation of closed polygons [604]

The generation of high-quality triangle and tetrahedral meshes is a very common
task in geometry processing. We provide wrappers in libigl to
[triangle](http://www.cs.cmu.edu/~quake/triangle.html) and
[Tetgen](http://wias-berlin.de/software/tetgen/).

A triangle mesh with a given boundary can be created with:

```cpp
igl::triangulate(V,E,H,V2,F2,"a0.005q");
```

where `E` is a set of boundary edges (#E by 2), `H` is a set of 2D positions of
points contained in holes of the triangulation (#H by 2) and (`V2`,`F2`) is the
generated triangulation. Additional parameters can be passed to `triangle`, to
control the quality: `"a0.005q"` enforces a bound on the maximal area of the
triangles and a minimal angle of 20 degrees. In [Example
604](604_Triangle/main.m), the interior of a square (excluded a smaller square
in its interior) is triangulated.

![Triangulation of the interior of a polygon.](images/604_Triangle.png)

## Tetrahedralization of closed surfaces [605]

Similarly, the interior of a closed manifold surface can be tetrahedralized
using the function `igl::tetrahedralize` which wraps the Tetgen library ([Example
605](605_Tetgen/main.c)):

```cpp
igl::tetrahedralize(V,F,"pq1.414", TV,TT,TF);
```

![Tetrahedralization of the interior of a surface mesh.](images/605_Tetgen.png)

## Baking ambient occlusion [606]

[Ambient occlusion](http://en.wikipedia.org/wiki/Ambient_occlusion) is a
rendering technique used to calculate the exposure of each point in a surface
to ambient lighting. It is usually encoded as a scalar (normalized between 0
and 1) associated with the vertice of a mesh.

Formally, ambient occlusion is defined as:

\\[ A_p = \frac{1}{\pi} \int_\omega V_{p,\omega}(n \cdot \omega) d\omega \\]

where $V_{p,\omega}$ is the visibility function at  p, defined to be zero if p
is occluded in the direction $\omega$ and one otherwise, and $d\omega$ is the
infinitesimal solid angle step of the integration variable $\omega$.

The integral is usually approximated by casting rays in random directions
around each vertex. This approximation can be computed using the function:

```cpp
igl::ambient_occlusion(V,F,V_samples,N_samples,500,AO);
```

that given a scene described in `V` and `F`, computes the ambient occlusion of
the points in `V_samples` whose associated normals are `N_samples`. The
number of casted rays can be controlled (usually at least 300-500 rays are
required to get a smooth result) and the result is returned in `AO`, as a
single scalar for each sample.

Ambient occlusion can be used to darken the surface colors, as shown in
[Example 606](606_AmbientOcclusion/main.c)

![A mesh rendered without (left) and with (right) ambient
occlusion.](images/606_AmbientOcclusion.png)

## Picking [607]

Picking vertices and faces using the mouse is very common in geometry
processing applications. While this might seem a simple operation, its
implementation is not straighforward. Libigl contains a function that solves this problem using the
[Embree](https://software.intel.com/en-us/articles/embree-photo-realistic-ray-tracing-kernels)
raycaster. Its usage is demonstrated in [Example 607](607_Picking/main.cpp):

```cpp
bool hit = igl::unproject_onto_mesh(
  Vector2f(x,y),
  F,
  viewer.core.view * viewer.core.model,
  viewer.core.proj,
  viewer.core.viewport,
  *ei,
  fid,
  vid);
```

This function casts a ray from the view plane in the view direction. Variables
`x` and `y` are
the mouse screen coordinates; `view`, `model`, `proj` are the view, model and
projection matrix respectively; `viewport` is the viewport in OpenGL format;
`ei`
contains a [Bounding Volume
Hierarchy](http://en.wikipedia.org/wiki/Bounding_volume_hierarchy) constructed
by Embree, and `fid` and `vid` are the picked face and vertex, respectively.

![([Example 607](607_Picking/main.cpp)) Picking via ray casting. The selected
vertices are colored in red.](images/607_Picking.png)

## Locally Injective Maps [608]

Extreme deformations or parametrizations with high-distortion might flip
elements.  This is undesirable in many applications, and it is possible to
avoid it by introducing a non-linear constraints that guarantees that the area
of every element remain positive.

Libigl can be used to compute Locally Injective Maps [#schuller_2013][] using a variety of
deformation energies. A simple deformation of a 2D grid is computed in [Example
608](608_LIM/main.cpp).

![A mesh (left) deformed using Laplacian editing (middle) and with Laplacian
editing plus the anti-flipping constraints (right).](images/608_LIM.png)

## Boolean operations on meshes [609]

Constructive solid geometry (CSG) is a technique to define a complex surface as
the result of a number of set operations on solid regions of space: union,
intersection, set difference, symmetric difference, complement. Typically, CSG
libraries represent the inputs and outputs to these operations _implicitly_:
the solid $A$ is defined as the open set of points $\mathbf{x}$ for which some
function $a(\mathbf{x})$ "returns true". The surface of this shape is the
_closure_ of all points $x$ in $A$.

With this sort of representation, boolean
operations are straightforward. For example, the union of solids $A$ and $B$
is simply

$A \cup B = \{\mathbf{x} \left.\right|
  a(\mathbf{x}) \text{ or } b(\mathbf{x})\},$

the intersection is

$A \cap B = \{\mathbf{x} \left.\right|
  a(\mathbf{x}) \text{ and } b(\mathbf{x})\},$

the difference $A$ _minus_ $B$ is

$A \setminus B = \{\mathbf{x} \left.\right|
  a(\mathbf{x}) \text{ and _not_ } b(\mathbf{x})\},$

and the symmetric difference (XOR) is

$A \triangle B = \{\mathbf{x} \left.\right|
  \text{either } a(\mathbf{x}) \text{ or } b(\mathbf{x}) \text{ but not both }\}.$

Stringing together many of these operations, one can design quite complex
shapes. A typical CSG library might only keep explicit _base-case_
representations of canonical shapes: half-spaces, quadrics, etc.

In libigl, we do currently _not_ have an implicit surface representation.
Instead we expect our users to be working with _explicit_ triangle mesh
_boundary representations_ of solid shapes. CSG operations are much hard to
compute robustly with boundary representations, but are nonetheless useful.

To compute a boolean operation on a triangle mesh with vertices `VA` and
triangles `FA` and another mesh `VB` and `FB`, libigl first computes a unified
mesh with vertices `V` and triangles `F` where all triangle-triangle
intersections have been "resolved". That is, edges and vertices are added
exactly at the intersection lines, so the resulting _non-manifold_ mesh `(V,F)`
has no self-intersections.

Then libigl _peals_ the outer hull [#attene_2014][] off this mesh recursively,
keeping track of the iteration parity and orientation flips for each layer.
For any boolean operation, these two pieces of information determine for each
triangle (1) if it should be included in the output, and (2) if its orientation
should be reversed before added to the output.

Calling libigl's boolean operations is simple. To compute the union of
`(VA,FA)` and `(VB,FB)` into a new mesh `(VC,FC)`, use:

```cpp
igl::mesh_boolean(VA,FA,VB,FB,MESH_BOOLEAN_TYPE_UNION,VC,FC);
```

The following figure shows each boolean operation on two meshes.

![The example [Boolean](609_Boolean/main.cpp) conducts
boolean operations on the _Cheburashka_ (red) and _Knight_ (green). From left
to right: union, intersection, set minus, symmetric difference (XOR),
"resolve". Bottom row reveals inner surfaces, darker color indicates
back-facing triangles.](images/cheburashka-knight-boolean.jpg)

The union, symmetric difference and "resolve" have the same outward
appearance, but differ in their treatment of internal structures. The union has
no internal surfaces: the triangles are not included in the output. The
symmetric difference is the same set of triangles as the "resolve", but
internal surfaces have been reversed in orientation, indicating that the solid
result of the operation. The "resolve" operation is not really a boolean
operation, it is simply the result of resolving all intersections and gluing
together coincident vertices, maintaining original triangle orientations.

Libigl also provides a wrapper `igl::mesh_boolean_cork` to the
[cork](https://github.com/gilbo/cork), which is typically faster, but is not
always robust.

# Miscellaneous [700]

Libigl contains a _wide_ variety of geometry processing tools and functions for
dealing with meshes and the linear algebra related to them: far too many to
discuss in this introductory tutorial. We've pulled out a couple of the
interesting functions in this chapter to highlight.

## Mesh Statistics [701]

Libigl contains various mesh statistics, including face angles, face areas and
the detection of singular vertices, which are vertices with more or less than 6
neighbours in triangulations or 4 in quadrangulations.

The example [Statistics](701_Statistics/main.cpp) computes these quantities and
does a basic statistic analysis that allows to estimate the isometry and
regularity of a mesh:

```bash
Irregular vertices:
136/2400 (5.67%)
Areas (Min/Max)/Avg_Area Sigma:
0.01/5.33 (0.87)
Angles in degrees (Min/Max) Sigma:
17.21/171.79 (15.36)
```

The first row contains the number and percentage of irregular vertices, which
is particularly important for quadrilateral meshes when they are used to define
subdivision surfaces: every singular point will result in a point of the
surface that is only C^1.

The second row reports the area of the minimal element, maximal element and the
standard deviation.  These numbers are normalized by the mean area, so in the
example above 5.33 max area means that the biggest face is 5 times larger than
the average face. An ideal isotropic mesh would have both min and max area
close to 1.

The third row measures the face angles, which should be close to 60 degrees (90
for quads) in a perfectly regular triangulation. For FEM purposes, the closer
the angles are to 60 degrees the more stable will the optimization be. In this
case, it is clear that the mesh is of bad quality and it will probably result
in artifacts if used for solving PDEs.

## Generalized Winding Number [702]

The problem of tetrahedralizing the interior of closed watertight surface mesh
is a difficult, but well-posed problem (see our [Tetgen wrappers][605]).  But
black-box tet-meshers like TetGen will _refuse_ input triangle meshes with
self-intersections, open boundaries, non-manifold edges from multiple connected
components.
The problem is two-fold: self-intersections present contradictory facet
constraints and self-intersections/open-boundaries/non-manifold edges make the
problem of determining inside from outside ill-posed without further
assumptions.

The first problem is _easily_ solved by "resolving" all self-intersections.
That is, meshing intersecting triangles so that intersects occur exactly at
edges and vertices. This is accomplished using `igl::selfintersect`.

TetGen can usually tetrahedralize the convex hull of this "resolved" mesh, and
then the problem becomes determining which of these tets are _inside_ the input
mesh and which are outside. That is, which should be kept and which should be
removed.

The "Generalized Winding Number" is a robust method for determined
inside and outside for troublesome meshes [#jacobson_2013][].  The generalized
winding number with respect to `(V,F)` at some point $\mathbf{p} \in
\mathcal{R}^3$ is defined as scalar function:

 $$w(\mathbf{p}) = \sum\limits_{f_i\in F} \frac{1}{4\pi}\Omega_{f_i}(\mathbf{p})$$

where $\Omega_{f_i}$ is the _solid angle_ subtended by $f_i$ (the ith face in
`F`) at the point $\mathbf{p}$. This solid angle contribution is a simple,
closed-form expression involving `atan2` and some dot-products.

If `(V,F)` _does_ form a closed watertight surface, then $w(\mathbf{p})=1$ if
$\mathbf{p}$ lies inside `(V,F)` and $w(\mathbf{p})=0$ if outside `(V,F)`.  If
`(V,F)` is closed but overlaps itself then $w(\mathbf{p})$ is an integer value
counting how many (signed) times `(V,F)` _wraps_ around $\mathbf{p}$.  Finally,
if `(V,F)` is not closed or not even manifold (but at least consistently
oriented), then $w(\mathbf{p})$ tends smoothly toward 1 as $\mathbf{p}$ is
_more_ inside `(V,F)`, and toward 0 as $\mathbf{p}$ is more outside.
 
![Example [702_WindingNumber](702_WindingNumber/main.cpp) computes the
generalized winding number function for a tetrahedral mesh inside a cat with
holes and self intersections (gold). The silver mesh is surface of the
extracted interior tets, and slices show the winding number function on all
tets in the convex hull: blue (~0), green (~1), yellow
(~2).](images/big-sigcat-winding-number.gif)

# Outlook for continuing development [future]

Libigl is in active development, and we plan to focus on the following features
in the next months:

* A better and more consistent **documentation**, plus extending this tutorial
  to cover more libigl features.

* Implement a **mixed-integer solver** which only uses Eigen to remove the
  dependency on CoMiSo.

* Improve the robustness and performance of the active set QP solver. In
  particular, handle linearly dependent constraints.

* Implement more mesh analysis functions, including structural analysis for
  masonry and _3D-printability_ analysis.

* Increase support for point clouds and general polygonal meshes.

* What would you like to see in libigl? [Contact
  us!](mailto:alecjacobson@gmail.com) or post a [feature
  request](https://github.com/libigl/libigl/issues/new).

We encourage you to contribute to the library and to report problems and bugs.
The best way to contribute new feature or bug fixes is to fork the libigl
repository and to open a [pull
request](https://help.github.com/articles/using-pull-requests) on [our github
repository](https://github.com/libigl/libigl).



[#attene_2014]:["Direct repair of
  self-intersecting
  meshes"](https://www.google.com/search?q=Direct+repair+of+self-intersecting+meshes),
  Marco Attene, 2014.
[#bommes_2009]:[Mixed-integer
quadrangulation](http://www-sop.inria.fr/members/David.Bommes/publications/miq.pdf),
David Bommes, Henrik Zimmer, Leif Kobbelt SIGGRAPH 2009
[#botsch_2004]: Matrio Botsch and Leif Kobbelt. ["An Intuitive Framework for
Real-Time Freeform
Modeling,"](https://www.google.com/search?q=An+Intuitive+Framework+for+Real-Time+Freeform+Modeling)
2004.
[#bouaziz_2012]:[Shape-Up: Shaping Discrete Geometry with
Projections](http://lgg.epfl.ch/publications/2012/shapeup.pdf) Sofien Bouaziz,
Mario Deuss, Yuliy Schwartzburg, Thibaut Weise, Mark Pauly
SGP 2012
[#chao_2010]: Isaac Chao, Ulrich Pinkall, Patrick Sanan, Peter Schröder.
["A Simple Geometric Model for Elastic
Deformations,"](https://www.google.com/search?q=A+Simple+Geometric+Model+for+Elastic+Deformations) 2010.
[#diamanti_2014]:[Designing N-PolyVector Fields with Complex
Polynomials](http://igl.ethz.ch/projects/complex-roots/) Olga Diamanti, Amir
Vaxman, Daniele Panozzo, Olga Sorkine-Hornung, SGP 2014
[#eck_2005]:[Multiresolution Analysis of Arbitrary
Meshes](http://research.microsoft.com/en-us/um/people/hoppe/mra.pdf), Matthias
Eck, Tony DeRose, Tom Duchamp, Hugues Hoppe, Michael Lounsbery, Werner
Stuetzle, SIGGRAPH 2005
[#jacobson_thesis_2013]: Alec Jacobson,
[_Algorithms and Interfaces for Real-Time Deformation of 2D and 3D
Shapes_](https://www.google.com/search?q=Algorithms+and+Interfaces+for+Real-Time+Deformation+of+2D+and+3D+Shapes),
2013.
[#jacobson_2013]: Alec Jacobson, Ladislav Kavan, and Olga Sorkine.
["Robust Inside-Outside Segmentation using Generalized Winding
Numbers,"](https://www.google.com/search?q=Robust+Inside-Outside+Segmentation+using+Generalized+Winding+Numbers) 2013.
[#jacobson_2012]: Alec Jacobson, Ilya Baran, Ladislav Kavan, Jovan Popović, and
Olga Sorkine. ["Fast Automatic Skinning
Transformations,"](https://www.google.com/search?q=Fast+Automatic+Skinning+Transformations) 2012.
[#jacobson_2011]: Alec Jacobson, Ilya Baran, Jovan Popović, and Olga Sorkine.
["Bounded Biharmonic Weights for Real-Time Deformation,"](https://www.google.com/search?q=Bounded+biharmonic+weights+for+real-time+deformation) 2011.
[#jacobson_mixed_2010]: Alec Jacobson, Elif Tosun, Olga Sorkine, and Denis
Zorin. ["Mixed Finite Elements for Variational Surface
Modeling,"](https://www.google.com/search?q=Mixed+Finite+Elements+for+Variational+Surface+Modeling) 2010.
[#kavan_2008]: Ladislav Kavan, Steven Collins, Jiri Zara, and Carol O'Sullivan.
["Geometric Skinning with Approximate Dual Quaternion
Blending,"](https://www.google.com/search?q=Geometric+Skinning+with+Approximate+Dual+Quaternion+Blending) 2008.
[#kazhdan_2012]: Michael Kazhdan, Jake Solomon, Mirela Ben-Chen,
["Can Mean-Curvature Flow Be Made
Non-Singular,"](https://www.google.com/search?q=Can+Mean-Curvature+Flow+Be+Made+Non-Singular) 2012.
[#knoppel_2013]:[Globally Optimal Direction
Fields](http://www.cs.columbia.edu/~keenan/Projects/GloballyOptimalDirectionFields/paper.pdf) Knöppel, Crane, Pinkall, Schröder SIGGRAPH 2013
[#levy_2002]: [Least Squares Conformal Maps, for Automatic Texture Atlas
Generation,](http://www.cs.jhu.edu/~misha/Fall09/Levy02.pdf) Bruno Lévy,
Sylvain Petitjean, Nicolas Ray, Jérome Maillot, SIGGRAPH 2002
[#levy_2008]:[N-Symmetry Direction Field
Design](http://alice.loria.fr/publications/papers/2008/DGF/NSDFD-TOG.pdf),
Nicolas Ray, Bruno Vallet, Wan Chiu Li, Bruno Lévy TOG 2008
[#liu_2008]: [A Local/Global Approach to Mesh
Parameterization](http://cs.harvard.edu/~sjg/papers/arap.pdf) Ligang Liu, Lei
Zhang, Yin Xu, Craig Gotsman, Steven J. Gortler SGP 2008
[#liu_2011]:[General Planar Quadrilateral Mesh Design Using Conjugate Direction
Field](http://research.microsoft.com/en-us/um/people/yangliu/publication/cdf.pdf ) Yang Liu, Weiwei Xu, Jun Wang, Lifeng Zhu, Baining Guo, Falai Chen, Guoping
Wang SIGGRAPH Asia 2011
[#mcadams_2011]: Alexa McAdams, Andrew Selle, Rasmus Tamstorf, Joseph Teran,
Eftychios Sifakis. ["Computing the Singular Value Decomposition of 3x3 matrices
with minimal branching and elementary floating point
operations,"](https://www.google.com/search?q=Computing+the+Singular+Value+Decomposition+of+3x3+matrices+with+minimal+branching+and+elementary+floating+point+operations)
2011.
[#meyer_2003]: Mark Meyer, Mathieu Desbrun, Peter Schröder and Alan H.  Barr,
["Discrete Differential-Geometry Operators for Triangulated
2-Manifolds,"](https://www.google.com/search?q=Discrete+Differential-Geometry+Operators+for+Triangulated+2-Manifolds)
2003.
[#mullen_2008]: [Spectral Conformal
Parameterization](http://www.geometry.caltech.edu/pubs/MTAD08.pdf), Patrick
Mullen, Yiying Tong, Pierre Alliez, Mathieu Desbrun, CGF 2008
[#panozzo_2010]: Daniele Panozzo, Enrico Puppo, Luigi Rocca, ["Efficient
Multi-scale Curvature and Crease
Estimation,"](https://www.google.com/search?q=Efficient+Multi-scale+Curvature+and+Crease+Estimation)
2010.
[#panozzo_2014]:[Frame Fields: Anisotropic and Non-Orthogonal Cross
Fields](http://www.inf.ethz.ch/personal/dpanozzo/papers/frame-fields-2014.pdf),
Daniele Panozzo, Enrico Puppo, Marco Tarini, Olga Sorkine-Hornung, SIGGRAPH,
2014
[#rustamov_2011]: Raid M. Rustamov, ["Multiscale Biharmonic
Kernels"](https://www.google.com/search?q=Multiscale+Biharmonic+Kernels), 2011.
[#schuller_2013]:[Locally Injective Mappings](http://igl.ethz.ch/projects/LIM/)
Christian Schüller, Ladislav Kavan, Daniele Panozzo, Olga Sorkine-Hornung,
SGP 2013
[#sharf_2007]: Andrei Sharf, Thomas Lewiner, Gil Shklarski, Sivan Toledo, and
Daniel Cohen-Or. ["Interactive topology-aware surface
reconstruction,"](https://www.google.com/search?q=Interactive+topology-aware+surface+reconstruction) 2007.
[#sorkine_2004]: Olga Sorkine, Yaron Lipman, Daniel Cohen-Or, Marc Alexa,
Christian Rössl and Hans-Peter Seidel. ["Laplacian Surface
Editing,"](https://www.google.com/search?q=Laplacian+Surface+Editing) 2004.
[#sorkine_2007]: Olga Sorkine and Marc Alexa, ["As-rigid-as-possible Surface
Modeling."](https://www.google.com/search?q=As-rigid-as-possible+Surface+Modeling) 2007.
