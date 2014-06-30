title: libigl Tutorial
author: Daniele Panozzo, Alec Jacobson and others
date: 20 June 2014
css: style.css
html header:   <script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
<link rel="stylesheet" href="http://yandex.st/highlightjs/7.3/styles/default.min.css">
<script src="http://yandex.st/highlightjs/7.3/highlight.min.js"></script>
<script>hljs.initHighlightingOnLoad();</script>

# libigl Tutorial notes
Libigl is an open source C++ library for geometry processing research and development.  Dropping the heavy data structures of tradition geometry libraries, libigl is a simple header-only library of encapsulated functions. This combines the rapid prototyping familiar to Matlab or Python programmers with the performance and versatility of C++.  The tutorial is a self-contained, hands-on introduction to libigl.  Via live coding and interactive examples, we demonstrate how to accomplish various common geometry processing tasks such as computation of differential quantities and operators, real-time deformation, global parametrization, numerical optimization and mesh repair.  Each section of these lecture notes links to a cross-platform example application.

# Table of Contents

* [Chapter 1: Introduction to libigl][100]
    * [101 Mesh representation][101]
    * [102 Plotting surfaces][102]
    * [103 Interaction with keyboard and mouse][103]
    * [104 Scalar field visualization][104]
    * [105 Overlays][105]
    * [106 Picking vertices and faces][106]
        * [libigl design principles][107]
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

* [Chapter 5: Parametrization][500]
    * [501 Harmonic parametrization][501]
    * [502 Least-Square Conformal Maps][502]
    * [503 As-Rigid-As-Possible][503]
    * [504 N-Rotationally symmetric tangent fields][504]
    * [505 Global, seamless integer-grid parametrization][505]
    * [506 Anisotropic remeshing using frame fields][506]
    * [507 N-PolyVector fields][507]
    * [508 Conjugate vector fields][508]
    * [509 Planarization][509]

* [Chapter 6: External libraries][600]
    * [601 State serialization][601]
    * [602 Mixing matlab code][602]
    * [603 Calling igl functions from matlab][603]
    * [604 Triangulation of closed polygons][604]
    * [605 Tetrahedralization of closed surfaces][605]
    * [606 Baking ambient occlusion][606]
    * [607 Locally Injective Maps][607]

* [Chapter 7: Outlook for continuing development][future]

# Chapter 1 [100]

We introduce libIGL with a series of self-contained examples. The purpose of each example is to showcase a feature of libIGL while applying to a practical problem in geometry processing. In this chapter, we will showcase the basic concepts of libigl and introduce a simple mesh viewer that allows to easily visualize surface mesh and its attributes. All the examples are cross-platform and can be compiled on MacOSX, Linux and Windows.

libigl can be downloaded from our [github repository](https://github.com/libigl/libigl) or cloned with git:
``` sh
git clone https://github.com/libigl/libigl.git
```

All examples depends on glfw, glew and anttweakbar. A copy
of the sourcecode of each library is provided together with libigl
and they can be precompiled using:

```sh
    sh compile_dependencies_macosx.sh (MACOSX)
    sh compile_dependencies_linux.sh (LINUX)
```
while precompiled binaries are provided for Visual Studio 2014 64bit.

You can use the CMakeLists.txt in the tutorial folder to build all the examples:

```sh
  cd tutorial
  mkdir build
  cd build
  cmake ../
  make
```

or you can use the CMakeLists.txt inside each example folder to build the examples independently.

For a few examples in Chapter 5, the [CoMiSo solver](http://www.graphics.rwth-aachen.de/software/comiso) has to be downloaded and compiled separately.

## Mesh representation [101]

libIGL uses the [Eigen](http://eigen.tuxfamily.org/) library to encode vector and matrices. We will review in this tutorial many of the basic operations that Eigen supports: If you want to get an idea of what operations are supported you can take a look at the [dense](http://eigen.tuxfamily.org/dox/group__QuickRefPage.html) and [sparse](http://eigen.tuxfamily.org/dox/group__SparseQuickRefPage.html) quick reference guides.

We encode a triangular mesh as a pair of matrices:

```cpp
Eigen::MatrixXd V;
Eigen::MatrixXi F;
```

**V** is a #N by 3 matrix which stores the coordinates of the vertices. Each row stores the coordinate of a vertex, with the x,y,z coordinates in the first, second and third column respectively. The matrix **F** stores the triangle connectivity: each line of **F** denotes a triangle whose 3 vertices are represented as indices pointing to vertex coordinates in **F**.

![A simple mesh made of 2 triangles and 4 vertices.](images/VF.png)

Note that the order of the vertex indices in F determines the orientation of the triangles and it should be consistent for the entire surface. As we will see later, additional properties of the mesh will be similarly stored as matrices. This simple representation has many advantages:

* it is memory efficient and cache friendly
* the use of indices instead of pointers greatly simplifies debuggind
* the data can be trivially read/written on disk

libIGL provides Input/Output functions to read and write common mesh formats.
The reading/writing functions are named read\*.h and write\*.h, respectively.

Reading a mesh from file requires a single igl function call:

```cpp
igl::readOFF("../shared/cube.off", V, F);
```

The functions read the mesh cube.off and fills the provided matrices V and F.
Similarly, to write a mesh to file (in OBJ format):

```cpp
igl::writeOBJ("cube.obj",V,F);
```

See [Example 101](101_FileIO/main.cpp) for the source code of a simple mesh converter from OFF to OBJ format.

## Plotting surfaces [102]

libigl contains an OpenGL viewer that can visualize surface and their properties.

The following code ([Example 102](102_DrawMesh/main.cpp)) is a basic skeleton that will be used over the entire tutorial. It is a standalone application that loads a mesh and visualize it.

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
  viewer.set_mesh(V, F);
  viewer.launch();
}
```

The function set_mesh assigns to the viewer the mesh that we want to plot, and the last line creates an opengl context and starts the draw loop. Additional properties can be plotted on the mesh, and it is also possible to extend the viewer with standard OpenGL code. Please see the documentation in  [Viewer.h](../include/igl/Viewer/Viewer.h) for more details.

![([Example 102](102_DrawMesh/main.cpp)) loads and draws a mesh.](images/102_DrawMesh.png)

## Interaction with keyboard and mouse [103]

Keyboard and mouse events triggers callbacks that can be registered in the viewer. The viewer supports the following callbacks:

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

A keyboard callback can be used to visualize multiple meshes or different stages of an algorithm, as demonstrated in [Example 103](103_Events/main.cpp). The keyboard callback changes the visualized mesh depending on the key pressed:

```cpp
bool key_down(igl::Viewer& viewer, unsigned char key, int modifier)
{
  if (key == '1')
  {
    viewer.clear();
    viewer.set_mesh(V1, F1);
  }
  else if (key == '2')
  {
    viewer.clear();
    viewer.set_mesh(V2, F2);
  }
  return false;
}
```
and it is registered in the viewer as follows:

```cpp
viewer.callback_key_down = &key_down;
```

Note that the mesh is cleared before using set_mesh. This has to be called every time the number of vertices or faces of the plotted mesh changes. Every callback returns a boolean value that tells the viewer if the event has been handled by the plugin, or if the viewer should process it normally. This is useful, for example, to disable the default mouse event handling if you want to control the camera directly in your code.

The viewer can be extended using plugins, which are classes that implements all the viewer's callbacks. See the class Viewer_plugin for more details.

## Scalar field visualization [104]

Colors and normals can be associated to both faces or normals using the set_colors function:

```cpp
viewer.set_colors(C);
```
**C** is a #C by 3 matrix with one RGB color per row, and as many rows as the number of faces **or** the number of vertices. Depending on the size of **C**, the viewer applies the color to faces or vertices.

Colors are commonly used to visualize scalar functions defined on a surface using a transfer functions, that maps a scalar value between 0 and 1 to a color scale. A simple example of a scalar field defined on a surface is the z coordinate of each point. We can extract this information from our mesh by taking the first column of V (which contains the stacked z coordiantes of all the vertices), and map it to colors using the igl::jet function:

```cpp
Eigen::VectorXd x = V.col(2);
igl::jet(x,true,C);
```

The first row extracts the third column from V and the second calls the libigl functions that converts a scalar field to colors. The second parameter of jet normalizes the scalar field to lie between 0 and 1 before applying the color scale.

![([Example 104](104_Colors/main.cpp)) igl::jet converts a scalar field to a color field.](images/104_Colors.png)

## Overlays [105]

In addition to the surface, the viewer supports the visualization of points, lines and text label that can be very helful while developing geometric processing algorithms. These additional informations can be drawn using the following functions:

```cpp
viewer.add_points(P,Eigen::RowVector3d(r,g,b));
```

Draws a point of color r,g,b for each row of P at the coordinates specified in each row of P, which is a #P by 3 matrix.

```cpp
viewer.add_edges(P1,P2,Eigen::RowVector3d(r,g,b);
```

Draws a line for each line of P1 and P2, which connects the point in P1 to the point in P2.

```cpp
viewer.add_label(p,str);
```

Draws a label containing the string str at the position p.

These functions are demonstrate in [Example 105](105_Overlays/main.cpp) where the bounding box of the mesh is plotted using lines and points. The bounding box of a mesh can be found using Eigen:

```cpp
Eigen::Vector3d m = V.colwise().minCoeff();
Eigen::Vector3d M = V.colwise().maxCoeff();
```

![([Example 105](105_Overlays/main.cpp)) The bounding box of a mesh is shown using overlays.](images/105_Overlays.png)

Using matrices to encode the mesh and its attributes allows to write short and efficient code for many operations, avoiding to write for loops.

## Picking [106]

Picking vertices and faces using the mouse is very common in geometry processing applications. While this might seem a simple operation, its implementation is quite involved. libigl contains a function that solves this problem using the [Embree](https://software.intel.com/en-us/articles/embree-photo-realistic-ray-tracing-kernels) raycaster. Its usage is demonstrated in [Example 106](106_Picking/main.cpp):

```cpp
bool hit = igl::unproject_in_mesh(
  Vector2f(x,y),
  F,
  viewer.view * viewer.model,
  viewer.proj,
  viewer.viewport,
  *ei,
  fid,
  vid);
```

This function casts a ray from the view plane in the view direction. x,y are the position of the mouse on screen; view,model,proj are the view, model and projection matrix respectively, viewport is the viewport in opengl format; ei contains a [Bounding Volume Hierarchy](http://en.wikipedia.org/wiki/Bounding_volume_hierarchy) constructed by Embree, and fid and vid are the picked face and vertex respectively.

This function is a good example of the design principles in libigl: the function takes very simple types, mostly matrix or vectors, and can be easily reused for many different tasks.
Not committing to heavy data structures, favors simplicity, ease of use and reusability.

# libigl design choices [107]

To conclude the introduction, we summarize the main design principles in libigl:

* No complex data types. Mostly matrices and vectors. This greatly favors code reusability and forces the authors to expose all the parameters used by the algorithm.  

* Minimal dependencies: we use external libraries only when necessary and we wrap them in a small set of functions.

* Header-only: it is straighforward to use our library since it is only one additional include directory in your project. (if you are worried about compilation speed, it is also possible to build the library as a [static library](../build/))

![([Example 106](106_Picking/main.cpp)) Picking via ray casting. The selected vertices are colored in red.](images/106_Picking.png)

# Chapter 2: Discrete Geometric Quantities and Operators
This chapter illustrates a few discrete quantities that libigl can compute on a
mesh. This also provides an introduction to basic drawing and coloring routines
in our example viewer. Finally, we construct popular discrete differential
geometry operators.

## Normals
Surface normals are a basic quantity necessary for rendering a surface. There
are a variety of ways to compute and store normals on a triangle mesh.

### Per-face
Normals are well defined on each triangle of a mesh as the vector orthogonal to
triangle's plane. These piecewise constant normals produce piecewise-flat
renderings: the surface appears non-smooth and reveals its underlying
discretization.

### Per-vertex
Storing normals at vertices, Phong or Gouraud shading will interpolate shading
inside mesh triangles to produce smooth(er) renderings. Most techniques for
computing per-vertex normals take an average of incident face normals. The
techniques vary with respect to their different weighting schemes. Uniform
weighting is heavily biased by the discretization choice, where as area-based
or angle-based weighting is more forgiving.

The typical half-edge style computation of area-based weights might look
something like this:

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

Without a half-edge data-structure it may seem at first glance that looping
over incident faces---and thus constructing the per-vertex normals---would be
inefficient. However, per-vertex normals may be _throwing_ each face normal to
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

Storing normals per-corner is an efficient an convenient way of supporting both
smooth and sharp (e.g. creases and corners) rendering. This format is common to
OpenGL and the .obj mesh file format. Often such normals are tuned by the mesh
designer, but creases and corners can also be computed automatically. Libigl
implements a simple scheme which computes corner normals as averages of
normals of faces incident on the corresponding vertex which do not deviate by a
specified dihedral angle (e.g. 20°).

![The `Normals` example computes per-face (left), per-vertex (middle) and
per-corner (right) normals](images/fandisk-normals.jpg)

## Gaussian Curvature
Gaussian curvature on a continuous surface is defined as the product of the
principal curvatures:

 $k_G = k_1 k_2.$

As an _intrinsic_ measure, it depends on the metric and
not the surface's embedding.

Intuitively, Gaussian curvature tells how locally spherical or _elliptic_ the
surface is ( $k_G>0$ ), how locally saddle-shaped or _hyperbolic_ the surface
is ( $k_G<0$ ), or how locally cylindrical or _parabolic_ ( $k_G=0$ ) the
surface is.

In the discrete setting, one definition for a ``discrete Gaussian curvature''
on a triangle mesh is via a vertex's _angular deficit_:

 $k_G(v_i) = 2π - \sum\limits_{j\in N(i)}θ_{ij},$

where $N(i)$ are the triangles incident on vertex $i$ and $θ_{ij}$ is the angle
at vertex $i$ in triangle $j$ [][#meyer_2003].

Just like the continuous analog, our discrete Gaussian curvature reveals
elliptic, hyperbolic and parabolic vertices on the domain.

![The `GaussianCurvature` example computes discrete Gaussian curvature and visualizes it in
pseudocolor.](images/bumpy-gaussian-curvature.jpg)

## Curvature Directions
The two principal curvatures $(k_1,k_2)$ at a point on a surface measure how much the
surface bends in different directions. The directions of maximum and minimum
(signed) bending are call principal directions and are always
orthogonal.

Mean curvature is defined simply as the average of principal curvatures:

 $H = \frac{1}{2}(k_1 + k_2).$

One way to extract mean curvature is by examining the Laplace-Beltrami operator
applied to the surface positions. The result is a so-called mean-curvature
normal:

  $-\Delta \mathbf{x} = H \mathbf{n}.$

It is easy to compute this on a discrete triangle mesh in libigl using the cotangent
Laplace-Beltrami operator [][#meyer_2003].

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
H = (HN.rowwise().squaredNorm()).array().sqrt();
```

Combined with the angle defect definition of discrete Gaussian curvature, one
can define principal curvatures and use least squares fitting to find
directions [][#meyer_2003].

Alternatively, a robust method for determining principal curvatures is via
quadric fitting [][#panozzo_2010]. In the neighborhood
around every vertex, a best-fit quadric is found and principal curvature values
and directions are sampled from this quadric. With these in tow, one can
compute mean curvature and Gaussian curvature as sums and products
respectively.

![The `CurvatureDirections` example computes principal curvatures via quadric
fitting and visualizes mean curvature in pseudocolor and principal directions
with a cross field.](images/fertility-principal-curvature.jpg)

This is an example of syntax highlighted code:

```cpp
#include <foo.html>
int main(int argc, char * argv[])
{
  return 0;
}
```

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
 \nabla \sum\limits_{i=1}^n \nabla \phi_i(\mathbf{x})\, f_i =
 \sum\limits_{i=1}^n \nabla \phi_i(\mathbf{x})\, f_i.$

This reveals that the gradient is a linear function of the vector of $f_i$
values. Because $\phi_i$ are linear in each triangle their gradient are
_constant_ in each triangle. Thus our discrete gradient operator can be written
as a matrix multiplication taking vertex values to triangle values:

 $\nabla f \approx \mathbf{G}\,\mathbf{f},$

where $\mathbf{f}$ is $n\times 1$ and $\mathbf{G}$ is an $md\times n$ sparse
matrix. This matrix $\mathbf{G}$ can be derived geometrically, e.g.
[ch. 2][#jacobson_thesis_2013].
Libigl's `gradMat`**Alec: check name** function computes $\mathbf{G}$ for
triangle and tetrahedral meshes:

![The `Gradient` example computes gradients of an input function on a mesh and
visualizes the vector field.](images/cheburashka-gradient.jpg)

## Laplacian

The discrete Laplacian is an essential geometry processing tool. Many
interpretations and flavors of the Laplace and Laplace-Beltrami operator exist.

In open Euclidean space, the _Laplace_ operator is the usual divergence of gradient
(or equivalently the Laplacian of a function is the trace of its Hessian):

 $\Delta f =
 \frac{\partial^2 f}{\partial x^2} +
 \frac{\partial^2 f}{\partial y^2} +
 \frac{\partial^2 f}{\partial z^2}.$

The _Laplace-Beltrami_ operator generalizes this to surfaces.

When considering piecewise-linear functions on a triangle mesh, a discrete Laplacian may
be derived in a variety of ways. The most popular in geometry processing is the
so-called ``cotangent Laplacian'' $\mathbf{L}$, arising simultaneously from FEM, DEC and
applying divergence theorem to vertex one-rings. As a linear operator taking
vertex values to vertex values, the Laplacian $\mathbf{L}$ is a $n\times n$
matrix with elements:

$L_{ij} = \begin{cases}j \in N(i) &\cot \alpha_{ij} + \cot \beta_{ij},\\
j \notin N(i) & 0,\\
i = j & -\sum\limits_{k\neq i} L_{ik},
\end{cases}$

where $N(i)$ are the vertices adjacent to (neighboring) vertex $i$, and
$\alpha_{ij},\beta_{ij}$ are the angles opposite edge ${ij}$.
This oft
produced formula leads to a typical half-edge style implementation for
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

Without a half-edge data-structure it may seem at first glance that looping
over one-rings, and thus constructing the Laplacian would be inefficient.
However, the Laplacian may be built by summing together contributions for each
triangle, much in spirit with its FEM discretization of the Dirichlet energy
(sum of squared gradients):

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
**Alec: cite Ariel reconstruction paper**.

The operator applied to mesh vertex positions amounts to smoothing by _flowing_
the surface along the mean curvature normal direction. This is equivalent to
minimizing surface area.

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

 $\nabla f \approx \mathbf{M}^{-1} \mathbf{L} \mathbf{f}.$

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
since the Laplacian is the divergence of gradient. Naturally, $\mathbf{G}^T$ is
$n \times md$ sparse matrix which takes vector values stored at triangle faces
to scalar divergence values at vertices.

# Chapter 3: Matrices and Linear Algebra
Libigl relies heavily on the Eigen library for dense and sparse linear algebra
routines. Besides geometry processing routines, libigl has a few linear algebra
routines which bootstrap Eigen and make Eigen feel even more like a high-level
algebra library like Matlab.

## Slice
A very familiar and powerful routine in Matlab is array slicing. This allows
reading from or writing to a possibly non-contiguous sub-matrix. Let's consider
the matlab code:

```matlab
B = A(R,C);
```

If `A` is a $m \times n$ matrix and `R` is a $j$-long list of row-indices
(between 1 and $m$) and `C` is a $k$-long list of column-indices, then as a
result `B` will be a $j \times k$ matrix drawing elements from `A` according to
`R` and `C`. In libigl, the same functionality is provided by the `slice`
function:

```cpp
VectorXi R,C;
MatrixXd A,B;
...
igl::slice(A,R,C,B);
```

`A` and `B` could also be sparse matrices.

Similarly, consider the matlab code:

```matlab
A(R,C) = B;
```

Now, the selection is on the left-hand side so the $j \times k$ matrix  `B` is
being _written into_ the submatrix of `A` determined by `R` and `C`. This
functionality is provided in libigl using `slice_into`:

```cpp
igl::slice_into(B,R,C,A);
```

![The example `Slice` shows how to use `igl::slice` to change the colors for triangles
on a mesh.](images/decimated-knight-slice-color.jpg)

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

Similarly, sorting entire rows can be accomplished in matlab using:

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
order.](images/decimated-knight-sort-color.jpg)


### Other Matlab-style functions
Libigl implements a variety of other routines with the same api and
functionality as common matlab functions.

- `igl::any_of` Whether any elements are non-zero (true)
- `igl::cat` Concatenate two matrices (especially useful for dealing with Eigen
  sparse matrices)
- `igl::ceil` Round entries up to nearest integer
- `igl::cumsum` Cumulative sum of matrix elements
- `igl::colon` Act like Matlab's `:`, similar to Eigen's `LinSpaced`
- `igl::cross` Cross product per-row
- `igl::dot` dot product per-row
- `igl::find` Find subscripts of non-zero entries
- `igl::floot` Round entries down to nearest integer
- `igl::histc` Counting occurrences for building a histogram
- `igl::hsv_to_rgb` Convert HSV colors to RGB (cf. Matlab's `hsv2rgb`)
- `igl::intersect` Set intersection of matrix elements.
- `igl::jet` Quantized colors along the rainbow.
- `igl::kronecker_product` Compare to Matlab's `kronprod`
- `igl::median` Compute the median per column
- `igl::mode` Compute the mode per column
- `igl::orth` Orthogonalization of a basis
- `igl::setdiff` Set difference of matrix elements
- `igl::speye` Identity as sparse matrix

## Laplace Equation
A common linear system in geometry processing is the Laplace equation:

 $∆z = 0$

subject to some boundary conditions, for example Dirichlet boundary conditions
(fixed value):

 $\left.z\right|_{\partial{S}} = z_{bc}$

In the discrete setting, this begins with the linear system:

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
 \mathbf{L}_{b}\end{array}\right) =
 \left(\begin{array}{c}
 \mathbf{0}_{in}\\
 \mathbf{*}_{b}\end{array}\right)$$

The bottom block of equations is no longer meaningful so we'll only consider
the top block:

 $$\left(\begin{array}{cc}
 \mathbf{L}_{in,in} & \mathbf{L}_{in,b}\end{array}\right)
 \left(\begin{array}{c}
 \mathbf{z}_{in}\\
 \mathbf{z}_{b}\end{array}\right) =
 \mathbf{0}_{in}$$

Where now we can move known values to the right-hand side:

 $$\mathbf{L}_{in,in}
 \mathbf{z}_{in} = -
 \mathbf{L}_{in,b}
 \mathbf{z}_{b}$$

Finally we can solve this equation for the unknown values at interior vertices
$\mathbf{z}_{in}$.

However, probably our vertices are not sorted. One option would be to sort `V`,
then proceed as above and then _unsort_ the solution `Z` to match `V`. However,
this solution is not very general.

With array slicing no explicit sort is needed. Instead we can _slice-out_
submatrix blocks ($\mathbf{L}_{in,in}$, $\mathbf{L}_{in,b}$, etc.) and follow
the linear algebra above directly. Then we can slice the solution _into_ the
rows of `Z` corresponding to the interior vertices.

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

Let's consider a general quadratic minimization problem subject to different
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
constraint right-hand sides:

```cpp
igl::min_quad_with_fixed_solve(mqwf,B,bc,Beq,Z);
```

The output `Z` is a $n \times 1$ vector of solutions with fixed values
correctly placed to match the mesh vertices `V`.

## Linear Equality Constraints
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
$1$
in the column corresponding to vertex $c$ and a $-1$ at $d$. The right-hand side
`Beq` is simply zero.

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
difference from
the straight quadratic
_minimization_ system, is that
this saddle problem system will not be positive definite. Thus, we must use a
different factorization technique (LDLT rather than LLT). Luckily, libigl's
`min_quad_with_fixed_precompute` automatically chooses the correct solver in
the presence of linear equality constraints.

![The example `LinearEqualityConstraints` first solves with just fixed value
constraints (left: 1 and -1 on the left hand and foot respectively), then
solves with an additional linear equality constraint (right: points on right
hand and foot constrained to be equal).](images/cheburashka-biharmonic-leq.jpg)

## Quadratic Programming

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

After deciding which constraints are active each iteration simple reduces to a
quadratic minimization subject to linear _equality_ constraints, and the method
from the previous section is invoked. This is repeated until convergence.

Currently the implementation is efficient for box constraints and sparse
non-overlapping linear inequality constraints.

Unlike alternative interior-point methods, the active set method benefits from
a warm-start (initial guess for the solution vector $\mathbf{z}$).

```cpp
igl::active_set_params as;
// Z is optional initial guess and output
igl::active_set(Q,B,b,bc,Aeq,Beq,Aieq,Bieq,lx,ux,as,Z);
```

![The example `QuadraticProgramming` uses an active set solver to optimize
discrete biharmonic kernels at multiple scales [#rustamov_2011][].](images/cheburashka-multiscale-biharmonic-kernels.jpg)

# Chapter 4: Shape Deformation
Modern mesh-based shape deformation methods satisfy user deformation
constraints at handles (selected vertices or regions on the mesh) and propagate
these handle deformations to the rest of shape _smoothly_ and _without removing
or distorting details_. Libigl provides implementations of a variety of
state-of-the-art deformation techniques, ranging from quadratic mesh-based
energy minimizers, to skinning methods, to non-linear elasticity-inspired
techniques.

## Biharmonic Deformation
The period of research between 2000 and 2010 produced a collection of
techniques that cast the problem of handle-based shape deformation as a
quadratic energy minimization problem or equivalently the solution to a linear
partial differential equation.

There are many flavors of these techniques, but
a prototypical subset are those that consider solutions to the bi-Laplace
equation, that is biharmonic functions [#botsch_2004][]. This fourth-order PDE provides
sufficient flexibility in boundary conditions to ensure $C^1$ continuity at
handle constraints (in the limit under refinement) [#jacobson_mixed_2010][].

### Biharmonic surfaces
Let us first begin our discussion of biharmonic _deformation_, by considering
biharmonic _surfaces_. We will casually define biharmonic surfaces as surface
whose _position functions_ are biharmonic with respect to some initial
parameterization:

 $\Delta \mathbf{x}' = 0$

and subject to some handle constraints, conceptualized as "boundary
conditions":

 $\mathbf{x}'_{b} = \mathbf{x}_{bc}.$

where $\mathbf{x}'$ is the unknown 3D position of a point on the surface. So we are
asking that the bi-Laplace of each of spatial coordinate functions to be zero.

In libigl, one can solve a biharmonic problem like this with `igl::harmonic`
and setting $k=2$ (_bi_-harmonic):

```cpp
// U_bc contains deformation of boundary vertices b
igl::harmonic(V,F,b,U_bc,2,U);
```

This produces smooth surfaces that interpolate the handle constraints, but all
original details on the surface will be _smoothed away_. Most obviously, if the
original surface is not already biharmonic, then giving all handles the identity
deformation (keeping them at their rest positions) will **not** reproduce the
original surface. Rather, the result will be the biharmonic surface that does
interpolate those handle positions.

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

 $\Delta \mathbf{d} = 0$

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

![The `BiharmonicDeformation` example deforms a statue's head as a _biharmonic
surface_ (top) and using a _biharmonic displacements_ (bottom).](images/max-biharmonic.jpg)

#### Relationship to "differential coordinates" and Laplacian surface editing
Biharmonic functions (whether positions or displacements) are solutions to the
bi-Laplace equation, but also minimizers of the "Laplacian energy". For
example, for displacements $\mathbf{d}$, the energy reads

 $\int\limits_S \|\Delta \mathbf{d}\|^2 dA.$

By linearity of the Laplace(-Beltrami) operator we can reexpress this energy in
terms of the original positions $\mathbf{x}$ and the unknown positions
$\mathbf{x}' = \mathbf{x} - \mathbf{d}$:

 $\int\limits_S \|\Delta (\mathbf{x}' - \mathbf{x})\|^2 dA = \int\limits_S \|\Delta \mathbf{x}' - \Delta \mathbf{x})\|^2 dA.$

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

![The `PolyharmonicDeformation` example deforms a flat domain (left) into a bump as a
solution to various $k$-harmonic PDEs.](images/bump-k-harmonic.jpg)

## Bounded Biharmonic Weights
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
are sent as vertex shader _attribtues_ and bone transformations are sent as
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
smoothness by minimizing a smoothness energy: the familiar Laplacian energy:

 $\sum\limits_{i = 1}^m \int_S (\Delta w_i)^2 dA$

subject to constraints which enforce interpolation of handle constraints:

 $w_i(\mathbf{x}) = \begin{cases} 1 & \text{ if } \mathbf{x} \in H_i\\ 0 & \text{ otherwise }
 \end{cases},$

where $H_i$ is the ith handle, and constraints which enforce non-negativity,
parition of unity and encourage sparsity:

 $0\le w_i \le 1$ and $\sum\limits_{i=1}^m w_i = 1.$

This is a quadratic programming problem and libigl solves it using its active
set solver or by calling out to Mosek.

![The example `BoundedBiharmonicWeights` computes weights for a tetrahedral
mesh given a skeleton (top) and then animates a linear blend skinning
deformation (bottom).](images/hand-bbw.jpg)

## Dual Quaternion Skinning
Even with high quality weights, linear blend skinning is limited. In
particular, it suffers from known artifacts stemming from blending rotations as
as matrices: a weight combination of rotation matrices is not necessarily a
rotation. Consider an equal blend between rotating by $-pi/2$ and by $pi/2$
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
transformation of bone $i$. The normalization forces the result of the linear blending
to again be a unit dual quaternion and thus also a rigid transformation.

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


# Chapter 5: Parametrization [500]

In computer graphics, we denote as parametrization of a surface a map from the surface to \\(\mathbf{R}^2\\). It is usually encoded by a new set of 2D coordinates for each vertex of the mesh and possibly also a new set of faces in one to one correspondence with the faces of the original surface. Note that this definition
is the *inverse* of the classical differential geometry parametrization.

A parametrization has many applications, ranging from texture mapping to surface remeshing. Many algorithms have been proposed, and they can be broadly characterized in four families:

1. **Single patch, fixed boundary**: these algorithm can parametrize a disk-like part of the surface given fixed 2D positions for its boundary. These algorithms are efficient and simple, but usually produce high-distortion maps due to the strong cosntraints on the border.

2. **Single patch, free border:** these algorithms allows the boundary to deform freely, reducing the distortion of the map. Care should be taken to prevent the border to self-intersect.

3. **Global parametrization**: these algorithms works on meshes with arbitrary genus. They cut the mesh in multiple patches that are then possible to flatten in a 2D domain. Usually the map is discontinuous on the seams.

4. **Global seamless parametrization**: similar to the global parametrization, but solved with a global solving strategy that hides the seams, making the parametrization "continuous", under specific assumptions that we will discuss later.

## Harmonic parametrization [501]

Harmonic parametrization is a single patch, fixed boundary parametrization algorithm that computes the 2D coordinates of the flattened mesh as two Harmonic functions.

The algorithm is divided in 3 steps:

* Detection of the boundary vertices

```cpp
Eigen::VectorXi bnd;
igl::boundary_loop(V,F,bnd);
```

* Map the boundary vertices to a circle

```cpp
Eigen::MatrixXd bnd_uv;
igl::map_vertices_to_circle(V,bnd,bnd_uv);
```

* Computation of harmonic functions for both the u and v coordinate on the plane, using the boundary positions as boundary constraints

```cpp
igl::harmonic(V,F,bnd,bnd_uv,1,V_uv);
```

bnd contains the indices of the boundary vertices, bnd_uv their position on the UV plane, and "1" denotes that we want to compute an harmonic function (2 will be for biharmonic, 3 for triharmonic, etc.). Note that each of the three functions is deisgned to be reusable in other parametrization algorithms.

A UV parametrization can be visualized in the viewer using the method:

```cpp
viewer.set_uv(V_uv);
```

which uses the UV coordinates to apply a procedural checkerboard texture to the mesh ([Example 501](501_HarmonicParam/main.cpp)).

![([Example 501](501_HarmonicParam/main.cpp)) Harmonic parametrization. (left) mesh with texture, (right) UV parametrization with texture](images/501_HarmonicParam.png)

###References:

[Multiresolution Analysis of Arbitrary Meshes](http://research.microsoft.com/en-us/um/people/hoppe/mra.pdf),
Matthias Eck, Tony DeRose, Tom Duchamp, Hugues Hoppe, Michael Lounsbery, Werner Stuetzle,
SIGGRAPH 2005

## Least-Square Conformal Maps [502]

Least-square conformal maps parametrization minimizes the conformal (angular) distortion of the generated parametrization. If does not need to have a fixed boundary.

LSCM minimizes the following energy:

\\[ E_{LSCM}(\mathbf{u},\mathbf{v}) = \int_X \frac{1}{2}| \nabla \mathbf{u}^{\perp} - \nabla \mathbf{v} |^2 dA \\]

which can be rewritten in matrix form as:

\\[ E_{LSCM}(\mathbf{u},\mathbf{v}) = \frac{1}{2} [\mathbf{u},\mathbf{v}]^t (L_c - 2A) [\mathbf{u},\mathbf{v}] \\]

where L_c is the cotangent laplacian matrix and A is a matrix such that \\( [\mathbf{u},\mathbf{v}]^t A  [\mathbf{u},\mathbf{v}] \\) is equal to the _vector area_ of the mesh.

Using libigl, this matrix energy can be written using a few lines of codes. The cotangent matrix can be computed using igl::cotmatrix:

```cpp
SparseMatrix<double> L;
igl::cotmatrix(V,F,L);
```

Note that we want to apply the laplacian matrix to the u and v coordinates at the same time, thus we need to extend the laplacian matrix taking the left Kronecker product with a 2x2 identity matrix:

```cpp
SparseMatrix<double> L_flat;
repdiag(L,2,L_flat);
```

The area matrix is computed with igl::vector_area_matrix:

```cpp
SparseMatrix<double> A;
igl::vector_area_matrix(F,A);
```

The final energy matrix is the sum of these two matrices. Note that in this case we don't need to fix the boundary, we only need to fix two arbitrary vertices to arbitrary positions to remove the null space of the energy and make the minimum unique. The full source code is provided in [Example 502](502_LSCMParam/main.cpp).


![([Example 502](502_LSCMParam/main.cpp)) LSCM parametrization. (left) mesh with texture, (right) UV parametrization with texture](images/502_LSCMParam.png)

####References:

[Least Squares Conformal Maps, for Automatic Texture Atlas Generation,](http://www.cs.jhu.edu/~misha/Fall09/Levy02.pdf)
Bruno Lévy, Sylvain Petitjean, Nicolas Ray, Jérome Maillot,
SIGGRAPH 2002

[Spectral Conformal Parameterization](http://www.geometry.caltech.edu/pubs/MTAD08.pdf),
Patrick Mullen, Yiying Tong, Pierre Alliez, Mathieu Desbrun,
CGF 2008

## As-Rigid-As-Possible parametrization [503]

As-Rigid-As-Possible parametrizationis a powerful single-patch, non-linear algorithm to compute a parametrization that strives to preserve distances (and thus angles). The idea is very similar to ARAP surface deformation: each triangle is mapped to the plane trying to preserve its original shape, up to a rigid 2x2 rotation.

The algorithm can be implemented reusing the functions discuss in the deformation chapter arap_precomputation and ara_solve. The only difference is that the optimization has to be done in 2D instead of 3D and that a starting point for the non-linear optimization is necessary. While for 3D deformation the original mesh is a perfect starting point, this is not the case for ARAP parametrization since the starting point must be a 2D mesh. In [Example 503](503_ARAPParam/main.cpp), we use Harmonic parametrization as a starting point for the ARAP parametrization: note that similarly to LSCM, the boundary is free to deform to minimize the distortion.

![([Example 503](502_ARAPParam/main.cpp)) As-Rigid-As-Possible parametrization. (left) mesh with texture, (right) UV parametrization with texture](images/503_ARAPParam.png)

### References
[A Local/Global Approach to Mesh Parameterization](http://cs.harvard.edu/~sjg/papers/arap.pdf)
Ligang Liu, Lei Zhang, Yin Xu, Craig Gotsman, Steven J. Gortler
SGP 2008

## N-Rotationally symmetric tangent fields [504]

The design of tangent fields is a basic tool used to design guidance fields for uniform quadrilateral and hexaedral remeshing. libigl contains an implementation of all the state- of-the-art to design algorithms for N-RoSy fields and their generalizations.

In libigl, tangent unit-length vector fields are piece-wise constant on the faces of a triangle mesh, and described by one or more vectors per-face. The function

```cpp
igl::nrosy(V,F,b,bc,b_soft,b_soft_weight,bc_soft,N,0.5,
           output_field,output_singularities);
```

creates a smooth vector field (N=1) starting from a sparse set of constrained faces, whose indices are listed in b and their constrained value is specified in bc. The functions supports soft_constraints (b_soft,b_soft_weight,bc_soft), and returns the interpolated field for each face of the triangle mesh (output_field) plus the singularities of the field (output_singularities).

![Design of a unit-lenght vector field](images/504_vector_field.png)

The singularities are vertices where the field vanishes, and they are highlighted in red. igl::nrosy can generate N-Rotation Symmetric fields, which are a generalization of vector fields where in every face the vector is defined up to a constant rotation of \\( 2\pi / N \\). As can be observed in the following figure, the singularities of fields generated with different N are in different positions and of a different kind.

![Design of a 2-,4- and 9-RoSy field](images/504_nrosy_field.png)

We demonstrate how to call and plot N-RoSy fields in [Example 504](504_NRosyDesign/main.cpp), where the degree of the field can be controlled by pressing the number keys.

### References

[N-Symmetry Direction Field Design](http://alice.loria.fr/publications/papers/2008/DGF/NSDFD-TOG.pdf),
Nicolas Ray, Bruno Vallet, Wan Chiu Li, Bruno Lévy
TOG 2008

[Mixed-integer quadrangulation](http://www-sop.inria.fr/members/David.Bommes/publications/miq.pdf),
David Bommes, Henrik Zimmer, Leif Kobbelt
SIGGRAPH 2009

#Global, seamless integer grid parametrization

The previous parametrization methods where focusing on generating parametrization of single patches, mainly aimed at texture mapping and baking of other surface properties like normals high-frequency details. Global, seamless parametrization aims at parametrizing complex shapes with a parametrization that is aligned with a given set of directions for the purpose of remeshing the surface. In libigl, we provide a reference  implementation of the pipeline of the  [MIQ](http://www-sop.inria.fr/members/David.Bommes/publications/miq.pdf) paper.

### Global, seamless integer-grid parametrization [505]

The first step involves the design of a 4-RoSy field (sometimes called cross field) that describes how the edges of the final quad remeshing should align. The field constraints are usually manually specified or extracted from curvature. In this example, we simply fix one face in a random direction.

![Initial cross field prescribing the edge alignment.](images/505_MIQ_1.png)

### Combing and cutting

Given the cross field, we now want to cut the surface so that it becomes homeorphic to a disk. While this can be done directly on the cross-field, we prefer to do this operation on its bisector field (a copy of the field rotated by 45 degrees) since it is more stable and generic.

We thus rotate the field,

![Bisector field.](images/505_MIQ_2.png)

and we remove the rotation ambiguity by assigning to each face a u and v direction, computed by diffusing this alignment from a random face.

![Combed bisector field.](images/505_MIQ_3.png)

You can imagine this process as combing an hairy surface: you'll be able to comb part of it, but at some point you will not be able to comb consistently the full surface ([Hairy ball theorem](http://en.wikipedia.org/wiki/Hairy_ball_theorem)). The discontinuites in the combing defines the cut graph:

![Cut graph.](images/505_MIQ_4.png)

Finally, we rotate the combed field by 45 degrees to undo the initial 45 degrees rotation:

![Combed cross field.](images/505_MIQ_5.png)

This cross field can be seen as the ideal gradient of the parametrization that we want to compute.

### Poisson parametrization

The mesh can be then cut along the seams and a parametrization is computed trying to find two scalar functions whose gradient matches the combed cross field. This is a classical Poisson problem, that is solved minimizing the following quadratic energy:

\\[ E(\mathbf{u},\mathbf{v}) = |\nabla \mathbf{u} - X_u|^2 + |\nabla \mathbf{v} - X_v|^2 \\]

where \\( X_u \\) and \\( X_u \\) denotes the combed cross field. Solving this problem generates a parametrization whose u and v isolines are aligned with the input cross field.

![Poisson parametrization.](images/505_MIQ_8.png)

We hide the seams by adding a set of integer constraints to the Poisson problem that aligns the isolines on both sides of each seam.

![Seamless Poisson parametrization.](images/505_MIQ_7.png)

Note that this parametrization can only be used for remeshing purposes, since it contains many overlaps.

![Seamless Poisson parametrization (in 2D).](images/505_MIQ_6.png)

A quad mesh can be extracted from this parametrization using
[libQEx](https://github.com/hcebke/libQEx) (not included in libigl).

The full pipeline is demonstrated in [Example 505](505_MIQ/main.cpp).

### References

[Mixed-integer quadrangulation](http://www-sop.inria.fr/members/David.Bommes/publications/miq.pdf),
David Bommes, Henrik Zimmer, Leif Kobbelt
SIGGRAPH 2009

## Anisotropic remeshing [506]

Anisotropic and non-uniform quad remeshing is important to concentrate the elements in the regions with more details. It is possible to extend the MIQ quad meshing framework to generate anisotropic quad meshes using a mesh deformation approach.

The input of the remeshing algorithm is now a sparse set of constraints that defines the shape and scale of the desired quad remeshing. This can be encoded as a frame-field, which is a pair of non-orthogonal and non-unit lenght vectors. The frame field can be interpolated by decomposing it in a 4-RoSy field and a unique affine transformation. The two parts can then be interpolated separately, using igl::nrosy for the cross field, and an harmonic interpolant for the affine part.

![Interpolation of a frame field. Colors on the vectors denote the desired scale. The red faces contains the frame field constraints.](images/506_FrameField_1.png)

After the interpolation, the surface is warped to transform each frame into an orthogonal and unit lenght cross (i.e. removing the scaling and skewness from the frame). This deformation defines a new embedding (and a new metric) for the surface.

![The surface is deformed to transform the frame field in a cross field.](images/506_FrameField_2.png)

The deformed surface can the be isotropically remeshed using the MIQ algorithm that has been presented in the previous section.

![The deformed surface is isotropically remeshed.](images/506_FrameField_3.png)

The UV coordinates of the deformed surface can then be used to transport the parametrization to the original surface, where the isolines will trace a quad mesh whose elements are similar to the shape prescribed in the input frame field.

![The global parametrization is lifted to the original surface to create the anisotropic quad meshing.](images/506_FrameField_4.png)

Our implementation ([Example 506](506_FrameField/main.cpp)) uses MIQ to generate the UV parametrization, but other algorithms could be applied: the only desiderata is that the generated quad mesh will be as isotropic as possible.

### References

[Frame Fields: Anisotropic and Non-Orthogonal Cross Fields](http://www.inf.ethz.ch/personal/dpanozzo/papers/frame-fields-2014.pdf),
Daniele Panozzo, Enrico Puppo, Marco Tarini, Olga Sorkine-Hornung,
SIGGRAPH, 2014

## N-PolyVector fields [507]

N-RoSy vector fields can be further generalized to represent arbitrary vector-sets, with arbitrary angles between them and with arbitrary lenghts. This generalization is called  N-PolyVector field, and libigl provides the function igl::n_polyvector to design them starting from a sparse set of constraints ([Example 507](507_PolyVectorField/main.cpp)).

![Interpolation of a N-PolyVector field from a sparse set of constraints.](images/507_PolyVectorField.png)

Globally Optimal Direction fields are a special case of Poly-Vector Fields: if the interpolation constraints are an N-RoSy field, then our algorithm generates a field that if normalized, is equivalent to a globally optimal direction field.

### References

[Designing N-PolyVector Fields with Complex Polynomials](http://igl.ethz.ch/projects/complex-roots/)
Olga Diamanti, Amir Vaxman, Daniele Panozzo, Olga Sorkine-Hornung,
SGP 2014

[Globally Optimal Direction Fields](http://www.cs.columbia.edu/~keenan/Projects/GloballyOptimalDirectionFields/paper.pdf)
Knöppel, Crane, Pinkall, Schröder
SIGGRAPH 2013

## Conjugate vector fields [508]

Two tangent vectors lying on a face of triangle mesh are conjugate if

\\[ k_1 (u^T d_1)(v^T d_1) + k_2(u^T d_2)(v^T d_2) = 0. \\]

This condition is very important in architectural geometry, since the faces of an infinitely dense quad mesh whose edges are aligned with a conjugate field are planar. Thus, creating a quad mesh following a Conjugate field generates quad meshes that are easier to planarize.

Finding a conjugate vector field that satisfies given directional constraints is a standard problem in architectural geometry, which can be tackled by warping a given PolyVector field to the closest conjugate field.

The algorithms alternates a global step, which enforces smoothness, with a local step that projects the field on every face to the closest conjugate field ([Example 508](508_ConjugateField/main.cpp)).

![A smooth 4-PolyVector field (left) is deformed to become a conjugate field (right).](images/508_ConjugateField.png)

### References

[Designing N-PolyVector Fields with Complex Polynomials](http://igl.ethz.ch/projects/complex-roots/)
Olga Diamanti, Amir Vaxman, Daniele Panozzo, Olga Sorkine-Hornung,
SGP 2014

[General Planar Quadrilateral Mesh Design Using Conjugate Direction Field](http://research.microsoft.com/en-us/um/people/yangliu/publication/cdf.pdf)
Yang Liu, Weiwei Xu, Jun Wang, Lifeng Zhu, Baining Guo, Falai Chen, Guoping Wang
SIGGRAPH Asia 2011

## Planarization [509]

A quad mesh can be transformed in a planar quad mesh using Shape-Up, that is a local/global approach
the uses the global step to enforce surface continuity and the local step to enforce planarity.

[Example 509](509_Planarization/main.cpp) planarizes a quad mesh until it satisfies a user-given planarity threshold.

![A non-planar quad mesh (left) is planarized using the libigl function igl::palanarize (right). The colors represent the planarity of the quads.](images/509_Planarization.png)

### References

[Shape-Up: Shaping Discrete Geometry with Projections](http://lgg.epfl.ch/publications/2012/shapeup.pdf)
Sofien Bouaziz, Mario Deuss, Yuliy Schwartzburg, Thibaut Weise, Mark Pauly
SGP 2012

# Chapter 6: External libraries [600]

An additional positive side effect of using matrices as basic types is that it is easy to exchange data between libigl and other softwares and libraries.

## State serialization [601]

Geometry processing applications often require a considerable amount of computational time and/or manual input. In order to make the development efficient it must be possible to serialize and deserialize the state of the application.

Having a good serialization framework allows to quickly start debugging just before the crash happens, avoiding to wait for the precomputation to take place every time. It also makes it easier to define unit testing that can be used to find bugs in interactive applications: if the input is slightly different every time the algorithm is executed, it is very difficult to find bugs.

Unfortunately, serialization is often not considered in geoemtry processing due to the extreme difficulty in serializing pointer-based data structures (like an helf-edge).

In libigl, serialization is simpler, since the majority of the functions use basic types, and pointers are used in very rare cases (usually to interface with external libraries). libigl provides an extremely easy to use XML serialization framework, that drastically reduces the overhead required to add serialization to your applications.

Assume that the state of your application is composed of a mesh and set of integer ids:

``` cpp
class State : public ::igl::XMLSerialization
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

A class can be made serializable by inheriting from ::igl::XMLSerialization and trivially implementing the InitSerialization method. Note that you don't have to care the types, Add is able to serialize all basic stl types, all Eigen types and any class inheriting from ::igl::XMLSerialization.

It is then possible to save the state to an xml file:

``` cpp
::igl::XMLSerializer serializer_save("601_Serialization");
serializer_save.Add(state,"State");
serializer_save.Save("temp.xml",true);
```

This code generates the following xml file (assuming V and F contains a simple mesh with two triangles, and ids contains the numbers 6 and 7):

``` xml
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

The xml file can then be loaded in a similar way:

``` cpp
State loaded_state;
::igl::XMLSerializer serializer_load("601_Serialization");
serializer_load.Add(loaded_state,"State");
serializer_load.Load("temp.xml");
```

This can also be used as a convenient interface to provide parameters to command line applications, since the xml files can be directly edited with a standard text editor.

We demonstrate the serialization framework in [Example 601](601_Serialization/main.cpp). We strongly suggest that you make the entire state of your application always serializable: this will save you a lot of troubles when you'll be making figures for a scientific publication. It is very common to have to do small changes to figures during the production of a paper, and being able to serialize the entire state just before you take screenshots will save you many painful hours before a submission deadline.

## Mixing matlab code [602]

libigl can be interfaced matlab, to offload some of the numerically heavy computation to a matlab script. This has the major advantage of allowing to develop efficient and complex UI in C++, while keeping the advantage of fast protototyping of matlab. In particular, using an external matlab script in a libigl application allows to change the algorithm in the matlab script without having to recompile the C++ part.

We demonstrate how to integrate matlab in a libigl application in [Example 602](602_Matlab/main.cpp). The example uses matlab to compute the Eigenfunctions of the discrete Laplacian operator, relying on libigl for mesh IO, visualization and for computing the Laplacian operator.

libigl can connect to an existing instance of matlab (or launching a new one on Linux/MacOSX) using:

``` cpp
igl::mlinit(&engine);
```

The cotangent laplacian is computed using igl::cotmatrix and uploaded to the matlab workspace:

``` cpp
igl::cotmatrix(V,F,L);
igl::mlsetmatrix(&engine,"L",L);
```

It is now possible to use any matlab function on the data. For example, we can see the sparsity pattern of L using spy:

``` cpp
igl::mleval(&engine,"spy(L)");
```

![The matlab spy function is called from a libigl-based application.](images/602_Matlab_1.png)

You can also do some computation and then return it back to the C++ application

``` cpp
igl::mleval(&engine,"[EV,~] = eigs(-L,10,'sm')");
igl::mlgetmatrix(&engine,"EV",EV);
```

and then use libigl functions to plot the eigenfunctions.

![4 Eigenfunctions of the Laplacian plotted in the libigl viewer.](images/602_Matlab_2.png)

## Calling igl functions from matlab [603]

It is also possible to call libigl functions from matlab, compiling them as MEX functions. This can be very useful to offload to C++ code the computationally intensive parts of a matlab application.

We provide a wrapper for igl::readOBJ in [Example 603](603_MEX/compileMEX.m). We plan to provide wrappers for all our functions in the future, if you are interested in this feature (or if you want to help implementing it) please let us know.

## Triangulation of closed polygons [604]

The generation of high-quality triangle and tetrahedral meshes is a very common task in geometry processing. We provide wrappers in libigl to triangle and tetegen.

A triangle mesh canb e cerated starting from a set of boundary edges using igl::triangulate.

``` cpp
igl::triangulate(V,E,H,V2,F2,"a0.005q");
```

where E is a set of boundary edges, H a set of 2D positions of points contained in holes of the triangulation and (V2,F2) is the generate triangulation. Additional parameters can be given to triangles, to control the quality: "a0.005q" puts a bound on the maximal area of the triangles and a minimal angle of 20 degrees. In Example [Example 604](604_Triangle/main.m), the interior of a square (excluded a smaller square in its interior) is triangulated.

![Triangulation of the interior of a polygon.](images/604_Triangle.png)

## Tetrahedralization of closed surfaces [605]

Similarly, the interior of a closed manifold surface can be tetrahedralized using the function igl::tetrahedralize which wraps the tetgen library ([Example 605](605_Tetgen/main.c)):

``` cpp
igl::tetrahedralize(V,F,"pq1.414", TV,TT,TF);
```

![Tetrahedralization of the interior of a surface mesh.](images/605_Tetgen.png)

## Baking ambient occlusion [606]

[Ambient occlusion](http://en.wikipedia.org/wiki/Ambient_occlusion) is a rendering technique used to calculate the exposure of each point in a surface to ambient lighting. It is usually encoded as a scalar (normalized between 0 and 1) associated with the vertice of a mesh.

Formally, ambient occlusion is defined as:

\\[ A_p = \frac{1}{\pi} \int_\omega V_{p,\omega}(n \cdot \omega) d\omega \\]

where \\( V_{p,\omega} \\) is the visibility function at  p, defined to be zero if p is occluded in the direction \\( \omega \\) and one otherwise, and \\( d\omega \\) is the infinitesimal solid angle step of the integration variable \\( \omega \\).

The integral is usually approximated by casting rays in random directions around each vertex. This approximation can be computed using the function:

``` cpp
igl::ambient_occlusion(V,F,V_samples,N_samples,500,AO);
```

that given a scene described in V,F, computes the ambient occlusion of the points in V_samples whose associated normals are N_samples. The number of casted rays can be controlled (usually at least 400-500 rays are required to get a smooth result) and the result is return in AO, as a single scalar for each sample.

Ambient occlusion can be used to darken the surface colors, as shown in [Example 606](606_AmbientOcclusion/main.c)

![A mesh rendered without (left) and with (right) ambient occlusion.](images/606_AmbientOcclusion.png)

## Locally Injective Maps [607]

Extreme deformations or parametrizations with high-distortion might flip elements.
This is undesirable in many applications, and it is possible to avoid it by introducing a non-linear contraints that guarantees that the area of every element remain positive.

libigl can be used to compute Locally Injective Maps using a variety of deformation energies. A simple deformation of a 2D grid is computed in [Example 607](607_LIM/main.c).

![A mesh (left) deformed using Laplacian editing (middle) and with Laplacian editing plus the anti-flipping conatraints (right).](images/607_LIM.png)

### References

[Locally Injective Mappings](http://igl.ethz.ch/projects/LIM/)
Christian Schüller, Ladislav Kavan, Daniele Panozzo, Olga Sorkine-Hornung,
SGP 2013

# Outlook for continuing development [future]

libigl is in active development, and we plan to focus on the following features in the next months:

* A better and more consistent documentation for all functions, plus exteding this tutorial to cover more features of libigl

* Include a robust, adaptive triangular remeshing algorithm. Currently, the only remeshing functions available are only able to create quadrilateral remeshings

* Generate matlab and python wrappers for all libigl functions

* Implement a mixed-integer solver which only uses Eigen to remove the dependency on CoMiSo

* Add a standalone BVH and a simple ray casting engine to make the dependency on Embree optional

We encourage you to contribute to the library and to report problems and bugs that you encounter while using it. The best way of contributing new feature or patches is to fork the libigl repository and to open a [pull request](https://help.github.com/articles/using-pull-requests) on [our github repository](https://github.com/libigl/libigl).



[#botsch_2004]: Matrio Botsch and Leif Kobbelt. "An Intuitive Framework for
Real-Time Freeform Modeling," 2004.
[#jacobson_thesis_2013]: Alec Jacobson,
_Algorithms and Interfaces for Real-Time Deformation of 2D and 3D Shapes_,
2013.
[#jacobson_2011]: Alec Jacobson, Ilya Baran, Jovan Popović, and Olga Sorkine.
["Bounded Biharmonic Weights for Real-Time Deformation,"](https://www.google.com/search?q=Bounded+biharmonic+weights+for+real-time+deformation) 2011.
[#jacobson_mixed_2010]: Alec Jacobson, Elif Tosun, Olga Sorkine, and Denis
Zorin. "Mixed Finite Elements for Variational Surface Modeling," 2010.
[#kavan_2008]: Ladislav Kavan, Steven Collins, Jiri Zara, and Carol O'Sullivan.
"Geometric Skinning with Approximate Dual Quaternion Blending," 2008.
[#kazhdan_2012]: Michael Kazhdan, Jake Solomon, Mirela Ben-Chen,
"Can Mean-Curvature Flow Be Made Non-Singular," 2012.
[#meyer_2003]: Mark Meyer, Mathieu Desbrun, Peter Schröder and Alan H.  Barr,
"Discrete Differential-Geometry Operators for Triangulated
2-Manifolds," 2003.
[#panozzo_2010]: Daniele Panozzo, Enrico Puppo, Luigi Rocca,
"Efficient Multi-scale Curvature and Crease Estimation," 2010.
[#rustamov_2011]: Raid M. Rustamov, "Multiscale Biharmonic Kernels", 2011.
[#sorkine_2004]: Olga Sorkine, Yaron Lipman, Daniel Cohen-Or, Marc Alexa,
Christian Rössl and Hans-Peter Seidel. "Laplacian Surface Editing," 2004.
