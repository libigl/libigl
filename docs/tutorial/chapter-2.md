# Chapter 2: Discrete Geometric Quantities and Operators

This chapter illustrates a few discrete quantities that libigl can compute on a
mesh and the libigl functions that construct popular discrete differential
geometry operators. It also provides an introduction to basic drawing and
coloring routines of our viewer.

## Normals
Surface normals are a basic quantity necessary for rendering a surface. There
are a variety of ways to compute and store normals on a triangle mesh. [Example
201]({{ repo_url }}/tutorial/201_Normals/main.cpp) demonstrates how to compute and visualize normals
with libigl.

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

![The `Normals` example computes per-face (left), per-vertex (middle) and per-corner (right) normals](images/fandisk-normals.jpg)

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
at vertex $i$ in triangle $j$ [^meyer_2003].

Just like the continuous analog, our discrete Gaussian curvature reveals
elliptic, hyperbolic and parabolic vertices on the domain, as demonstrated in [Example 202]({{ repo_url }}/tutorial/202_GaussianCurvature/main.cpp).

![The `GaussianCurvature` example computes discrete Gaussian curvature and visualizes it in pseudocolor.](images/bumpy-gaussian-curvature.jpg)

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
cotangent Laplace-Beltrami operator [^meyer_2003].

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
directions [^meyer_2003].

Alternatively, a robust method for determining principal curvatures is via
quadric fitting [^panozzo_2010]. In the neighborhood around every vertex, a
best-fit quadric is found and principal curvature values and directions are
analytically computed on this quadric ([Example
203]({{ repo_url }}/tutorial/203_curvatureDirections/main.cpp)).

![The `CurvatureDirections` example computes principal curvatures via quadric fitting and visualizes mean curvature in pseudocolor and principal directions with a cross field.](images/fertility-principal-curvature.jpg)

## Gradient
Scalar functions on a surface can be discretized as a piecewise linear function
with values defined at each mesh vertex:

 $f(\mathbf{x}) \approx \sum\limits_{i=1}^n \phi_i(\mathbf{x})\, f_i,$

where $\phi_i$ is a piecewise linear hat function defined by the mesh so that
for each triangle $\phi_i$ is _the_ linear function which is one only at
vertex $i$ and zero at the other corners.

![Hat function $\phi_i$ is one at vertex $i$, zero at all other vertices, and linear on incident triangles.](images/hat-function.jpg)

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
ch. 2[^jacobson_thesis_2013].
Libigl's `grad` function computes $\mathbf{G}$ for
triangle and tetrahedral meshes ([Example 204]({{ repo_url }}/tutorial/204_Gradient/main.cpp)):

![The `Gradient` example computes gradients of an input function on a mesh and visualizes the vector field.](images/cheburashka-gradient.jpg)

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
[^sharf_2007].

The operator applied to mesh vertex positions amounts to smoothing by _flowing_
the surface along the mean curvature normal direction ([Example 205]({{ repo_url }}/tutorial/205_Laplacian/main.cpp)). Note that this is equivalent to minimizing surface area.

![The `Laplacian` example computes conformalized mean curvature flow using the cotangent Laplacian [^kazhdan_2012].](images/cow-curvature-flow.jpg)

### Mass matrix
The mass matrix $\mathbf{M}$ is another $n \times n$ matrix which takes vertex
values to vertex values. From an FEM point of view, it is a discretization of
the inner-product: it accounts for the area around each vertex. Consequently,
$\mathbf{M}$ is often a diagonal matrix, such that $M_{ii}$ is the barycentric
or voronoi area around vertex $i$ in the mesh [^meyer_2003]. The inverse of
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

## Geodesic

The discrete geodesic distance between two points is the length of the shortest path between then restricted to the surface. For triangle meshes, such a path is made of a set of segments which can be either edges of the mesh or crossing a triangle.

Libigl includes a wrapper for the exact geodesic algorithm [^mitchell_1987] developed by Danil Kirsanov (https://code.google.com/archive/p/geodesic/), exposing it through an Eigen-based API. The function 
```cpp
igl::exact_geodesic(V,F,VS,FS,VT,FT,d);
```
computes the closest geodesic distances of each vertex in VT or face in FT, from the source vertices VS or faces FS of the input mesh V,F. The output is writted in the vector d, which lists first the distances for the vertices in VT, and then for the faces in FT. For example, if you want to compute the distance from the vertex with id ```vid```, to all vertices of F you can use:
```cpp
Eigen::VectorXi VS,FS,VT,FT;
// The selected vertex is the source
VS.resize(1);
VS << vid;
// All vertices are the targets
VT.setLinSpaced(V.rows(),0,V.rows()-1);
Eigen::VectorXd d;
igl::exact_geodesic(V,F,VS,FS,VT,FT,d);
```

![[Example 206]({{ repo_url }}/tutorial/206_GeodesicDistance/main.cpp) allows to interactively pick the source vertex and displays the distance using a periodic color pattern.](images/geodesicdistance.jpg)

## References

[^jacobson_thesis_2013]: Alec Jacobson, [_Algorithms and Interfaces for Real-Time Deformation of 2D and 3D Shapes_](https://www.google.com/search?q=Algorithms+and+Interfaces+for+Real-Time+Deformation+of+2D+and+3D+Shapes), 2013.
[^kazhdan_2012]: Michael Kazhdan, Jake Solomon, Mirela Ben-Chen, [Can Mean-Curvature Flow Be Made Non-Singular](https://www.google.com/search?q=Can+Mean-Curvature+Flow+Be+Made+Non-Singular), 2012.
[^meyer_2003]: Mark Meyer, Mathieu Desbrun, Peter Schröder and Alan H.  Barr, [Discrete Differential-Geometry Operators for Triangulated 2-Manifolds](https://www.google.com/search?q=Discrete+Differential-Geometry+Operators+for+Triangulated+2-Manifolds), 2003.
[^mitchell_1987]: Joseph S. B. Mitchell, David M. Mount, Christos H. Papadimitriou. [The Discrete Geodesic Problem](https://www.google.com/search?q=The+Discrete+Geodesic+Problem), 1987
[^panozzo_2010]: Daniele Panozzo, Enrico Puppo, Luigi Rocca, [Efficient Multi-scale Curvature and Crease Estimation](https://www.google.com/search?q=Efficient+Multi-scale+Curvature+and+Crease+Estimation), 2010.
[^sharf_2007]: Andrei Sharf, Thomas Lewiner, Gil Shklarski, Sivan Toledo, and Daniel Cohen-Or. [Interactive topology-aware surface reconstruction](https://www.google.com/search?q=Interactive+topology-aware+surface+reconstruction), 2007.
