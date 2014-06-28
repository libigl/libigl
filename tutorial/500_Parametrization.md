title: libigl Tutorial
author: Daniele Panozzo, Alec Jacobson and others
date: 20 June 2014
css: style.css
html header:   <script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
<link rel="stylesheet" href="http://yandex.st/highlightjs/7.3/styles/default.min.css">
<script src="http://yandex.st/highlightjs/7.3/highlight.min.js"></script>
<script>hljs.initHighlightingOnLoad();</script>

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

creates a smooth vector field (N=1) starting from a sparse set of constrained faces, whose indices are listed in b and their constrained value is specified in bc. The functions supports soft_constraints (b_soft,b_soft_weight,bc_soft), and returns the interpolated field for each face of the triangle mesh (output_field) plus the singularities of the field (output_field).

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

## Gradient field design (up to rotation and trianslation) [505]

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

[Frame Fields: Anisotropic and Non-Orthogonal Cross Fields],
Daniele Panozzo, Enrico Puppo, Marco Tarini, Olga Sorkine-Hornung,
SIGGRAPH, 2014

## N-PolyVector fields [507]

* further generalization to arbitrary rosy, same interface

* globally optimal and keenan optimal field are a subset of them

## Conjugate vector fields [508]

* they can be used to encode conjugate field -> planar meshing

* global/local approach

## Planarization [509]

* given a mesh from conjugate directions, enforce planarity with a local/global approach
* useful for architecture
