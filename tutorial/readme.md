xhtml header:   <script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
css: style.css

# Introduction

TODO

# Index

* **100_FileIO**: Example of reading/writing mesh files
* **101_Serialization**: Example of using the XML serialization framework
* **102_DrawMesh**: Example of plotting a mesh
* [202 Gaussian Curvature](#gaus)

# Compilation Instructions

All examples depends on glfw, glew and anttweakbar. A copy
of the sourcecode of each library is provided together with libigl
and they can be precompiled using:

    sh compile_macosx.sh (MACOSX)
    sh compile_linux.sh (LINUX)
    compile_windows.bat (Visual Studio 2012)

Every example can be compiled by using the cmake file provided in its folder.
On Linux and MacOSX, you can use the provided bash script:

    sh ../compile_example.sh

# Chapter 2: Discrete Geometric Quantities and Operators
This chapter illustrates a few discrete quantities that libigl can compute on a
mesh. This also provides an introduction to basic drawing and coloring routines
in our example viewer. Finally, we construct popular discrete differential
geometry operators.

## <a id=gaus></a> Gaussian Curvature
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
at vertex $i$ in triangle $j$.

Just like the continuous analog, our discrete Gaussian curvature reveals
elliptic, hyperbolic and parabolic vertices on the domain.

![The `GaussianCurvature` example computes discrete Gaussian curvature and visualizes it in
pseudocolor.](images/bumpy-gaussian-curvature.jpg)
