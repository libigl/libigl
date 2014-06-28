title: libigl Tutorial
author: Daniele Panozzo, Alec Jacobson and others
date: 20 June 2014
css: style.css
html header:   <script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
<link rel="stylesheet" href="http://yandex.st/highlightjs/7.3/styles/default.min.css">
<script src="http://yandex.st/highlightjs/7.3/highlight.min.js"></script>
<script>hljs.initHighlightingOnLoad();</script>

* [Chapter 1: Introduction to libigl][100]
    * [Mesh representation][101]
    * [Plotting surfaces][102]
    * [Interaction with keyboard and mouse][103]
    * [Scalar field visualization][104]
    * [Overlays][105]
    * [Picking vertices and faces][106]
    * [libigl design principles][107]

# Chapter 1 [100]

We introduce libIGL with a series of self-contained examples. The purpose of each example is to showcase a feature of libIGL while applying to a practical problem in geometry processing. In this chapter, we will showcase the basic concepts of libigl and introduce a simple mesh viewer that allows to easily visualize surface mesh and its attributes. All the examples are cross-platform and can be compiled on MacOSX, Linux and Windows.

All dependencies for the compilation of these examples are contained in libigl (external folder), with the exception of Eigen, which should be downloaded and unpacked in the folder containing the libigl root folder.

All examples depends on glfw, glew and anttweakbar. A copy
of the sourcecode of each library is provided together with libigl
and they can be precompiled using:
```sh
    sh compile_dependencies_macosx.sh (MACOSX)
    sh compile_dependencies_linux.sh (LINUX)
```
Precompiled binaries are provided for Visual Studio 2014 64bit.

Use the cmake file in the tutorial folder to build all the examples:
```sh
  cd tutorial
  mkdir build
  cd build
  cmake ../
  make
```

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
    viewer.clear_mesh();
    viewer.set_mesh(V1, F1);
  }
  else if (key == '2')
  {
    viewer.clear_mesh();
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
