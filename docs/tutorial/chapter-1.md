title: libigl Tutorial
author: Daniele Panozzo and Alec Jacobson
date: 07 November 2015
css: style.css
html header:   <script type="text/javascript" src="http://cdn.mathjax.org/mathjax/latest/MathJax.js?config=TeX-AMS-MML_HTMLorMML"></script>
<link rel="stylesheet" href="http://yandex.st/highlightjs/7.3/styles/default.min.css">
<script src="http://yandex.st/highlightjs/7.3/highlight.min.js"></script>
<script>hljs.initHighlightingOnLoad();</script>

# libigl tutorial notes

#### originally presented by Daniele Panozzo and Alec Jacobson at SGP Graduate School 2014

![](images/libigl-logo.jpg)

Libigl is an open source C++ library for geometry processing research and development.  Dropping the heavy data structures of tradition geometry libraries, libigl is a simple header-only library of encapsulated functions. This combines the rapid prototyping familiar to Matlab or Python programmers with the performance and versatility of C++.  The tutorial is a self-contained, hands-on introduction to libigl.  Via interactive, step-by-step examples, we demonstrate how to accomplish common geometry processing tasks such as computation of differential quantities and operators, real-time deformation, parametrization, numerical optimization and remeshing. Each section of the lecture notes links to a cross-platform example application.


# Chapter 1 [chapter1:introductiontolibigl]

We introduce libigl with a series of self-contained examples. The purpose of
each example is to showcase a feature of libigl while applying to a practical
problem in geometry processing. In this chapter, we will present the basic
concepts of libigl and introduce a simple mesh viewer that allows to
visualize a surface mesh and its attributes. All the tutorial examples are
cross-platform and can be compiled on MacOSX, Linux and Windows.

## [libigl design principles](#libigldesignprinciples) [libigldesignprinciples]

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
  library](../optional/))

4. **Function encapsulation.** Every function (including its full
  implementation) is contained in a pair of .h/.cpp files with the same name of
  the function.


### Downloading libigl
libigl can be downloaded from our [github
repository](https://github.com/libigl/libigl) or cloned with git:

```bash
git clone --recursive https://github.com/libigl/libigl.git
```

The core libigl functionality only depends on the C++ Standard Library and
Eigen.

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

*Note for linux users*: Many linux distributions do not include gcc and the basic development tools
in their default installation. On Ubuntu, you need to install the following packages:

```bash
sudo apt-get install git
sudo apt-get install build-essential
sudo apt-get install cmake
sudo apt-get install libx11-dev
sudo apt-get install mesa-common-dev libgl1-mesa-dev libglu1-mesa-dev
sudo apt-get install libxrandr-dev
sudo apt-get install libxi-dev
sudo apt-get install libxmu-dev
sudo apt-get install libblas-dev
sudo apt-get install libxinerama-dev
sudo apt-get install libxcursor-dev
```
*Note for windows users*: libigl only supports the Microsoft Visual Studio 2015 compiler in 64bit mode. It will not work with a 32bit build and it will not work
with older versions of visual studio.

A few examples in Chapter 5 requires the [CoMiSo
solver](http://www.graphics.rwth-aachen.de/software/comiso). We provide a
mirror of CoMISo that works out of the box with libigl. To install it:

```bash
cd libigl/external
git clone --recursive https://github.com/libigl/CoMISo.git
```

You can then build the tutorials again and it libigl will automatically find and
compile CoMISo.

*Note 1*: CoMISo is distributed under the GPL3 license, it does impose restrictions on commercial usage.

*Note 2*: CoMISo requires a blas implementation. We use the built-in blas in macosx and linux, and we bundle a precompiled binary for VS2015 64 bit. Do NOT compile the tutorials
in 32 bit on windows.

### libigl example project

We provide a [blank project example](https://github.com/libigl/libigl-example-project) showing how to use libigl and cmake. Feel free and encouraged to copy or fork this project as a way of starting a new personal project using libigl.

## [Mesh representation](#meshrepresentation) [meshrepresentation]

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
igl::readOFF(TUTORIAL_SHARED_PATH "/cube.off", V, F);
```

The function reads the mesh cube.off and it fills the provided `V` and `F` matrices.
Similarly, a mesh can be written in an OBJ file using:

```cpp
igl::writeOBJ("cube.obj",V,F);
```

[Example 101](101_FileIO/main.cpp) contains a simple mesh
converter from OFF to OBJ format.

## [Visualizing surfaces](#visualizingsurfaces) [visualizingsurfaces]

Libigl provides an glfw-based OpenGL 3.2 viewer to visualize surfaces, their
properties and additional debugging information.

The following code ([Example 102](102_DrawMesh/main.cpp)) is a basic skeleton
for all the examples that will be used in the tutorial.
It is a standalone application that loads a mesh and uses the viewer to
render it.

```cpp
#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>

Eigen::MatrixXd V;
Eigen::MatrixXi F;

int main(int argc, char *argv[])
{
  // Load a mesh in OFF format
  igl::readOFF(TUTORIAL_SHARED_PATH "/bunny.off", V, F);

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  viewer.data().set_mesh(V, F);
  viewer.launch();
}
```

The function `set_mesh` copies the mesh into the viewer.
`Viewer.launch()`  creates a window, an OpenGL context and it starts the draw loop.
Additional properties can be plotted on the mesh (as we will see later),
and it is possible to extend the viewer with standard OpenGL code.
Please see the documentation in
[Viewer.h](../include/igl/opengl/glfw/Viewer.h) for more details.

![([Example 102](102_DrawMesh/main.cpp)) loads and draws a
mesh.](images/102_DrawMesh.png)

## [Interaction with keyboard and mouse](#interactionwithkeyboardandmouse) [interactionwithkeyboardandmouse]

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
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  if (key == '1')
  {
    viewer.data().clear();
    viewer.data().set_mesh(V1, F1);
    viewer.core.align_camera_center(V1,F1);
  }
  else if (key == '2')
  {
    viewer.data().clear();
    viewer.data().set_mesh(V2, F2);
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
[Viewer_plugin](../include/igl/opengl/glfw/ViewerPlugin.h) for more details.

## [Scalar field visualization](#scalarfieldvisualization) [scalarfieldvisualization]

Colors and normals can be associated to faces or vertices using the
set_colors function:

```cpp
viewer.data().set_colors(C);
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

## [Overlays](#overlays) [overlays]

In addition to plotting the surface, the viewer supports the visualization of points, lines and text labels: these overlays can be very helpful while developing geometric processing algorithms to plot debug information.

```cpp
viewer.data().add_points(P,Eigen::RowVector3d(r,g,b));
```

Draws a point of color r,g,b for each row of P. The point is placed at the coordinates specified in each row of P, which is a #P by 3 matrix.

```cpp
viewer.data().add_edges(P1,P2,Eigen::RowVector3d(r,g,b);
```

Draws a line of color r,g,b for each row of P1 and P2, which connects the 3D point in to the point in P2. Both P1 and P2 are of size #P by 3.

```cpp
viewer.data().add_label(p,str);
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

## [Viewer Menu](#viewermenu) [viewermenu]

As of latest version, the viewer uses a new menu and completely replaces
[AntTweakBar](http://anttweakbar.sourceforge.net/doc/) and
[nanogui](https://github.com/wjakob/nanogui) with [Dear ImGui](https://github.com/ocornut/imgui). To extend the default menu of the
viewer and to expose more user defined variables you have to implement a custom interface, as in [Example 106](106_ViewerMenu/main.cpp):
```cpp
// Add content to the default menu window
menu.callback_draw_viewer_menu = [&]()
{
  // Draw parent menu content
  menu.draw_viewer_menu();

  // Add new group
  if (ImGui::CollapsingHeader("New Group", ImGuiTreeNodeFlags_DefaultOpen))
  {
    // Expose variable directly ...
    ImGui::InputFloat("float", &floatVariable, 0, 0, 3);

    // ... or using a custom callback
    static bool boolVariable = true;
    if (ImGui::Checkbox("bool", &boolVariable))
    {
      // do something
      std::cout << "boolVariable: " << std::boolalpha << boolVariable << std::endl;
    }

    // Expose an enumeration type
    enum Orientation { Up=0, Down, Left, Right };
    static Orientation dir = Up;
    ImGui::Combo("Direction", (int *)(&dir), "Up\0Down\0Left\0Right\0\0");

    // We can also use a std::vector<std::string> defined dynamically
    static int num_choices = 3;
    static std::vector<std::string> choices;
    static int idx_choice = 0;
    if (ImGui::InputInt("Num letters", &num_choices))
    {
      num_choices = std::max(1, std::min(26, num_choices));
    }
    if (num_choices != (int) choices.size())
    {
      choices.resize(num_choices);
      for (int i = 0; i < num_choices; ++i)
        choices[i] = std::string(1, 'A' + i);
      if (idx_choice >= num_choices)
        idx_choice = num_choices - 1;
    }
    ImGui::Combo("Letter", &idx_choice, choices);

    // Add a button
    if (ImGui::Button("Print Hello", ImVec2(-1,0)))
    {
      std::cout << "Hello\n";
    }
  }
};
```

If you need a separate new menu window implement:

```cpp
// Draw additional windows
menu.callback_draw_custom_window = [&]()
{
  // Define next window position + size
  ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10), ImGuiSetCond_FirstUseEver);
  ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiSetCond_FirstUseEver);
  ImGui::Begin(
      "New Window", nullptr,
      ImGuiWindowFlags_NoSavedSettings
  );

  // Expose the same variable directly ...
  ImGui::PushItemWidth(-80);
  ImGui::DragFloat("float", &floatVariable, 0.0, 0.0, 3.0);
  ImGui::PopItemWidth();

  static std::string str = "bunny";
  ImGui::InputText("Name", str);

  ImGui::End();
};
```

![([Example 106](106_ViewerMenu/main.cpp)) The UI of the viewer can be easily
customized.](images/106_ViewerMenu.png)

## [Multiple Meshes](#multiplemeshes) [multiplemeshes]

Libigl's `igl::opengl::glfw::Viewer` provides basic support for rendering
multiple meshes.

Which mesh is _selected_ is controlled via the `viewer.selected_data_index`
field. By default it his is set to `0`, so in the typical case of a single mesh
`viewer.data()` returns the `igl::ViewerData` corresponding to the one
and only mesh.

![([Example 107](107_MultipleMeshes/main.cpp)) The `igl::opengl::glfw::Viewer`
can render multiple meshes, each with its own attributes like
colors.](images/multiple-meshes.png)
