# libigl - A simple C++ geometry processing library

This detailed documentation browser is automatically generated from the comments
in libigl header (.h) files.

In general, each [libigl function](./namespaceigl.html#func-members) (e.g., `igl::func`) will be defined in a
correspondingly named header file (e.g., `#include <igl/func.h>`).

The _core_ library only depends on the standard template library (`std::`) and 
Eigen. These functions reside directly the [`igl::` namespace](./namespaceigl.html)

Functions with further dependencies reside in a corresonding sub-namespace. For
example, the function `igl::spectra::lscm` depends on the Spectra library so it
resides in the [`igl::spectra::` namespace](./namespaceigl_1_1spectra.html).

Functions which depend on external code under a copyleft license reside in the
[`igl::copyleft::` namepsace](file:///Users/alecjacobson/Repos/libigl/dox/namespaceigl_1_1copyleft.html).


Most libigl functions are templated over the Eigen matrix inputs and outputs.
Callers can choose their own scalar types (e.g., `double`/`float`) and storage
orders (`Eigen::ColMajor`/`Eigen::RowMajor`). Libigl can be used as a:

 - **header only library** (via CMake, make sure
   `LIBIGL_USE_STATIC_LIBRARY=OFF`) and insure that `IGL_STATIC_LIBRARY` is
   _not_ defined_ when compiling --- easiest if you're new to libigl, or 
 - **static library** (`LIBIGL_USE_STATIC_LIBRARY=ON` → `IGL_STATIC_LIBRARY` is
  defined) --- speeds up repeated compilation.

The libigl static library is filled with _explicit template instantiations_ for
common Eigen inputs and outputs. If the library doesn't contain your types, you
may get some form of linker error (e.g., `Undefined symbols for architecture`,
`undefined reference to` or `unresolved external symbol`).

You can fix this by:

1. Switching to header only mode for your project,
2. Making  a file in your project to compile the missing templates. E.g., `my_templates.cpp` 
```cpp
#ifdef IGL_STATIC_LIBRARY
#undef IGL_STATIC_LIBRARY
#endif
#include <igl/per_vertex_normals.h>
template void igl::per_vertex_normals<Eigen::Matrix<float, -1, -1, 0, -1, -1>, Eigen::Matrix<int, -1, -1, 0, -1, -1>, Eigen::Matrix<float, -1, -1, 0, -1, -1> >(Eigen::MatrixBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> > const&, Eigen::MatrixBase<Eigen::Matrix<int, -1, -1, 0, -1, -1> > const&, Eigen::PlainObjectBase<Eigen::Matrix<float, -1, -1, 0, -1, -1> >&); 
```
3. [Submit a PR](https://github.com/libigl/libigl/pulls) containing your missing template to the development branch of libigl
4. Change your input/output types to match existing templates (e.g., `Eigen::MatrixXd`).




https://libigl.github.io/

https://github.com/libigl/libigl/
