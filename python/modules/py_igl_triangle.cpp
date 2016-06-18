//#include <Eigen/Geometry>
//#include <Eigen/Dense>
//#include <Eigen/Sparse>


#include "../python_shared.h"

#include <igl/triangle/triangulate.h>


void python_export_igl_triangle(py::module &me) {

  py::module m = me.def_submodule(
    "triangle", "Wrappers for libigl functions that use triangle");

  #include "../py_igl/triangle/py_triangulate.cpp"

}
