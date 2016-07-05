//#include <Eigen/Geometry>
//#include <Eigen/Dense>
//#include <Eigen/Sparse>


#include "../../python_shared.h"

#include <igl/copyleft/cgal/mesh_boolean.h>


void python_export_igl_cgal(py::module &me) {

  py::module m = me.def_submodule(
    "cgal", "Wrappers for libigl functions that use cgal");

  #include "../../py_igl/copyleft/cgal/py_mesh_boolean.cpp"

}
