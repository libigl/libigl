//#include <Eigen/Geometry>
//#include <Eigen/Dense>
//#include <Eigen/Sparse>


#include "../../python_shared.h"

#include <igl/copyleft/cgal/mesh_boolean.h>
#include <igl/copyleft/cgal/remesh_self_intersections.h>
#include <igl/copyleft/cgal/RemeshSelfIntersectionsParam.h>


void python_export_igl_cgal(py::module &me) {

  py::module m = me.def_submodule(
    "cgal", "Wrappers for libigl functions that use cgal");

  #include "../../py_igl/copyleft/cgal/py_mesh_boolean.cpp"
  #include "../../py_igl/copyleft/cgal/py_remesh_self_intersections.cpp"
  #include "../../py_igl/copyleft/cgal/py_RemeshSelfIntersectionsParam.cpp"

}
