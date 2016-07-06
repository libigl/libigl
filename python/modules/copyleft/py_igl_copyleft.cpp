//#include <Eigen/Geometry>
//#include <Eigen/Dense>
//#include <Eigen/Sparse>


#include "../../python_shared.h"

#include <igl/copyleft/marching_cubes.h>
#include <igl/copyleft/swept_volume.h>


void python_export_igl_copyleft(py::module &me) {

  py::module m = me.def_submodule(
    "copyleft", "Wrappers for libigl functions that are copyleft");

  #include "../../py_igl/copyleft/py_marching_cubes.cpp"
  #include "../../py_igl/copyleft/py_swept_volume.cpp"

}
