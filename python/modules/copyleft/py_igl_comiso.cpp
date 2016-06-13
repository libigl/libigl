#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/Sparse>


#include "../../python_shared.h"

#include <igl/copyleft/comiso/nrosy.h>
#include <igl/copyleft/comiso/miq.h>

void python_export_igl_comiso(py::module &me) {

  py::module m = me.def_submodule(
    "comiso", "Wrappers for libigl functions that use comiso");

  #include "../../py_igl/copyleft/comiso/py_nrosy.cpp"
  #include "../../py_igl/copyleft/comiso/py_miq.cpp"

}
