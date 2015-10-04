#include <Eigen/Geometry>
#include <Eigen/Dense>
#include <Eigen/Sparse>


#include "python.h"

#include <igl/comiso/nrosy.h>
#include <igl/comiso/miq.h>

void python_export_igl_comiso(py::module &me) {

  py::module m = me.def_submodule(
    "comiso", "Wrappers for libigl functions that use comiso");

  #include "py_igl/comiso/py_nrosy.cpp"
  #include "py_igl/comiso/py_miq.cpp"

}
