//#include <Eigen/Geometry>
//#include <Eigen/Dense>
//#include <Eigen/Sparse>


#include "../python_shared.h"

#include <igl/bbw.h>


void python_export_igl_bbw(py::module &me) {

  py::module m = me.def_submodule(
    "bbw", "Wrappers for libigl functions that use bbw");

  #include "../py_igl/py_bbw.cpp"

}
