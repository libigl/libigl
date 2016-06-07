//#include <Eigen/Geometry>
//#include <Eigen/Dense>
//#include <Eigen/Sparse>


#include "../../python_shared.h"

#include <igl/copyleft/tetgen/tetrahedralize.h>


void python_export_igl_tetgen(py::module &me) {

  py::module m = me.def_submodule(
    "tetgen", "Wrappers for libigl functions that use tetgen");

  #include "../../py_igl/copyleft/tetgen/py_tetrahedralize.cpp"

}
