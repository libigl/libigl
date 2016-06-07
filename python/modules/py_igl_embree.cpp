//#include <Eigen/Geometry>
//#include <Eigen/Dense>
//#include <Eigen/Sparse>


#include "../python_shared.h"

#include <igl/embree/ambient_occlusion.h>


void python_export_igl_embree(py::module &me) {

  py::module m = me.def_submodule(
    "embree", "Wrappers for libigl functions that use embree");

  #include "../py_igl/embree/py_ambient_occlusion.cpp"

}
