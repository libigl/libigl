//#include <Eigen/Geometry>
//#include <Eigen/Dense>
//#include <Eigen/Sparse>


#include "../python_shared.h"

#include <igl/embree/ambient_occlusion.h>
#include <igl/embree/reorient_facets_raycast.h>
#include <igl/embree/line_mesh_intersection.h>


void python_export_igl_embree(py::module &me) {

  py::module m = me.def_submodule(
    "embree", "Wrappers for libigl functions that use embree");

  #include "../py_igl/embree/py_ambient_occlusion.cpp"
  #include "../py_igl/embree/py_reorient_facets_raycast.cpp"
  #include "../py_igl/embree/py_line_mesh_intersection.cpp"

}
