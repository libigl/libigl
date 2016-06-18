
#include "../python_shared.h"

#include <igl/png/readPNG.h>
#include <igl/png/writePNG.h>


void python_export_igl_png(py::module &me) {

  py::module m = me.def_submodule(
    "png", "Wrappers for libigl functions that use png");

  #include "../py_igl/png/py_readPNG.cpp"
  #include "../py_igl/png/py_writePNG.cpp"

}
