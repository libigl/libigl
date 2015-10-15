#include "python.h"
#include <sstream>
#include <string>
#include <fstream>

extern void python_export_vector(py::module &);
extern void python_export_igl(py::module &);

#ifdef PY_VIEWER
extern void python_export_igl_viewer(py::module &);
#endif

#ifdef PY_COMISO
extern void python_export_igl_comiso(py::module &);
#endif

PYTHON_PLUGIN(igl) {
    py::init_threading();
    py::module m("igl", "Python wrappers for libigl");

    python_export_vector(m);
    python_export_igl(m);


    #ifdef PY_VIEWER
    python_export_igl_viewer(m);
    #endif

    #ifdef PY_COMISO
    python_export_igl_comiso(m);
    #endif

    return m.ptr();
}
