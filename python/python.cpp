#include "python.h"
#include <sstream>
#include <string>
#include <fstream>

extern void python_export_vector(py::module &);
extern void python_export_igl(py::module &);

PYTHON_PLUGIN(igl) {
    py::init_threading();
    py::module m("igl", "Python wrappers for libigl");

    python_export_vector(m);
    python_export_igl(m);

    return m.ptr();
}
