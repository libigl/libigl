#include <Eigen/Dense>

#include "python_shared.h"
#include "modules/py_typedefs.h"

% for f in functions:
#include <igl/${f}.h>
% endfor


void python_export_igl(py::module &m)
{
#include "modules/py_typedefs.cpp"

% for f in functions:
#include "py_igl/py_${f}.cpp"
% endfor

}
