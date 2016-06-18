#include <Eigen/Dense>

#include "python_shared.h"

% for f in functions:
#include <igl/${f}.h>
% endfor


void python_export_igl(py::module &m)
{
% for f in functions:
#include "py_igl/py_${f}.cpp"
% endfor

}
