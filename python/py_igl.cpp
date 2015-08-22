#include <Eigen/Dense>

#include "python.h"

#include <igl/readOFF.h>
#include <igl/writeOBJ.h>
#include <igl/per_face_normals.h>
#include <igl/per_corner_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/gaussian_curvature.h>
#include <igl/jet.h>

void python_export_igl(py::module &m)
{
#include "py_igl/py_readOFF.cpp"
#include "py_igl/py_writeOBJ.cpp"
#include "py_igl/py_per_face_normals.cpp"
#include "py_igl/py_per_corner_normals.cpp"
#include "py_igl/py_per_vertex_normals.cpp"
#include "py_igl/py_gaussian_curvature.cpp"
#include "py_igl/py_jet.cpp"
}
