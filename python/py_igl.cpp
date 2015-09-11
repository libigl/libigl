#include <Eigen/Dense>

#include "python.h"

#include <igl/readOFF.h>
#include <igl/writeOBJ.h>
#include <igl/per_face_normals.h>
#include <igl/per_corner_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/gaussian_curvature.h>
#include <igl/jet.h>
#include <igl/read_triangle_mesh.h>
#include <igl/cotmatrix.h>
#include <igl/massmatrix.h>
#include <igl/invert_diag.h>
#include <igl/principal_curvature.h>
#include <igl/parula.h>
#include <igl/readDMAT.h>
#include <igl/grad.h>
#include <igl/avg_edge_length.h>
#include <igl/barycenter.h>
#include <igl/doublearea.h>

void python_export_igl(py::module &m)
{
#include "py_igl/py_readOFF.cpp"
#include "py_igl/py_writeOBJ.cpp"
#include "py_igl/py_per_face_normals.cpp"
#include "py_igl/py_per_corner_normals.cpp"
#include "py_igl/py_per_vertex_normals.cpp"
#include "py_igl/py_gaussian_curvature.cpp"
#include "py_igl/py_jet.cpp"
#include "py_igl/py_read_triangle_mesh.cpp"
#include "py_igl/py_cotmatrix.cpp"
#include "py_igl/py_massmatrix.cpp"
#include "py_igl/py_invert_diag.cpp"
#include "py_igl/py_principal_curvature.cpp"
#include "py_igl/py_parula.cpp"
#include "py_igl/py_readDMAT.cpp"
#include "py_igl/py_grad.cpp"
#include "py_igl/py_avg_edge_length.cpp"
#include "py_igl/py_barycenter.cpp"
#include "py_igl/py_doublearea.cpp"

}
