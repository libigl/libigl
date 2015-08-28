#include "python.h"
#include <sstream>
#include <string>
#include <fstream>

void assert_is_VectorXd(const std::string name, const Eigen::MatrixXd& v)
{
  if (v.cols() != 1)
    throw std::runtime_error(name + " must be a column vector.");
}

void assert_is_RowVectorXd(const std::string name, const Eigen::MatrixXd& v)
{
  if (v.rows() != 1)
    throw std::runtime_error(name + " must be a row vector.");
}

void assert_is_Vector3d(const std::string name, const Eigen::MatrixXd& v)
{
  if ((v.cols() != 1) || (v.rows() != 3))
    throw std::runtime_error(name + " must be a column vector with 3 entries.");
}

void assert_is_RowVector3d(const std::string name, const Eigen::MatrixXd& v)
{
  if ((v.cols() != 3) || (v.rows() != 1))
    throw std::runtime_error(name + " must be a row vector with 3 entries.");
}

extern void python_export_vector(py::module &);
extern void python_export_igl(py::module &);
extern void python_export_igl_viewer(py::module &);

PYTHON_PLUGIN(igl) {
    py::init_threading();
    py::module m("igl", "Python wrappers for libigl");

    python_export_vector(m);
    python_export_igl(m);
    python_export_igl_viewer(m);

    return m.ptr();
}
