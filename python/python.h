#pragma once

#include <pybind/pybind.h>
#include <pybind/operators.h>
#include <pybind/complex.h>
#include <pybind/numpy.h>
#include <pybind/stl.h>
#include <pybind/functional.h>

#include "py_doc.h"

#include <Eigen/Dense>

void assert_is_VectorXd(const std::string name, const Eigen::MatrixXd& v);
void assert_is_RowVectorXd(const std::string name, const Eigen::MatrixXd& v);
void assert_is_Vector3d(const std::string name, const Eigen::MatrixXd& v);
void assert_is_RowVector3d(const std::string name, const Eigen::MatrixXd& v);
void assert_is_Vector4d(const std::string name, const Eigen::MatrixXd& v);
void assert_is_RowVector4d(const std::string name, const Eigen::MatrixXd& v);
void assert_is_Matrix4d(const std::string name, const Eigen::MatrixXd& v);


namespace py = pybind;
