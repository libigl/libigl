// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

// Wrap the data class, no properties are exposed since it is not necessary
py::class_<igl::min_quad_with_fixed_data<double> > min_quad_with_fixed_data(m, "min_quad_with_fixed_data");

min_quad_with_fixed_data
.def(py::init<>());

m.def("min_quad_with_fixed_precompute", []
(
  const Eigen::SparseMatrix<double>& A,
  const Eigen::MatrixXi& known,
  const Eigen::SparseMatrix<double>& Aeq,
  const bool pd,
  igl::min_quad_with_fixed_data<double> & data
)
{
  assert_is_VectorX("known",known);
  return igl::min_quad_with_fixed_precompute(A,known,Aeq,pd,data);
}, __doc_igl_min_quad_with_fixed,
py::arg("A"), py::arg("known"), py::arg("Aeq"), py::arg("pd"), py::arg("data"));

m.def("min_quad_with_fixed_solve", []
(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd& B,
  const Eigen::MatrixXd& Y,
  const Eigen::MatrixXd & Beq,
  Eigen::MatrixXd& Z,
  Eigen::MatrixXd& sol
)
{
  assert_is_VectorX("B",B);
  assert_is_VectorX("Y",Y);
  assert_is_VectorX("Beq",Beq);
  return igl::min_quad_with_fixed_solve(data,B,Y,Beq,Z,sol);
}, __doc_igl_min_quad_with_fixed,
py::arg("data"), py::arg("B"), py::arg("Y"), py::arg("Beq"), py::arg("Z"), py::arg("sol"));

m.def("min_quad_with_fixed_solve", []
(
  const igl::min_quad_with_fixed_data<double> & data,
  const Eigen::MatrixXd& B,
  const Eigen::MatrixXd& Y,
  const Eigen::MatrixXd & Beq,
  Eigen::MatrixXd& Z
)
{
  assert_is_VectorX("B",B);
  assert_is_VectorX("Y",Y);
  assert_is_VectorX("Beq",Beq);
  return igl::min_quad_with_fixed_solve(data,B,Y,Beq,Z);
}, __doc_igl_min_quad_with_fixed,
py::arg("data"), py::arg("B"), py::arg("Y"), py::arg("Beq"), py::arg("Z"));

m.def("min_quad_with_fixed", []
(
  const Eigen::SparseMatrix<double>& A,
  const Eigen::MatrixXd& B,
  const Eigen::MatrixXi& known,
  const Eigen::MatrixXd& Y,
  const Eigen::SparseMatrix<double>& Aeq,
  const Eigen::MatrixXd& Beq,
  const bool pd,
  Eigen::MatrixXd& Z
)
{
  assert_is_VectorX("B",B);
  assert_is_VectorX("known",known);
  assert_is_VectorX("Y",Y);
  assert_is_VectorX("Beq",Beq);
  return igl::min_quad_with_fixed(A,B,known,Y,Aeq,Beq,pd,Z);
}, __doc_igl_min_quad_with_fixed,
py::arg("A"), py::arg("B"), py::arg("known"), py::arg("Y"), py::arg("Aeq"), py::arg("Beq"), py::arg("pd"), py::arg("Z"));
