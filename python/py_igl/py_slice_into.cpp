// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("slice_into", []
(
  const Eigen::SparseMatrix<double>& X,
  const Eigen::MatrixXi& R,
  const Eigen::MatrixXi& C,
  Eigen::SparseMatrix<double>& Y
)
{
  assert_is_VectorX("R",R);
  assert_is_VectorX("C",C);
  return igl::slice_into(X,R,C,Y);
}, __doc_igl_slice_into,
py::arg("X"), py::arg("R"), py::arg("C"), py::arg("Y"));

m.def("slice_into", []
(
  const Eigen::MatrixXd& X,
  const Eigen::MatrixXi& R,
  const Eigen::MatrixXi& C,
  Eigen::MatrixXd& Y
)
{
  assert_is_VectorX("R",R);
  assert_is_VectorX("C",C);
  return igl::slice_into(X,R,C,Y);
}, __doc_igl_slice_into,
py::arg("X"), py::arg("R"), py::arg("C"), py::arg("Y"));

m.def("slice_into", []
(
  const Eigen::MatrixXd& X,
  const Eigen::MatrixXi& R,
  const int& dim,
  Eigen::MatrixXd& Y
)
{
  assert_is_VectorX("R",R);
  return igl::slice_into(X,R,dim,Y);
}, __doc_igl_slice_into,
py::arg("X"), py::arg("R"), py::arg("dim"), py::arg("Y"));

m.def("slice_into", []
(
  const Eigen::MatrixXd& X,
  const Eigen::MatrixXi& R,
  Eigen::MatrixXd& Y
)
{
  assert_is_VectorX("R",R);
  return igl::slice_into(X,R,Y);
}, __doc_igl_slice_into,
py::arg("X"), py::arg("R"), py::arg("Y"));

// int

m.def("slice_into", []
(
  const Eigen::SparseMatrix<int>& X,
  const Eigen::MatrixXi& R,
  const Eigen::MatrixXi& C,
  Eigen::SparseMatrix<int>& Y
)
{
  assert_is_VectorX("R",R);
  assert_is_VectorX("C",C);
  return igl::slice_into(X,R,C,Y);
}, __doc_igl_slice_into,
py::arg("X"), py::arg("R"), py::arg("C"), py::arg("Y"));

m.def("slice_into", []
(
  const Eigen::MatrixXi& X,
  const Eigen::MatrixXi& R,
  const Eigen::MatrixXi& C,
  Eigen::MatrixXi& Y
)
{
  assert_is_VectorX("R",R);
  assert_is_VectorX("C",C);
  return igl::slice_into(X,R,C,Y);
}, __doc_igl_slice_into,
py::arg("X"), py::arg("R"), py::arg("C"), py::arg("Y"));

m.def("slice_into", []
(
  const Eigen::MatrixXi& X,
  const Eigen::MatrixXi& R,
  const int& dim,
  Eigen::MatrixXi& Y
)
{
  assert_is_VectorX("R",R);
  return igl::slice_into(X,R,dim,Y);
}, __doc_igl_slice_into,
py::arg("X"), py::arg("R"), py::arg("dim"), py::arg("Y"));

m.def("slice_into", []
(
  const Eigen::MatrixXi& X,
  const Eigen::MatrixXi& R,
  Eigen::MatrixXi& Y
)
{
  assert_is_VectorX("R",R);
  return igl::slice_into(X,R,Y);
}, __doc_igl_slice_into,
py::arg("X"), py::arg("R"), py::arg("Y"));
