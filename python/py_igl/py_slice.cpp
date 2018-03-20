// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
// // Double
//
m.def("slice", []
(
  const Eigen::SparseMatrix<double>& X,
  const Eigen::MatrixXi& R,
  const Eigen::MatrixXi& C,
  Eigen::SparseMatrix<double>& Y
)
{
  assert_is_VectorX("R",R);
  assert_is_VectorX("C",C);
  return igl::slice(X,R,C,Y);
}, __doc_igl_slice,
py::arg("X"), py::arg("R"), py::arg("C"), py::arg("Y"));

m.def("slice", []
(
  const Eigen::SparseMatrix<double>& X,
  const Eigen::MatrixXi& R,
  const int& dim,
  Eigen::SparseMatrix<double>& Y
)
{
  assert_is_VectorX("R",R);
  return igl::slice(X,R,dim,Y);
}, __doc_igl_slice,
py::arg("X"), py::arg("R"), py::arg("dim"), py::arg("Y"));

m.def("slice", []
(
  const Eigen::MatrixXd& X,
  const Eigen::MatrixXi& R,
  const Eigen::MatrixXi& C,
  Eigen::MatrixXd& Y
)
{
  assert_is_VectorX("R",R);
  assert_is_VectorX("C",C);
  return igl::slice(X,R,C,Y);
}, __doc_igl_slice,
py::arg("X"), py::arg("R"), py::arg("C"), py::arg("Y"));

m.def("slice", []
(
  const Eigen::MatrixXd& X,
  const Eigen::MatrixXi& R,
  const int dim,
  Eigen::MatrixXd& Y
)
{
  assert_is_VectorX("R",R);
  return igl::slice(X,R,dim,Y);
}, __doc_igl_slice,
py::arg("X"), py::arg("R"), py::arg("dim"), py::arg("Y"));

m.def("slice", []
(
  const Eigen::MatrixXd& X,
  const Eigen::MatrixXi& R,
  Eigen::MatrixXd& Y
)
{
  assert_is_VectorX("R",R);
  return igl::slice(X,R,Y);
}, __doc_igl_slice,
py::arg("X"), py::arg("R"), py::arg("Y"));

m.def("slice", []
(
  const Eigen::MatrixXd& X,
  const Eigen::MatrixXi& R
)
{
  assert_is_VectorX("R",R);
  return igl::slice(X,R);
}, __doc_igl_slice,
py::arg("X"), py::arg("A"));

m.def("slice", []
(
  const Eigen::MatrixXd& X,
  const Eigen::MatrixXi& R,
  const int& dim
)
{
  assert_is_VectorX("R",R);
  return igl::slice(X,R,dim);
}, __doc_igl_slice,
py::arg("X"), py::arg("R"), py::arg("dim"));

// int
m.def("slice", []
(
  const Eigen::SparseMatrix<int>& X,
  const Eigen::MatrixXi& R,
  const Eigen::MatrixXi& C,
  Eigen::SparseMatrix<int>& Y
)
{
  assert_is_VectorX("R",R);
  assert_is_VectorX("C",C);
  return igl::slice(X,R,C,Y);
}, __doc_igl_slice,
py::arg("X"), py::arg("R"), py::arg("C"), py::arg("Y"));

m.def("slice", []
(
  const Eigen::SparseMatrix<int>& X,
  const Eigen::MatrixXi& R,
  const int& dim,
  Eigen::SparseMatrix<int>& Y
)
{
  assert_is_VectorX("R",R);
  return igl::slice(X,R,dim,Y);
}, __doc_igl_slice,
py::arg("X"), py::arg("R"), py::arg("dim"), py::arg("Y"));

m.def("slice", []
(
  const Eigen::MatrixXi& X,
  const Eigen::MatrixXi& R,
  const Eigen::MatrixXi& C,
  Eigen::MatrixXi& Y
)
{
  assert_is_VectorX("R",R);
  assert_is_VectorX("C",C);
  return igl::slice(X,R,C,Y);
}, __doc_igl_slice,
py::arg("X"), py::arg("R"), py::arg("C"), py::arg("Y"));


m.def("slice", []
(
  const Eigen::MatrixXi& X,
  const Eigen::MatrixXi& R,
  Eigen::MatrixXi& Y
)
{
  assert_is_VectorX("R",R);
  return igl::slice(X,R,Y);
}, __doc_igl_slice,
py::arg("X"), py::arg("R"), py::arg("Y"));

m.def("slice", []
(
  const Eigen::MatrixXi& X,
  const Eigen::MatrixXi& R
)
{
  assert_is_VectorX("R",R);
  return igl::slice(X,R);
}, __doc_igl_slice,
py::arg("X"), py::arg("R"));

m.def("slice", []
(
  const Eigen::MatrixXi& X,
  const Eigen::MatrixXi& R,
  const int& dim
)
{
  assert_is_VectorX("R",R);
  return igl::slice(X,R,dim);
}, __doc_igl_slice,
py::arg("X"), py::arg("R"), py::arg("dim"));
