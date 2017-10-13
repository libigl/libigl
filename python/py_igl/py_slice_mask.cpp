// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("slice_mask", []
(
  const Eigen::MatrixXd& X,
  const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> & R,
  const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> & C,
  Eigen::MatrixXd& Y
)
{
  assert_is_VectorX("R",R);
  assert_is_VectorX("C",C);
  return igl::slice_mask(X, R, C, Y);
}, __doc_igl_slice_mask,
py::arg("X"), py::arg("R"), py::arg("C"), py::arg("Y"));

m.def("slice_mask", []
(
  const Eigen::MatrixXd& X,
  const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> & R,
  const int dim,
  Eigen::MatrixXd& Y
)
{
  assert_is_VectorX("R",R);
  return igl::slice_mask(X, R, dim, Y);
}, __doc_igl_slice_mask,
py::arg("X"), py::arg("R"), py::arg("dim"), py::arg("Y"));

m.def("slice_mask", []
(
  const Eigen::MatrixXi& X,
  const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> & R,
  const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> & C,
  Eigen::MatrixXi& Y
)
{
  assert_is_VectorX("R",R);
  assert_is_VectorX("C",C);
  return igl::slice_mask(X, R, C, Y);
}, __doc_igl_slice_mask,
py::arg("X"), py::arg("R"), py::arg("C"), py::arg("Y"));

m.def("slice_mask", []
(
  const Eigen::MatrixXi& X,
  const Eigen::Matrix<bool, Eigen::Dynamic, Eigen::Dynamic> & R,
  const int dim,
  Eigen::MatrixXi& Y
)
{
  assert_is_VectorX("R",R);
  return igl::slice_mask(X, R, dim, Y);
}, __doc_igl_slice_mask,
py::arg("X"), py::arg("R"), py::arg("dim"), py::arg("Y"));

//m.def("slice_mask", []
//(
//  const Eigen::MatrixXd& X,
//  Eigen::Array<bool, Eigen::Dynamic, 1> & R,
//  Eigen::Array<bool, Eigen::Dynamic, 1> & C
//)
//{
//  return igl::slice_mask(X, R, C);
//}, __doc_igl_slice_mask,
//py::arg("X"), py::arg("R"), py::arg("C"));

//m.def("slice_mask", []
//(
//  const Eigen::MatrixXd& X,
//  Eigen::Array<bool, Eigen::Dynamic, 1> & R,
//  int dim
//)
//{
//  return igl::slice_mask(X, R, dim);
//}, __doc_igl_slice_mask,
//py::arg("X"), py::arg("R"), py::arg("dim"));

