// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
// COMPLETE BINDINGS ========================


m.def("lbs_matrix", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXd& W,
  Eigen::MatrixXd& M
)
{
  return igl::lbs_matrix(V, W, M);
}, __doc_igl_lbs_matrix,
py::arg("V"), py::arg("W"), py::arg("M"));

m.def("lbs_matrix_column", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXd& W,
  Eigen::MatrixXd& M
)
{
  return igl::lbs_matrix_column(V, W, M);
}, __doc_igl_lbs_matrix_column,
py::arg("V"), py::arg("W"), py::arg("M"));

m.def("lbs_matrix_column", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXd& W,
  const Eigen::MatrixXi& WI,
  Eigen::MatrixXd& M
)
{
  return igl::lbs_matrix_column(V, W, WI, M);
}, __doc_igl_lbs_matrix_column,
py::arg("V"), py::arg("W"), py::arg("WI"), py::arg("M"));





// INCOMPLETE BINDINGS ========================


m.def("lbs_matrix_column", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXd& W,
  Eigen::SparseMatrix<double>& M
)
{
  return igl::lbs_matrix_column(V, W, M);
}, __doc_igl_lbs_matrix_column,
py::arg("V"), py::arg("W"), py::arg("M"));

m.def("lbs_matrix_column", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXd& W,
  const Eigen::MatrixXi& WI,
  Eigen::SparseMatrix<double>& M
)
{
  return igl::lbs_matrix_column(V, W, WI, M);
}, __doc_igl_lbs_matrix_column,
py::arg("V"), py::arg("W"), py::arg("WI"), py::arg("M"));

