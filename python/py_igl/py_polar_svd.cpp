// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("polar_svd", []
(
  const Eigen::MatrixXd& A,
  Eigen::MatrixXd& R,
  Eigen::MatrixXd& T,
  Eigen::MatrixXd& U,
  Eigen::MatrixXd& S,
  Eigen::MatrixXd& V
)
{
  return igl::polar_svd(A, R, T, U, S, V);
}, __doc_igl_polar_svd,
py::arg("A"), py::arg("R"), py::arg("T"), py::arg("U"), py::arg("S"), py::arg("V"));


m.def("polar_svd", []
(
  const Eigen::MatrixXd& A,
  Eigen::MatrixXd& R,
  Eigen::MatrixXd& T
)
{
  return igl::polar_svd(A, R, T);
}, __doc_igl_polar_svd,
py::arg("A"), py::arg("R"), py::arg("T"));

