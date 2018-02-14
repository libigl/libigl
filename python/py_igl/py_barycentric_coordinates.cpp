// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


m.def("barycentric_coordinates", []
(
  const Eigen::MatrixXd& P,
  const Eigen::MatrixXd& A,
  const Eigen::MatrixXd& B,
  const Eigen::MatrixXd& C,
  const Eigen::MatrixXd& D,
  Eigen::MatrixXd& L
)
{
  return igl::barycentric_coordinates(P, A, B, C, D, L);
}, __doc_igl_barycentric_coordinates,
py::arg("P"), py::arg("A"), py::arg("B"), py::arg("C"), py::arg("D"), py::arg("L"));

m.def("barycentric_coordinates", []
(
  const Eigen::MatrixXd& P,
  const Eigen::MatrixXd& A,
  const Eigen::MatrixXd& B,
  const Eigen::MatrixXd& C,
  Eigen::MatrixXd& L
)
{
  return igl::barycentric_coordinates(P, A, B, C, L);
}, __doc_igl_barycentric_coordinates,
py::arg("P"), py::arg("A"), py::arg("B"), py::arg("C"), py::arg("L"));

