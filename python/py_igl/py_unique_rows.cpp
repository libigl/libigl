// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("unique_rows", []
(
  const Eigen::MatrixXd& A,
  Eigen::MatrixXd& C,
  Eigen::MatrixXi& IA,
  Eigen::MatrixXi& IC
)
{
  return igl::unique_rows(A,C,IA,IC);
}, __doc_igl_unique,
py::arg("A"), py::arg("C"), py::arg("IA"), py::arg("IC"));

m.def("unique_rows", []
(
  const Eigen::MatrixXi& A,
  Eigen::MatrixXi& C,
  Eigen::MatrixXi& IA,
  Eigen::MatrixXi& IC
)
{
  return igl::unique_rows(A,C,IA,IC);
}, __doc_igl_unique,
py::arg("A"), py::arg("C"), py::arg("IA"), py::arg("IC"));

