// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("local_basis", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXd& B1,
  Eigen::MatrixXd& B2,
  Eigen::MatrixXd& B3
)
{
  return igl::local_basis(V,F,B1,B2,B3);
}, __doc_igl_local_basis,
py::arg("V"), py::arg("F"), py::arg("B1"), py::arg("B2"), py::arg("B3"));
