// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("harmonic", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXi& b,
  const Eigen::MatrixXd& bc,
  const int k,
  Eigen::MatrixXd& W
)
{
  assert_is_VectorX("b",b);
  return igl::harmonic(V,F,b,bc,k,W);
}, __doc_igl_harmonic,
py::arg("V"), py::arg("F"), py::arg("b"), py::arg("bc"), py::arg("k"), py::arg("W"));
