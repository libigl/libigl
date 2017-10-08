// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

m.def("planarize_quad_mesh", []
(
  const Eigen::MatrixXd& Vin,
  const Eigen::MatrixXi& F,
  const int maxIter,
  const double & threshold,
  Eigen::MatrixXd& Vout
)
{
  return igl::planarize_quad_mesh(Vin, F, maxIter, threshold, Vout);
}, __doc_igl_planarize_quad_mesh,
py::arg("Vin"), py::arg("F"), py::arg("maxIter"), py::arg("threshold"), py::arg("Vout"));

