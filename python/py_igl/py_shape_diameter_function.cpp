// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

m.def("shape_diameter_function", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& P,
  const Eigen::MatrixXd& N,
  const int num_samples,
  Eigen::MatrixXd& S
)
{
  return igl::shape_diameter_function(V, F, P, N, num_samples, S);
}, __doc_igl_shape_diameter_function,
py::arg("V"), py::arg("F"), py::arg("P"), py::arg("N"), py::arg("num_samples"), py::arg("S"));

m.def("shape_diameter_function", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const bool per_face,
  const int num_samples,
  Eigen::MatrixXd& S
)
{
  return igl::shape_diameter_function(V, F, per_face, num_samples, S);
}, __doc_igl_shape_diameter_function,
py::arg("V"), py::arg("F"), py::arg("per_face"), py::arg("num_samples"), py::arg("S"));

