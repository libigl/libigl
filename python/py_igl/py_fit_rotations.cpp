// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("fit_rotations", []
(
  const Eigen::MatrixXd& S,
  const bool single_precision,
  Eigen::MatrixXd& R
)
{
  return igl::fit_rotations(S, single_precision, R);
}, __doc_igl_fit_rotations,
py::arg("S"), py::arg("single_precision"), py::arg("R"));


m.def("fit_rotations_planar", []
(
  const Eigen::MatrixXd& S,
  Eigen::MatrixXd& R
)
{
  return igl::fit_rotations_planar(S, R);
}, __doc_igl_fit_rotations_planar,
py::arg("S"), py::arg("R"));

