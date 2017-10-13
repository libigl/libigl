// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

m.def("directed_edge_orientations", []
(
  const Eigen::MatrixXd& C,
  const Eigen::MatrixXi& E,
  RotationList& Q
)
{
  return igl::directed_edge_orientations(C, E, Q);
}, __doc_igl_directed_edge_orientations,
py::arg("C"), py::arg("E"), py::arg("Q"));

