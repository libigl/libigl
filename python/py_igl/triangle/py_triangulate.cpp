// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


m.def("triangulate", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& E,
  const Eigen::MatrixXd& H,
  const std::string flags,
  Eigen::MatrixXd& V2,
  Eigen::MatrixXi& F2
)
{
  return igl::triangle::triangulate(V, E, H, flags, V2, F2);
}, __doc_igl_triangle_triangulate,
py::arg("V"), py::arg("E"), py::arg("H"), py::arg("flags"), py::arg("V2"), py::arg("F2"));

