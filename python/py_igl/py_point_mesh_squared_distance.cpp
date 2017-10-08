// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("point_mesh_squared_distance", []
(
  const Eigen::MatrixXd& P,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& Ele,
  Eigen::MatrixXd& sqrD,
  Eigen::MatrixXi& I,
  Eigen::MatrixXd& C
)
{
//  assert_is_VectorX("I",I);
  return igl::point_mesh_squared_distance(P, V, Ele, sqrD, I, C);
}, __doc_igl_point_mesh_squared_distance,
py::arg("P"), py::arg("V"), py::arg("Ele"), py::arg("sqrD"), py::arg("I"), py::arg("C"));

