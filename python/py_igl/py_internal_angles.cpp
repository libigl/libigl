// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


m.def("internal_angles", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXd& K
)
{
  return igl::internal_angles(V, F, K);
}, __doc_igl_internal_angles,
py::arg("V"), py::arg("F"), py::arg("K"));

m.def("internal_angles_using_squared_edge_lengths", []
(
  const Eigen::MatrixXd& L_sq,
  Eigen::MatrixXd& K
)
{
  return igl::internal_angles_using_squared_edge_lengths(L_sq, K);
}, __doc_igl_internal_angles,
py::arg("L_sq"), py::arg("K"));

//m.def("internal_angles_using_edge_lengths", []
//(
//  const Eigen::MatrixXd& L,
//  Eigen::MatrixXd& K
//)
//{
//  return igl::internal_angles_using_edge_lengths(L, K);
//}, __doc_igl_internal_angles,
//py::arg("L"), py::arg("K"));


