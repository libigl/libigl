// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

m.def("triangle_triangle_adjacency", []
(
  const Eigen::MatrixXi& F,
  Eigen::MatrixXi& TT,
  Eigen::MatrixXi& TTi
)
{
  return igl::triangle_triangle_adjacency(F, TT, TTi);
}, __doc_igl_triangle_triangle_adjacency,
py::arg("F"), py::arg("TT"), py::arg("TTi"));

m.def("triangle_triangle_adjacency", []
(
  const Eigen::MatrixXi& F,
  Eigen::MatrixXi& TT
)
{
  return igl::triangle_triangle_adjacency(F, TT);
}, __doc_igl_triangle_triangle_adjacency,
py::arg("F"), py::arg("TT"));
