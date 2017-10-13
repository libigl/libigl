// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("boundary_loop", []
(
  const Eigen::MatrixXi& F,
  Eigen::MatrixXi& L
)
{
  Eigen::VectorXi T;
  igl::boundary_loop(F,T);
  L = T;
}, __doc_igl_boundary_loop,
py::arg("F"), py::arg("L"));

m.def("boundary_loop", []
(
  const Eigen::MatrixXi& F,
  std::vector<std::vector<int> >& L
)
{
  return igl::boundary_loop(F,L);
}, __doc_igl_boundary_loop,
py::arg("F"), py::arg("L"));

m.def("boundary_loop", []
(
  const Eigen::MatrixXi& F,
  std::vector<int>& L
)
{
  return igl::boundary_loop(F,L);
}, __doc_igl_boundary_loop,
py::arg("F"), py::arg("L"));


