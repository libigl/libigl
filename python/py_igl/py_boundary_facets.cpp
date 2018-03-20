// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("boundary_facets", []
(
  const Eigen::MatrixXi& T,
  Eigen::MatrixXi& F
)
{
  return igl::boundary_facets(T,F);
}, __doc_igl_boundary_facets,
py::arg("T"), py::arg("F"));

m.def("boundary_facets", []
(
  const Eigen::MatrixXi& T
)
{
  Eigen::MatrixXi F;
  igl::boundary_facets(T,F);
  return F;
}, __doc_igl_boundary_facets,
py::arg("T"));

m.def("boundary_facets", []
(
  const std::vector<std::vector<int> > & T,
  std::vector<std::vector<int> > & F
)
{
  return igl::boundary_facets(T,F);
}, __doc_igl_boundary_facets,
py::arg("T"), py::arg("F"));
