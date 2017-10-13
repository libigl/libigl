// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("jet", []
(
  const Eigen::MatrixXd& Z,
  const bool normalize,
  Eigen::MatrixXd& C
)
{
  assert_is_VectorX("Z",Z);
  return igl::jet(Z,normalize,C);
}, __doc_igl_jet,
py::arg("Z"), py::arg("normalize"), py::arg("C"));

m.def("jet", []
(
  const Eigen::MatrixXd& Z,
  const double min_Z,
  const double max_Z,
  Eigen::MatrixXd& C
)
{
  assert_is_VectorX("Z",Z);
  return igl::jet(Z,min_Z,max_Z,C);
}, __doc_igl_jet,
py::arg("Z"), py::arg("min_Z"), py::arg("max_Z"), py::arg("C"));
