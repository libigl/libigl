// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


m.def("boundary_conditions", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& Ele,
  const Eigen::MatrixXd& C,
  const Eigen::MatrixXi& P,
  const Eigen::MatrixXi& BE,
  const Eigen::MatrixXi& CE,
  Eigen::MatrixXi& b,
  Eigen::MatrixXd& bc
)
{
  assert_is_VectorX("P", P);
  Eigen::VectorXi Pv;
  if (P.size() != 0)
    Pv = P;
  Eigen::VectorXi bv;
  igl::boundary_conditions(V, Ele, C, Pv, BE, CE, bv, bc);
  b = bv;
}, __doc_igl_boundary_conditions,
py::arg("V"), py::arg("Ele"), py::arg("C"), py::arg("P"), py::arg("BE"), py::arg("CE"), py::arg("b"), py::arg("bc"));

