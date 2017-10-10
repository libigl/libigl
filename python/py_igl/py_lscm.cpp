// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("lscm", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXi& b,
  const Eigen::MatrixXd& bc,
  Eigen::MatrixXd& V_uv
)
{
  assert_is_VectorX("b",b);
  return igl::lscm(V,F,b,bc,V_uv);
}, __doc_igl_lscm,
py::arg("V"), py::arg("F"), py::arg("b"), py::arg("bc"), py::arg("V_uv"));
