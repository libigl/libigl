// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("comb_cross_field", []
(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  const Eigen::MatrixXd &PD1in,
  const Eigen::MatrixXd &PD2in,
  Eigen::MatrixXd &PD1out,
  Eigen::MatrixXd &PD2out
)
{
  return igl::comb_cross_field(V,F,PD1in,PD2in,PD1out,PD2out);
}, __doc_igl_comb_cross_field,
py::arg("V"), py::arg("F"), py::arg("PD1in"), py::arg("PD2in"), py::arg("PD1out"), py::arg("PD2out"));
