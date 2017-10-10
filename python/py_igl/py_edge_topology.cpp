// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("edge_topology", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXi& EV,
  Eigen::MatrixXi& FE,
  Eigen::MatrixXi& EF
)
{
  return igl::edge_topology(V, F, EV, FE, EF);
}, __doc_igl_edge_lengths,
py::arg("V"), py::arg("F"), py::arg("EV"), py::arg("FE"), py::arg("EF"));

