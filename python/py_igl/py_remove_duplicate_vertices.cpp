// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("remove_duplicate_vertices", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const double epsilon,
  Eigen::MatrixXd& SV,
  Eigen::MatrixXi& SVI,
  Eigen::MatrixXi& SVJ,
  Eigen::MatrixXi& SF
)
{
  return igl::remove_duplicate_vertices(V, F, epsilon, SV, SVI, SVJ, SF);
}, __doc_igl_remove_duplicate_vertices,
py::arg("V"), py::arg("F"), py::arg("epsilon"), py::arg("SV"), py::arg("SVI"), py::arg("SVJ"), py::arg("SF"));
