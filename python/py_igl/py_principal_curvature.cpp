// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("principal_curvature", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXd& PD1,
  Eigen::MatrixXd& PD2,
  Eigen::MatrixXd& PV1,
  Eigen::MatrixXd& PV2,
  unsigned radius,
  bool useKring
)
{
  return igl::principal_curvature(V,F,PD1,PD2,PV1,PV2,radius,useKring);
}, __doc_igl_principal_curvature,
py::arg("V"), py::arg("F"), py::arg("PD1"), py::arg("PD2"), py::arg("PV1"), py::arg("PV2"), py::arg("radius") = 5, py::arg("useKring") = true);
