// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("marching_tets", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& T,
  const Eigen::MatrixXd& plane,
  Eigen::MatrixXd& U,
  Eigen::MatrixXi& G,
  Eigen::MatrixXi& J,
  Eigen::SparseMatrix<double>& BC
)
{
  assert_is_VectorX("plane", plane);
  Eigen::VectorXd planev;
  if (plane.size() != 0)
    planev = plane;
  Eigen::VectorXi Jv;
  igl::marching_tets(V, T, planev, U, G, Jv, BC);
  J = Jv;
}, __doc_igl_marching_tets,
py::arg("V"), py::arg("T"), py::arg("plane"), py::arg("U"), py::arg("G"), py::arg("J"), py::arg("BC"));

