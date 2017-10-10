// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
// COMPLETE BINDINGS ========================


m.def("winding_number", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& O,
  Eigen::MatrixXd& W
)
{
  Eigen::VectorXd Wv;
  igl::winding_number(V, F, O, Wv);
  W = Wv;
}, __doc_igl_winding_number,
py::arg("V"), py::arg("F"), py::arg("O"), py::arg("W"));


