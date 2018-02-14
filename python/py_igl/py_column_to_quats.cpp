// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("column_to_quats", []
(
  const Eigen::MatrixXd& Q,
  RotationList& vQ
)
{
  return igl::column_to_quats(Q, vQ);
}, __doc_igl_column_to_quats,
py::arg("Q"), py::arg("vQ"));

