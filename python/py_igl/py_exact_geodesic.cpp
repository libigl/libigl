// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Zhongshi Jiang <jiangzs@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


m.def("exact_geodesic", []
(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    const Eigen::MatrixXi &VS,
    const Eigen::MatrixXi &FS,
    const Eigen::MatrixXi &VT,
    const Eigen::MatrixXi &FT,
    Eigen::MatrixXd &D
)
{
  return igl::exact_geodesic(V, F, VS,FS,VT,FT, D);
}, __doc_igl_exact_geodesic,
py::arg("V"), py::arg("F"), py::arg("VS"), py::arg("FS"), py::arg("VT"), py::arg("FT"), py::arg("D"));

