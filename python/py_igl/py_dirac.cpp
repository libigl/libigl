// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Zhongshi Jiang <jiangzs@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


m.def("dirac", []
(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    Eigen::SparseMatrix<double> &D,
    Eigen::SparseMatrix<double> &DA
)
{
  return igl::dirac(V, F, D, DA);
}, __doc_igl_dirac,
py::arg("V"), py::arg("F"), py::arg("D"), py::arg("DA"));

