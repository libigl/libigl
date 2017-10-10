// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
py::enum_<igl::MassMatrixType>(m, "MassMatrixType")
    .value("MASSMATRIX_TYPE_BARYCENTRIC", igl::MASSMATRIX_TYPE_BARYCENTRIC)
    .value("MASSMATRIX_TYPE_VORONOI", igl::MASSMATRIX_TYPE_VORONOI)
    .value("MASSMATRIX_TYPE_FULL", igl::MASSMATRIX_TYPE_FULL)
    .value("MASSMATRIX_TYPE_DEFAULT", igl::MASSMATRIX_TYPE_DEFAULT)
    .value("NUM_MASSMATRIX_TYPE", igl::NUM_MASSMATRIX_TYPE)
    .export_values();

m.def("massmatrix", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const igl::MassMatrixType type,
  Eigen::SparseMatrix<double>& M
)
{
  return igl::massmatrix(V,F,type,M);
}, __doc_igl_massmatrix,
py::arg("V"), py::arg("F"), py::arg("type"), py::arg("M"));
