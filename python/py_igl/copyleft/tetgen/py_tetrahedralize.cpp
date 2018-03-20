// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

m.def("tetrahedralize", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const std::string switches,
  Eigen::MatrixXd& TV,
  Eigen::MatrixXi& TT,
  Eigen::MatrixXi& TF
)
{
  return igl::copyleft::tetgen::tetrahedralize(V, F, switches, TV, TT, TF);
}, __doc_igl_copyleft_tetgen_tetrahedralize,
py::arg("V"), py::arg("F"), py::arg("switches"), py::arg("TV"), py::arg("TT"), py::arg("TF"));

m.def("tetrahedralize", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXi& VM,
  const Eigen::MatrixXi& FM,
  const std::string switches,
  Eigen::MatrixXd& TV,
  Eigen::MatrixXi& TT,
  Eigen::MatrixXi& TF,
  Eigen::MatrixXi& TM
)
{
  return igl::copyleft::tetgen::tetrahedralize(V, F, VM, FM, switches, TV, TT, TF, TM);
}, __doc_igl_copyleft_tetgen_tetrahedralize,
py::arg("V"), py::arg("F"), py::arg("VM"), py::arg("FM"), py::arg("switches"), py::arg("TV"), py::arg("TT"), py::arg("TF"), py::arg("TM"));

