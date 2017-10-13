// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("doublearea", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXd& dblA
)
{
  return igl::doublearea(V,F,dblA);
}, __doc_igl_doublearea,
py::arg("V"), py::arg("F"), py::arg("dblA"));

m.def("doublearea", []
(
  const Eigen::MatrixXd& A,
  const Eigen::MatrixXd& B,
  const Eigen::MatrixXd& C,
  Eigen::MatrixXd& D
)
{
  return igl::doublearea(A,B,C,D);
}, __doc_igl_doublearea,
py::arg("A"), py::arg("B"), py::arg("C"), py::arg("D"));

m.def("doublearea_single", []
(
  const Eigen::MatrixXd& A,
  const Eigen::MatrixXd& B,
  const Eigen::MatrixXd& C
)
{
  return igl::doublearea_single(A,B,C);
}, __doc_igl_doublearea_single,
py::arg("A"), py::arg("B"), py::arg("C"));

m.def("doublearea", []
(
  const Eigen::MatrixXd& l,
  Eigen::MatrixXd& dblA
)
{
  return igl::doublearea(l,dblA);
}, __doc_igl_doublearea,
py::arg("l"), py::arg("dblA"));

m.def("doublearea_quad", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXd& dblA
)
{
  return igl::doublearea_quad(V,F,dblA);
}, __doc_igl_doublearea_quad,
py::arg("V"), py::arg("F"), py::arg("dblA"));
