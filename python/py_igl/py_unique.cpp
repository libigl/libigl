// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("unique", []
(
  const Eigen::MatrixXd& A,
  Eigen::MatrixXd& C,
  Eigen::MatrixXi& IA,
  Eigen::MatrixXi& IC
)
{
  return igl::unique(A,C,IA,IC);
}, __doc_igl_unique,
py::arg("A"), py::arg("C"), py::arg("IA"), py::arg("IC"));

m.def("unique", []
(
  const Eigen::MatrixXd& A,
  Eigen::MatrixXd& C
)
{
  return igl::unique(A,C);
}, __doc_igl_unique,
py::arg("A"), py::arg("C"));

//m.def("unique", []
//(
//  const std::vector<double> & A,
//  std::vector<double> & C,
//  std::vector<size_t> & IA,
//  std::vector<size_t> & IC
//)
//{
//  return igl::unique(A,C,IA,IC);
//}, __doc_igl_unique,
//py::arg("A"), py::arg("C"), py::arg("IA"), py::arg("IC"));

//m.def("unique", []
//(
//  const std::vector<double> & A,
//  std::vector<double> & C
//)
//{
//  return igl::unique(A,C);
//}, __doc_igl_unique,
//py::arg("A"), py::arg("C"));


// int


m.def("unique", []
(
  const Eigen::MatrixXi& A,
  Eigen::MatrixXi& C,
  Eigen::MatrixXi& IA,
  Eigen::MatrixXi& IC
)
{
  return igl::unique(A,C,IA,IC);
}, __doc_igl_unique,
py::arg("A"), py::arg("C"), py::arg("IA"), py::arg("IC"));

m.def("unique", []
(
  const Eigen::MatrixXi& A,
  Eigen::MatrixXi& C
)
{
  return igl::unique(A,C);
}, __doc_igl_unique,
py::arg("A"), py::arg("C"));

//m.def("unique", []
//(
//  const std::vector<int> & A,
//  std::vector<int> & C,
//  std::vector<size_t> & IA,
//  std::vector<size_t> & IC
//)
//{
//  return igl::unique(A,C,IA,IC);
//}, __doc_igl_unique,
//py::arg("A"), py::arg("C"), py::arg("IA"), py::arg("IC"));

//m.def("unique", []
//(
//  const std::vector<int> & A,
//  std::vector<int> & C
//)
//{
//  return igl::unique(A,C);
//}, __doc_igl_unique,
//py::arg("A"), py::arg("C"));
