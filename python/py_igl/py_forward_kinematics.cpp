// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

//m.def("forward_kinematics", []
//(
//  const Eigen::MatrixXd& C,
//  const Eigen::MatrixXi& BE,
//  const Eigen::MatrixXi& P,
//  std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond> > & dQ,
//  std::vector<Eigen::Vector3d> & dT,
//  std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond> > & vQ,
//  std::vector<Eigen::Vector3d> & vT
//)
//{
//  return igl::forward_kinematics(C, BE, P, dQ, dT, vQ, vT);
//}, __doc_igl_forward_kinematics,
//py::arg("C"), py::arg("BE"), py::arg("P"), py::arg("dQ"), py::arg("dT"), py::arg("vQ"), py::arg("vT"));

m.def("forward_kinematics", []
(
  const Eigen::MatrixXd& C,
  const Eigen::MatrixXi& BE,
  const Eigen::MatrixXi& P,
  const RotationList& dQ,
  RotationList& vQ,
  py::list vT
)
{
  std::vector<Eigen::Vector3d> vTl;
  igl::forward_kinematics(C, BE, P, dQ, vQ, vTl);
  for (auto item : vTl) {
    py::object obj = py::cast(Eigen::MatrixXd(item));
    vT.append(obj);
  }
}, __doc_igl_forward_kinematics,
py::arg("C"), py::arg("BE"), py::arg("P"), py::arg("dQ"), py::arg("vQ"), py::arg("vT"));

//m.def("forward_kinematics", []
//(
//  const Eigen::MatrixXd& C,
//  const Eigen::MatrixXi& BE,
//  const Eigen::MatrixXi& P,
//  std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond> > & dQ,
//  std::vector<Eigen::Vector3d> & dT,
//  Eigen::MatrixXd& T
//)
//{
//  return igl::forward_kinematics(C, BE, P, dQ, dT, T);
//}, __doc_igl_forward_kinematics,
//py::arg("C"), py::arg("BE"), py::arg("P"), py::arg("dQ"), py::arg("dT"), py::arg("T"));

//m.def("forward_kinematics", []
//(
//  const Eigen::MatrixXd& C,
//  const Eigen::MatrixXi& BE,
//  const Eigen::MatrixXi& P,
//  std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond> > & dQ,
//  Eigen::MatrixXd& T
//)
//{
//  return igl::forward_kinematics(C, BE, P, dQ, T);
//}, __doc_igl_forward_kinematics,
//py::arg("C"), py::arg("BE"), py::arg("P"), py::arg("dQ"), py::arg("T"));

