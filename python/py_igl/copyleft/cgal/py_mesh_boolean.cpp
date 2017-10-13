// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
// COMPLETE BINDINGS ========================

m.def("mesh_boolean", []
(
  const Eigen::MatrixXd& VA,
  const Eigen::MatrixXi& FA,
  const Eigen::MatrixXd& VB,
  const Eigen::MatrixXi& FB,
  igl::MeshBooleanType & type,
  Eigen::MatrixXd& VC,
  Eigen::MatrixXi& FC,
  Eigen::MatrixXi& J
)
{
  return igl::copyleft::cgal::mesh_boolean(VA, FA, VB, FB, type, VC, FC, J);
}, __doc_igl_copyleft_cgal_mesh_boolean,
py::arg("VA"), py::arg("FA"), py::arg("VB"), py::arg("FB"), py::arg("type"), py::arg("VC"), py::arg("FC"), py::arg("J"));


m.def("mesh_boolean", []
(
  const Eigen::MatrixXd& VA,
  const Eigen::MatrixXi& FA,
  const Eigen::MatrixXd& VB,
  const Eigen::MatrixXi& FB,
  const std::string & type_str,
  Eigen::MatrixXd& VC,
  Eigen::MatrixXi& FC,
  Eigen::MatrixXi& J
)
{
  return igl::copyleft::cgal::mesh_boolean(VA, FA, VB, FB, type_str, VC, FC, J);
}, __doc_igl_copyleft_cgal_mesh_boolean,
py::arg("VA"), py::arg("FA"), py::arg("VB"), py::arg("FB"), py::arg("type_str"), py::arg("VC"), py::arg("FC"), py::arg("J"));

m.def("mesh_boolean", []
(
  const Eigen::MatrixXd& VA,
  const Eigen::MatrixXi& FA,
  const Eigen::MatrixXd& VB,
  const Eigen::MatrixXi& FB,
  const igl::MeshBooleanType & type,
  Eigen::MatrixXd& VC,
  Eigen::MatrixXi& FC
)
{
  return igl::copyleft::cgal::mesh_boolean(VA, FA, VB, FB, type, VC, FC);
}, __doc_igl_copyleft_cgal_mesh_boolean,
py::arg("VA"), py::arg("FA"), py::arg("VB"), py::arg("FB"), py::arg("type"), py::arg("VC"), py::arg("FC"));



// INCOMPLETE BINDINGS ========================




//m.def("mesh_boolean", []
//(
//  const Eigen::MatrixXd& VA,
//  const Eigen::MatrixXd& FA,
//  const Eigen::MatrixXd& VB,
//  const Eigen::MatrixXd& FB,
//  std::function<int (const Eigen::Matrix<int, 1, Eigen::Dynamic>)> & wind_num_op,
//  std::function<int (const int, const int)> & keep,
//  Eigen::MatrixXd& VC,
//  Eigen::MatrixXd& FC,
//  Eigen::MatrixXd& J
//)
//{
//  return igl::copyleft::cgal::mesh_boolean(VA, FA, VB, FB, wind_num_op, keep, VC, FC, J);
//}, __doc_igl_copyleft_cgal_mesh_boolean,
//py::arg("VA"), py::arg("FA"), py::arg("VB"), py::arg("FB"), py::arg("wind_num_op"), py::arg("keep"), py::arg("VC"), py::arg("FC"), py::arg("J"));

//m.def("mesh_boolean", []
//(
//  std::vector<DerivedV> & Vlist,
//  std::vector<DerivedF> & Flist,
//  std::function<int (const Eigen::Matrix<int, 1, Eigen::Dynamic>)> & wind_num_op,
//  std::function<int (const int, const int)> & keep,
//  Eigen::MatrixXd& VC,
//  Eigen::MatrixXd& FC,
//  Eigen::MatrixXd& J
//)
//{
//  return igl::copyleft::cgal::mesh_boolean(Vlist, Flist, wind_num_op, keep, VC, FC, J);
//}, __doc_igl_copyleft_cgal_mesh_boolean,
//py::arg("Vlist"), py::arg("Flist"), py::arg("wind_num_op"), py::arg("keep"), py::arg("VC"), py::arg("FC"), py::arg("J"));

//m.def("mesh_boolean", []
//(
//  const Eigen::MatrixXd& VV,
//  const Eigen::MatrixXd& FF,
//  const Eigen::MatrixXd& sizes,
//  std::function<int (const Eigen::Matrix<int, 1, Eigen::Dynamic>)> & wind_num_op,
//  std::function<int (const int, const int)> & keep,
//  Eigen::MatrixXd& VC,
//  Eigen::MatrixXd& FC,
//  Eigen::MatrixXd& J
//)
//{
//  return igl::copyleft::cgal::mesh_boolean(VV, FF, sizes, wind_num_op, keep, VC, FC, J);
//}, __doc_igl_copyleft_cgal_mesh_boolean,
//py::arg("VV"), py::arg("FF"), py::arg("sizes"), py::arg("wind_num_op"), py::arg("keep"), py::arg("VC"), py::arg("FC"), py::arg("J"));



