// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

py::enum_<igl::SignedDistanceType>(m, "SignedDistanceType")
    .value("SIGNED_DISTANCE_TYPE_PSEUDONORMAL", igl::SIGNED_DISTANCE_TYPE_PSEUDONORMAL)
    .value("SIGNED_DISTANCE_TYPE_WINDING_NUMBER", igl::SIGNED_DISTANCE_TYPE_WINDING_NUMBER)
    .value("SIGNED_DISTANCE_TYPE_DEFAULT", igl::SIGNED_DISTANCE_TYPE_DEFAULT)
    .value("SIGNED_DISTANCE_TYPE_UNSIGNED", igl::SIGNED_DISTANCE_TYPE_UNSIGNED)
    .value("NUM_SIGNED_DISTANCE_TYPE", igl::NUM_SIGNED_DISTANCE_TYPE)
    .export_values();


m.def("signed_distance", []
(
  const Eigen::MatrixXd& P,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const igl::SignedDistanceType sign_type,
  Eigen::MatrixXd& S,
  Eigen::MatrixXi& I,
  Eigen::MatrixXd& C,
  Eigen::MatrixXd& N
)
{
  Eigen::VectorXd Sv;
  Eigen::VectorXi Iv;
  igl::signed_distance(P, V, F, sign_type, Sv, Iv, C, N);
  S = Sv;
  I = Iv;
}, __doc_igl_signed_distance,
py::arg("P"), py::arg("V"), py::arg("F"), py::arg("sign_type"), py::arg("S"), py::arg("I"), py::arg("C"), py::arg("N"));

//m.def("signed_distance_pseudonormal", []
//(
//  const AABB<Eigen::MatrixXd, 3> & tree,
//  const Eigen::MatrixXd& V,
//  const Eigen::MatrixXi& F,
//  const Eigen::MatrixXd& FN,
//  const Eigen::MatrixXd& VN,
//  const Eigen::MatrixXd& EN,
//  const Eigen::MatrixXi& EMAP,
//  const Eigen::MatrixXd& q
//)
//{
//  assert_is_VectorX("q", q);
//  assert_is_VectorX("EMAP",EMAP);
//  return igl::signed_distance_pseudonormal(tree, V, F, FN, VN, EN, EMAP, q);
//}, __doc_igl_signed_distance_pseudonormal,
//py::arg("tree"), py::arg("V"), py::arg("F"), py::arg("FN"), py::arg("VN"), py::arg("EN"), py::arg("EMAP"), py::arg("q"));

m.def("signed_distance_pseudonormal", []
(
  const Eigen::MatrixXd& P,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const igl::AABB<Eigen::MatrixXd, 3> & tree,
  const Eigen::MatrixXd& FN,
  const Eigen::MatrixXd& VN,
  const Eigen::MatrixXd& EN,
  const Eigen::MatrixXi& EMAP,
  Eigen::MatrixXd& S,
  Eigen::MatrixXi& I,
  Eigen::MatrixXd& C,
  Eigen::MatrixXd& N
)
{
  assert_is_VectorX("EMAP", EMAP);
  Eigen::VectorXi EMAPv;
  if (EMAP.size() != 0)
    EMAPv = EMAP;
  Eigen::VectorXd Sv;
  Eigen::VectorXi Iv;
  igl::signed_distance_pseudonormal(P, V, F, tree, FN, VN, EN, EMAPv, Sv, Iv, C, N);
  S = Sv;
  I = Iv;
}, __doc_igl_signed_distance_pseudonormal,
py::arg("P"), py::arg("V"), py::arg("F"), py::arg("tree"), py::arg("FN"), py::arg("VN"), py::arg("EN"), py::arg("EMAP"), py::arg("S"), py::arg("I"), py::arg("C"), py::arg("N"));

//m.def("signed_distance_pseudonormal", []
//(
//  const AABB<Eigen::MatrixXd, 3> & tree,
//  const Eigen::MatrixXd& V,
//  const Eigen::MatrixXi& F,
//  const Eigen::MatrixXd& FN,
//  const Eigen::MatrixXd& VN,
//  const Eigen::MatrixXd& EN,
//  const Eigen::MatrixXi & EMAP,
//  const Eigen::MatrixXd & q,
//  double & s,
//  double & sqrd,
//  int & i,
//  Eigen::MatrixXd & c,
//  Eigen::MatrixXd & n
//)
//{
//  assert_is_VectorX("EMAP",EMAP);
//  assert_is_VectorX("q",q);
//  return igl::signed_distance_pseudonormal(tree, V, F, FN, VN, EN, EMAP, q, s, sqrd, i, c, n);
//}, __doc_igl_signed_distance_pseudonormal,
//py::arg("tree"), py::arg("V"), py::arg("F"), py::arg("FN"), py::arg("VN"), py::arg("EN"), py::arg("EMAP"), py::arg("q"), py::arg("s"), py::arg("sqrd"), py::arg("i"), py::arg("c"), py::arg("n"));

//m.def("signed_distance_pseudonormal", []
//(
//  const AABB<Eigen::MatrixXd, 2> & tree,
//  const Eigen::MatrixXd& V,
//  const Eigen::MatrixXi& F,
//  const Eigen::MatrixXd& FN,
//  const Eigen::MatrixXd& VN,
//  const Eigen::MatrixXd & q,
//  double & s,
//  double & sqrd,
//  int & i,
//  Eigen::MatrixXd & c,
//  Eigen::MatrixXd & n
//)
//{
//  assert_is_VectorX("q",q);
//  return igl::signed_distance_pseudonormal(tree, V, F, FN, VN, q, s, sqrd, i, c, n);
//}, __doc_igl_signed_distance_pseudonormal,
//py::arg("tree"), py::arg("V"), py::arg("F"), py::arg("FN"), py::arg("VN"), py::arg("q"), py::arg("s"), py::arg("sqrd"), py::arg("i"), py::arg("c"), py::arg("n"));

//m.def("signed_distance_winding_number", []
//(
//  AABB<Eigen::MatrixXd, 3> & tree,
//  const Eigen::MatrixXd& V,
//  const Eigen::MatrixXi& F,
//  igl::WindingNumberAABB<Eigen::Vector3d> & hier,
//  Eigen::RowVector3d & q
//)
//{
//  return igl::signed_distance_winding_number(tree, V, F, hier, q);
//}, __doc_igl_signed_distance_winding_number,
//py::arg("tree"), py::arg("V"), py::arg("F"), py::arg("hier"), py::arg("q"));

