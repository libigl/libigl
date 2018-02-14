// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("nrosy", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXi& b,
  const Eigen::MatrixXd& bc,
  const Eigen::MatrixXi& b_soft,
  const Eigen::MatrixXd& w_soft,
  const Eigen::MatrixXd& bc_soft,
  const int N,
  const double soft,
  Eigen::MatrixXd& R,
  Eigen::MatrixXd& S
)
{
  assert_is_VectorX("b",b);
  assert_is_VectorX("b_soft",b_soft);
  assert_is_VectorX("w_soft",w_soft);

  Eigen::VectorXi bt;
  if (b.size() != 0)
    bt = b;

  Eigen::VectorXi b_softt;
  if (b_soft.size() != 0)
    b_softt = b_soft;

  Eigen::VectorXd w_softt;
  if (w_soft.size() != 0)
    w_softt = w_soft;

  Eigen::VectorXd St;
  igl::copyleft::comiso::nrosy(V,F,bt,bc,b_softt,w_softt,bc_soft,N,soft,R,St);
  S = St;

}, __doc_igl_copyleft_comiso_nrosy,
py::arg("V"), py::arg("F"), py::arg("b"), py::arg("bc"), py::arg("b_soft"), py::arg("w_soft"), py::arg("bc_soft"), py::arg("N"), py::arg("soft"), py::arg("R"), py::arg("S"));

m.def("nrosy", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXi& b,
  const Eigen::MatrixXd& bc,
  const int N,
  Eigen::MatrixXd& R,
  Eigen::MatrixXd& S
)
{
  assert_is_VectorX("b",b);

  Eigen::VectorXi bt;
  if (b.size() != 0)
    bt = b;

  Eigen::VectorXd St;
  igl::copyleft::comiso::nrosy(V,F,bt,bc,N,R,St);
  S = St;
}, __doc_igl_copyleft_comiso_nrosy,
py::arg("V"), py::arg("F"), py::arg("b"), py::arg("bc"), py::arg("N"), py::arg("R"), py::arg("S"));
