// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
// Wrap the params struct
py::class_<igl::active_set_params > active_set_params(m, "active_set_params");

active_set_params
.def("__init__", [](igl::active_set_params &m)
{
    new (&m) igl::active_set_params();
    m.Auu_pd = false;
    m.max_iter = 100;
    m.inactive_threshold = igl::DOUBLE_EPS;
    m.constraint_threshold = igl::DOUBLE_EPS;
    m.solution_diff_threshold = igl::DOUBLE_EPS;
})
.def_readwrite("Auu_pd", &igl::active_set_params::Auu_pd)
.def_readwrite("max_iter", &igl::active_set_params::max_iter)
.def_readwrite("inactive_threshold", &igl::active_set_params::inactive_threshold)
.def_readwrite("constraint_threshold", &igl::active_set_params::constraint_threshold)
.def_readwrite("solution_diff_threshold", &igl::active_set_params::solution_diff_threshold)
.def_readwrite("Auu_pd", &igl::active_set_params::Auu_pd)
;

m.def("active_set", []
(
  const Eigen::SparseMatrix<double>& A,
  const Eigen::MatrixXd& B,
  const Eigen::MatrixXi& known,
  const Eigen::MatrixXd& Y,
  const Eigen::SparseMatrix<double>& Aeq,
  const Eigen::MatrixXd& Beq,
  const Eigen::SparseMatrix<double>& Aieq,
  const Eigen::MatrixXd& Bieq,
  const Eigen::MatrixXd& lx,
  const Eigen::MatrixXd& ux,
  const igl::active_set_params& params,
  Eigen::MatrixXd& Z
)
{
  assert_is_VectorX("B",B);
  assert_is_VectorX("known",known);
  assert_is_VectorX("Y",Y);
  assert_is_VectorX("Beq",Beq);
  assert_is_VectorX("Bieq",Bieq);
  assert_is_VectorX("Z",Z);

  return igl::active_set(A,B,known,Y,Aeq,Eigen::VectorXd(Beq),Aieq,Bieq,lx,ux,params,Z);
}, __doc_igl_active_set,
py::arg("A"), py::arg("B"), py::arg("known"), py::arg("Y"), py::arg("Aeq"), py::arg("Beq")
, py::arg("Aieq"), py::arg("Bieq"), py::arg("lx"), py::arg("ux"), py::arg("params"), py::arg("Z"));
