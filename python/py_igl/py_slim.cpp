// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Zhongshi Jiang <jiangzs@nyu.edu>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

m.def("slim", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& V_init,
  igl::SLIMData::SLIM_ENERGY slim_energy,
  Eigen::VectorXi& b,
  Eigen::MatrixXd& bc,
  double soft_p,
  int iter_num,
  Eigen::MatrixXd& V_out
)
{
  return igl::slim_solve(V,F,V_init, slim_energy, b, bc, soft_p, iter_num, V_out);
}, __doc_igl_slim, py::arg("V"), py::arg("F"), py::arg("V_init"), 
py::arg("slim_energy"), py::arg("b"), py::arg("bc"), 
py::arg("soft_p"), py::arg("iter_num"), py::arg("V_out")
);