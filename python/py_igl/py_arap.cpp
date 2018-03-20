// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
py::class_<igl::ARAPData> ARAPData(m, "ARAPData");

ARAPData
.def(py::init<>())
.def_readwrite("n", &igl::ARAPData::n)
.def_readwrite("energy", &igl::ARAPData::energy)
.def_property("G",
[](const igl::ARAPData& data) {return Eigen::MatrixXi(data.G);},
[](igl::ARAPData& data, const Eigen::MatrixXi& G)
{
  assert_is_VectorX("G",G);
  data.G = Eigen::VectorXi(G.cast<int>());
})
.def_readwrite("with_dynamics", &igl::ARAPData::with_dynamics)
.def_readwrite("f_ext", &igl::ARAPData::f_ext)
.def_readwrite("h", &igl::ARAPData::h)
.def_readwrite("vel", &igl::ARAPData::vel)
.def_readwrite("ym", &igl::ARAPData::ym)
.def_readwrite("max_iter", &igl::ARAPData::max_iter)
.def_readwrite("K", &igl::ARAPData::K)
.def_readwrite("M", &igl::ARAPData::M)
.def_readwrite("CSM", &igl::ARAPData::CSM)
// .def_readwrite("solver_data", &igl::ARAPData::solver_data)
.def_readwrite("dim", &igl::ARAPData::dim)
.def_property("b",
[](const igl::ARAPData& data) {return Eigen::MatrixXi(data.b);},
[](igl::ARAPData& data, const Eigen::MatrixXi& b)
{
  assert_is_VectorX("b",b);
  data.b = Eigen::VectorXi(b.cast<int>());
})
;

m.def("arap_precomputation", []
(
  const Eigen::MatrixXd & V,
  const Eigen::MatrixXi & F,
  const int dim,
  const Eigen::MatrixXi& b,
  igl::ARAPData & data
)
{
  assert_is_VectorX("b",b);
  Eigen::VectorXi bt;
  if (b.size() != 0)
    bt = b;

  return igl::arap_precomputation(V,F,dim,bt,data);
}, __doc_igl_arap_precomputation,
py::arg("V"), py::arg("F"), py::arg("dim"), py::arg("b"), py::arg("data"));

m.def("arap_solve", []
(
  const Eigen::MatrixXd & bc,
  igl::ARAPData & data,
  Eigen::MatrixXd& U
)
{
  return igl::arap_solve(bc,data,U);
}, __doc_igl_arap_solve,
py::arg("bc"), py::arg("data"), py::arg("U"));
