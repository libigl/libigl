// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Amrollah Seifoddini <a.seif67@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

// Wrap the data class, no properties are exposed since it is not necessary
py::class_<igl::HeatGeodesicsData<double> > heat_geodesics_data(m, "heat_geodesics_data");

heat_geodesics_data.def(py::init<>());

m.def("heat_geodesics_precompute", []
(
    const Eigen::MatrixXd &V,
    const Eigen::MatrixXi &F,
    igl::HeatGeodesicsData<double> &data
)
{
  return igl::heat_geodesics_precompute(V, F, data);
}, __doc_igl_heat_geodesics_precompute,
py::arg("V"), py::arg("F"), py::arg("data"));


m.def("heat_geodesics_solve", []
(
    const igl::HeatGeodesicsData<double> &data,
    const Eigen::MatrixXi &gamma,
    Eigen::MatrixXd &D
)
{
assert_is_VectorX("D", D);
assert_is_VectorX("gamma", gamma);
Eigen::VectorXd vD;
igl::heat_geodesics_solve(data, Eigen::VectorXi(gamma), vD);
D.resize(vD.size(), 1);
for (int i = 0; i < vD.size(); i++)
{
    D(i, 0) = vD[i];
}
return true;
}, __doc_igl_heat_geodesics_solve,
py::arg("data"), py::arg("gamma"), py::arg("D"));

