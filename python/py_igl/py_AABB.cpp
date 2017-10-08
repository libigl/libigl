// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
py::class_<igl::AABB<Eigen::MatrixXd,3> > AABB(m, "AABB");

AABB
.def(py::init<>())
.def(py::init<const igl::AABB<Eigen::MatrixXd,3>& >())
.def("init",[](igl::AABB<Eigen::MatrixXd,3>& tree, const Eigen::MatrixXd& V, const Eigen::MatrixXi& Ele)
{
    return tree.init(V, Ele, Eigen::Matrix<double, Eigen::Dynamic, 3>(), Eigen::Matrix<double, Eigen::Dynamic, 3>(), Eigen::VectorXi(), 0); 
})
.def("squared_distance", [](const igl::AABB<Eigen::MatrixXd,3>& tree, const Eigen::MatrixXd& V, const Eigen::MatrixXi& Ele, const Eigen::MatrixXd& P, Eigen::MatrixXd& sqrD, Eigen::MatrixXi& I, Eigen::MatrixXd& C)
{
    return tree.squared_distance(V, Ele, P, sqrD, I, C);
})
;
