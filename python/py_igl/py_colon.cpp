// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("colon", []
(
  const double low,
  const double step,
  const double high,
  Eigen::MatrixXd& I
)
{
  Eigen::Matrix<double,Eigen::Dynamic,1> temp;
  igl::colon<double>(low,step,high,temp);
  I = temp;
}, __doc_igl_colon,
py::arg("low"), py::arg("step"), py::arg("high"), py::arg("I"));

m.def("colon", []
(
  const double low,
  const double high,
  Eigen::MatrixXd& I
)
{
  Eigen::Matrix<double,Eigen::Dynamic,1> temp;
  igl::colon<double>(low,high,temp);
  I = temp;
}, __doc_igl_colon,
py::arg("low"), py::arg("high"), py::arg("I"));

m.def("colon", []
(
  const double& low,
  const double& high
)
{
  return Eigen::MatrixXd(igl::colon<double>(low,high));
}, __doc_igl_colon,
py::arg("low"), py::arg("high"));


m.def("coloni", []
(
  const int low,
  const int step,
  const int high,
  Eigen::MatrixXi& I
)
{
  Eigen::Matrix<int,Eigen::Dynamic,1> temp;
  igl::colon<int>(low,step,high,temp);
  I = temp;
}, __doc_igl_colon,
py::arg("low"), py::arg("step"), py::arg("high"), py::arg("I"));

m.def("coloni", []
(
  const int low,
  const int high,
  Eigen::MatrixXi& I
)
{
  Eigen::Matrix<int,Eigen::Dynamic,1> temp;
  igl::colon<int>(low,high,temp);
  I = temp;
}, __doc_igl_colon,
py::arg("low"), py::arg("high"), py::arg("I"));

m.def("coloni", []
(
  const int& low,
  const int& high
)
{
  return Eigen::MatrixXi(igl::colon<int>(low,high));
}, __doc_igl_colon,
py::arg("low"), py::arg("high"));
