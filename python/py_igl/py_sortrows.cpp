// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("sortrows", []
(
  const Eigen::MatrixXd& X,
  const bool ascending,
  Eigen::MatrixXd& Y,
  Eigen::MatrixXi& I
)
{
  return igl::sortrows(X,ascending,Y,I);
}, __doc_igl_sortrows,
py::arg("X"), py::arg("ascending"), py::arg("Y"), py::arg("I"));

m.def("sortrows", []
(
  const Eigen::MatrixXi& X,
  const bool ascending,
  Eigen::MatrixXi& Y,
  Eigen::MatrixXi& I
)
{
  return igl::sortrows(X,ascending,Y,I);
}, __doc_igl_sortrows,
py::arg("X"), py::arg("ascending"), py::arg("Y"), py::arg("I"));
