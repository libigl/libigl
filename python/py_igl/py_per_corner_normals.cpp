// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("per_corner_normals", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const double corner_threshold,
  Eigen::MatrixXd& CN
)
{
  return igl::per_corner_normals(V,F,corner_threshold,CN);
}, __doc_igl_per_corner_normals,
py::arg("V"), py::arg("F"), py::arg("corner_threshold"), py::arg("CN"));

m.def("per_corner_normals", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& FN,
  const double corner_threshold,
  Eigen::MatrixXd& CN
)
{
  return igl::per_corner_normals(V,F,FN,corner_threshold,CN);
}, __doc_igl_per_corner_normals,
py::arg("V"), py::arg("F"), py::arg("FN"), py::arg("corner_threshold"), py::arg("CN"));

m.def("per_corner_normals", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& FN,
  const double corner_threshold,
  const std::vector<std::vector<int> >& VF,
  Eigen::MatrixXd& CN
)
{
  return igl::per_corner_normals(V,F,FN,VF,corner_threshold,CN);
}, __doc_igl_per_corner_normals,
py::arg("V"), py::arg("F"), py::arg("FN"), py::arg("corner_threshold"), py::arg("VF"), py::arg("CN"));
