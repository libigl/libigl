// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("readOBJ", []
(
  const std::string str,
  Eigen::MatrixXd& V,
  Eigen::MatrixXd& TC,
  Eigen::MatrixXd& CN,
  Eigen::MatrixXi& F,
  Eigen::MatrixXi& FTC,
  Eigen::MatrixXi& FN
)
{
  return igl::readOBJ(str,V,TC,CN,F,FTC,FN);
}, __doc_igl_readOBJ,
py::arg("str"), py::arg("V"), py::arg("TC"), py::arg("CN"), py::arg("F"), py::arg("FTC"), py::arg("FN"));

m.def("readOBJ", []
(
  const std::string str,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F
)
{
  return igl::readOBJ(str,V,F);
}, __doc_igl_readOBJ,
py::arg("str"), py::arg("V"), py::arg("F"));
