// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("readOFF", []
(
  const std::string str,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F
)
{
  return igl::readOFF(str,V,F);
}, __doc_igl_readOFF,
py::arg("str"), py::arg("V"), py::arg("F"));

m.def("readOFF", []
(
  const std::string str,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F,
  Eigen::MatrixXd& N
)
{
  return igl::readOFF(str,V,F,N);
}, __doc_igl_readOFF,
py::arg("str"), py::arg("V"), py::arg("F"), py::arg("N"));
