// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("writeOBJ", []
(
  const std::string str,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& CN,
  const Eigen::MatrixXi& FN,
  const Eigen::MatrixXd& TC,
  const Eigen::MatrixXi& FTC
)
{
  return igl::writeOBJ(str,V,F,CN,FN,TC,FTC);
}, __doc_igl_writeOBJ,
py::arg("str"), py::arg("V"), py::arg("F"), py::arg("CN"), py::arg("FN"), py::arg("TC"), py::arg("FTC"));

m.def("writeOBJ", []
(
  const std::string str,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F
)
{
  return igl::writeOBJ(str,V,F);
}, __doc_igl_writeOBJ,
py::arg("str"), py::arg("V"), py::arg("F"));
