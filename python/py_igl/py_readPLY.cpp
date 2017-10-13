// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("readPLY", []
(
  const std::string str,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F,
  Eigen::MatrixXd& N,
  Eigen::MatrixXd& UV
)
{
  return igl::readPLY(str,V,F,N,UV);
}, __doc_igl_readPLY,
py::arg("str"), py::arg("V"), py::arg("F"), py::arg("N"), py::arg("UV"));
