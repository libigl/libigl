// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("readMESH", []
(
  const std::string mesh_file_name,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& T,
  Eigen::MatrixXi& F
)
{
  return igl::readMESH(mesh_file_name, V, T, F);
}, __doc_igl_readMESH,
py::arg("mesh_file_name"), py::arg("V"), py::arg("T"), py::arg("F"));


m.def("readMESH", []
(
  const std::string mesh_file_name,
  std::vector<std::vector<double> > & V,
  std::vector<std::vector<int> > & T,
  std::vector<std::vector<int> > & F
)
{
  return igl::readMESH(mesh_file_name, V, T, F);
}, __doc_igl_readMESH,
py::arg("mesh_file_name"), py::arg("V"), py::arg("T"), py::arg("F"));



