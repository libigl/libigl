// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("writeMESH", []
(
  const std::string mesh_file_name,
  const std::vector<std::vector<double> > & V,
  const std::vector<std::vector<int> > & T,
  const std::vector<std::vector<int> > & F
)
{
  return igl::writeMESH(mesh_file_name, V, T, F);
}, __doc_igl_writeMESH,
py::arg("mesh_file_name"), py::arg("V"), py::arg("T"), py::arg("F"));

m.def("writeMESH", []
(
  const std::string str,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXd& T,
  const Eigen::MatrixXi& F
)
{
  return igl::writeMESH(str, V, T, F);
}, __doc_igl_writeMESH,
py::arg("str"), py::arg("V"), py::arg("T"), py::arg("F"));

