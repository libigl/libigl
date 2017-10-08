// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("read_triangle_mesh", []
(
  const std::string str,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F
)
{
  return igl::read_triangle_mesh(str,V,F);
}, __doc_igl_read_triangle_mesh,
py::arg("str"), py::arg("V"), py::arg("F"));

m.def("read_triangle_mesh", []
(
  const std::string str,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F,
  std::string & dir,
  std::string & base,
  std::string & ext,
  std::string & name
)
{
  return igl::read_triangle_mesh(str,V,F,dir,base,ext,name);
}, __doc_igl_read_triangle_mesh,
py::arg("str"), py::arg("V"), py::arg("F"), py::arg("dir"), py::arg("base"), py::arg("ext"), py::arg("name"));

m.def("read_triangle_mesh", []
(
  const std::string str,
  std::vector<std::vector<double> >& V,
  std::vector<std::vector<int> >& F
)
{
  return igl::read_triangle_mesh(str,V,F);
}, __doc_igl_read_triangle_mesh,
py::arg("str"), py::arg("V"), py::arg("F"));
