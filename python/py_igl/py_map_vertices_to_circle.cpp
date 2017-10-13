// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("map_vertices_to_circle", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& bnd,
  Eigen::MatrixXd& UV
)
{
  assert_is_VectorX("bnd",bnd);
  return igl::map_vertices_to_circle(V,bnd,UV);
}, __doc_igl_map_vertices_to_circle,
py::arg("V"), py::arg("bnd"), py::arg("UV"));
