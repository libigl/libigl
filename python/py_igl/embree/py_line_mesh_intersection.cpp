// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


m.def("line_mesh_intersection", []
(
  const Eigen::MatrixXd& V_source,
  const Eigen::MatrixXd& N_source,
  const Eigen::MatrixXd& V_target,
  const Eigen::MatrixXi& F_target
)
{
  return igl::embree::line_mesh_intersection(V_source, N_source, V_target, F_target);
}, __doc_igl_embree_line_mesh_intersection,
py::arg("V_source"), py::arg("N_source"), py::arg("V_target"), py::arg("F_target"));
