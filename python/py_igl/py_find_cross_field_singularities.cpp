// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("find_cross_field_singularities", []
(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  const Eigen::MatrixXi &Handle_MMatch,
  Eigen::MatrixXi &isSingularity,
  Eigen::MatrixXi &singularityIndex
)
{
  return igl::find_cross_field_singularities(V,F,Handle_MMatch,isSingularity,singularityIndex);
}, __doc_igl_find_cross_field_singularities,
py::arg("V"), py::arg("F"), py::arg("Handle_MMatch"), py::arg("isSingularity"), py::arg("singularityIndex"));

m.def("find_cross_field_singularities", []
(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  const Eigen::MatrixXd &PD1,
  const Eigen::MatrixXd &PD2,
  Eigen::MatrixXi &isSingularity,
  Eigen::MatrixXi &singularityIndex,
  bool isCombed
)
{
  return igl::find_cross_field_singularities(V,F,PD1,PD2,isSingularity,singularityIndex,isCombed);
}, __doc_igl_find_cross_field_singularities,
py::arg("V"), py::arg("F"), py::arg("PD1"), py::arg("PD2"), py::arg("isSingularity"), py::arg("singularityIndex"),  py::arg("isCombed") = false);
