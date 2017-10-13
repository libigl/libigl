// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


m.def("ambient_occlusion", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& P,
  const Eigen::MatrixXd& N,
  const int num_samples,
  Eigen::MatrixXd& S
)
{
  return igl::embree::ambient_occlusion(V, F, P, N, num_samples, S);
}, __doc_igl_embree_ambient_occlusion,
py::arg("V"), py::arg("F"), py::arg("P"), py::arg("N"), py::arg("num_samples"), py::arg("S"));

