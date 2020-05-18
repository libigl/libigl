// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2019 Thomas Davies thomasryantrd@gmail.com
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

py::class_<igl::FastWindingNumberBVH> FastWindingNumberBVH(m,"FastWindingNumberBVH");
FastWindingNumberBVH.def(py::init<>());

//we only expose precomputation of tree 
// for fast_winding_number signing in signed distance
m.def("fast_winding_number", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const int order,
  igl::FastWindingNumberBVH & fwn_bvh
)
{
  igl::fast_winding_number(V.cast<float>(), F, order, fwn_bvh);
}, __doc_igl_signed_distance,
py::arg("V"), py::arg("F"), py::arg("order"), py::arg("fwn_bvh"));