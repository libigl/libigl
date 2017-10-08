// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


m.def("marching_cubes", []
(
  const Eigen::MatrixXd& values,
  const Eigen::MatrixXd& points,
  const unsigned int x_res,
  const unsigned int y_res,
  const unsigned int z_res,
  Eigen::MatrixXd& vertices,
  Eigen::MatrixXi& faces
)
{
  assert_is_VectorX("values", values);
  Eigen::VectorXd valuesv;
  if (values.size() != 0)
    valuesv = values;
  return igl::copyleft::marching_cubes(valuesv, points, x_res, y_res, z_res, vertices, faces);
}, __doc_igl_copyleft_marching_cubes,
py::arg("values"), py::arg("points"), py::arg("x_res"), py::arg("y_res"), py::arg("z_res"), py::arg("vertices"), py::arg("faces"));

