// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("rotate_vectors", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXd& A,
  const Eigen::MatrixXd& B1,
  const Eigen::MatrixXd& B2
)
{
  assert_is_VectorX("A",A);
  return igl::rotate_vectors(V,A,B1,B2);
}, __doc_igl_rotate_vectors,
py::arg("V"), py::arg("A"), py::arg("B1"), py::arg("B2"));
