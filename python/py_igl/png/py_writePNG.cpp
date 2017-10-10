// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


m.def("writePNG", []
(
  const Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> & R,
  const Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> & G,
  const Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> & B,
  const Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> & A,
  const std::string png_file
)
{
  return igl::png::writePNG(R, G, B, A, png_file);
}, __doc_igl_png_writePNG,
py::arg("R"), py::arg("G"), py::arg("B"), py::arg("A"), py::arg("png_file"));

