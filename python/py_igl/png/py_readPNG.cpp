// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

m.def("readPNG", []
(
  const std::string png_file,
  Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> & R,
  Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> & G,
  Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> & B,
  Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> & A
)
{
  return igl::png::readPNG(png_file, R, G, B, A);
}, __doc_igl_png_readPNG,
py::arg("png_file"), py::arg("R"), py::arg("G"), py::arg("B"), py::arg("A"));

