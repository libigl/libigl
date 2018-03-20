// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.


m.def("reorient_facets_raycast", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  int rays_total,
  int rays_minimum,
  bool facet_wise,
  bool use_parity,
  bool is_verbose,
  Eigen::MatrixXi& I,
  Eigen::MatrixXi& C
)
{
  Eigen::VectorXi Iv;
  Eigen::VectorXi Cv;
  igl::embree::reorient_facets_raycast(V, F, rays_total, rays_minimum, facet_wise, use_parity, is_verbose, Iv, Cv);
  I = Iv;
  C = Cv;
}, __doc_igl_embree_reorient_facets_raycast,
py::arg("V"), py::arg("F"), py::arg("rays_total"), py::arg("rays_minimum"), py::arg("facet_wise"), py::arg("use_parity"), py::arg("is_verbose"), py::arg("I"), py::arg("C"));

m.def("reorient_facets_raycast", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXi& FF,
  Eigen::MatrixXi& I
)
{
  Eigen::VectorXi Iv;
  igl::embree::reorient_facets_raycast(V, F, FF, Iv);
  I = Iv;
}, __doc_igl_embree_reorient_facets_raycast,
py::arg("V"), py::arg("F"), py::arg("FF"), py::arg("I"));

