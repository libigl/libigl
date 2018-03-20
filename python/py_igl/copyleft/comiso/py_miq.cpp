// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
m.def("miq", []
(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  const Eigen::MatrixXd &PD1,
  const Eigen::MatrixXd &PD2,
  Eigen::MatrixXd &UV,
  Eigen::MatrixXi &FUV,
  double scale,
  double stiffness,
  bool direct_round,
  int iter,
  int local_iter,
  bool DoRound,bool SingularityRound
  //  std::vector<int> round_vertices,
  //  std::vector<std::vector<int> > hard_features
)
{
  std::vector<int> round_vertices;
  std::vector<std::vector<int> > hard_features;

  igl::copyleft::comiso::miq(V,F,PD1,PD2,UV,FUV,scale,stiffness,direct_round,iter,local_iter,DoRound, SingularityRound, round_vertices, hard_features);
}, __doc_igl_copyleft_comiso_miq,
py::arg("V"), py::arg("F"), py::arg("PD1"), py::arg("PD2"), py::arg("UV"), py::arg("FUV"), py::arg("scale") = 30.0, py::arg("stiffness") = 5.0, py::arg("direct_round") = false, py::arg("iter") = 5, py::arg("local_iter") = 5, py::arg("DoRound") = true, py::arg("SingularityRound") = true
// , py::arg("round_vertices"), py::arg("hard_features")
);

m.def("miq", []
(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  const Eigen::MatrixXd &PD1_combed,
  const Eigen::MatrixXd &PD2_combed,
  const Eigen::MatrixXi &MMatch,
  const Eigen::MatrixXi &Singular,
  const Eigen::MatrixXi &Seams,
  Eigen::MatrixXd &UV,
  Eigen::MatrixXi &FUV,
  double GradientSize,
  double Stiffness,
  bool DirectRound,
  int iter,
  int localIter, bool DoRound, bool SingularityRound
  // std::vector<int> roundVertices,
  // std::vector<std::vector<int> > hardFeatures
)
{
  assert_is_VectorX("Singular",Singular);

  std::vector<int> roundVertices;
  std::vector<std::vector<int> > hardFeatures;

  igl::copyleft::comiso::miq(V,F,PD1_combed,PD2_combed,MMatch,Singular,Seams,UV,FUV,GradientSize,Stiffness,DirectRound,iter,localIter,DoRound, SingularityRound, roundVertices, hardFeatures);
}, __doc_igl_copyleft_comiso_miq,
py::arg("V"), py::arg("F"), py::arg("PD1_combed"), py::arg("PD2_combed"),
py::arg("MMatch"), py::arg("Singular"), py::arg("Seams"),
py::arg("UV"), py::arg("FUV"), py::arg("GradientSize") = 30.0, py::arg("Stiffness") = 5.0, py::arg("DirectRound") = false, py::arg("iter") = 5, py::arg("localIter") = 5, py::arg("DoRound") = true, py::arg("SingularityRound") = true
// , py::arg("roundVertices"), py::arg("hardFeatures")
);
