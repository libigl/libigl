// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
// COMPLETE BINDINGS ========================






// INCOMPLETE BINDINGS ========================


//m.def("collapse_edge", []
//(
//  int e,
//  Eigen::RowVectorXd & p,
//  Eigen::MatrixXd& V,
//  Eigen::MatrixXi& F,
//  Eigen::MatrixXi& E,
//  Eigen::MatrixXi& EMAP,
//  Eigen::MatrixXi& EF,
//  Eigen::MatrixXi& EI,
//  int & e1,
//  int & e2,
//  int & f1,
//  int & f2
//)
//{
//  return igl::collapse_edge(e, p, V, F, E, EMAP, EF, EI, e1, e2, f1, f2);
//}, __doc_igl_collapse_edge,
//py::arg("e"), py::arg("p"), py::arg("V"), py::arg("F"), py::arg("E"), py::arg("EMAP"), py::arg("EF"), py::arg("EI"), py::arg("e1"), py::arg("e2"), py::arg("f1"), py::arg("f2"));

//m.def("collapse_edge", []
//(
//  int e,
//  Eigen::RowVectorXd & p,
//  Eigen::MatrixXd& V,
//  Eigen::MatrixXi& F,
//  Eigen::MatrixXi& E,
//  Eigen::MatrixXi& EMAP,
//  Eigen::MatrixXi& EF,
//  Eigen::MatrixXi& EI
//)
//{
//  return igl::collapse_edge(e, p, V, F, E, EMAP, EF, EI);
//}, __doc_igl_collapse_edge,
//py::arg("e"), py::arg("p"), py::arg("V"), py::arg("F"), py::arg("E"), py::arg("EMAP"), py::arg("EF"), py::arg("EI"));

//m.def("collapse_edge", []
//(
//  std::function<void (const int, const Eigen::MatrixXd &, const Eigen::MatrixXi &, const Eigen::MatrixXi &, const Eigen::VectorXi &, const Eigen::MatrixXi &, const Eigen::MatrixXi &, double &, Eigen::RowVectorXd &)> & cost_and_placement,
//  Eigen::MatrixXd& V,
//  Eigen::MatrixXi& F,
//  Eigen::MatrixXi& E,
//  Eigen::MatrixXi& EMAP,
//  Eigen::MatrixXi& EF,
//  Eigen::MatrixXi& EI,
//  std::set<std::pair<double, int> > & Q,
//  std::vector<std::set<std::pair<double, int> >::iterator> & Qit,
//  Eigen::MatrixXd& C
//)
//{
//  return igl::collapse_edge(cost_and_placement, V, F, E, EMAP, EF, EI, Q, Qit, C);
//}, __doc_igl_collapse_edge,
//py::arg("cost_and_placement"), py::arg("V"), py::arg("F"), py::arg("E"), py::arg("EMAP"), py::arg("EF"), py::arg("EI"), py::arg("Q"), py::arg("Qit"), py::arg("C"));

//m.def("collapse_edge", []
//(
//  std::function<void (const int, const Eigen::MatrixXd &, const Eigen::MatrixXi &, const Eigen::MatrixXi &, const Eigen::VectorXi &, const Eigen::MatrixXi &, const Eigen::MatrixXi &, double &, Eigen::RowVectorXd &)> & cost_and_placement,
//  Eigen::MatrixXd& V,
//  Eigen::MatrixXi& F,
//  Eigen::MatrixXi& E,
//  Eigen::MatrixXi& EMAP,
//  Eigen::MatrixXi& EF,
//  Eigen::MatrixXi& EI,
//  std::set<std::pair<double, int> > & Q,
//  std::vector<std::set<std::pair<double, int> >::iterator> & Qit,
//  Eigen::MatrixXd& C,
//  int & e,
//  int & e1,
//  int & e2,
//  int & f1,
//  int & f2
//)
//{
//  return igl::collapse_edge(cost_and_placement, V, F, E, EMAP, EF, EI, Q, Qit, C, e, e1, e2, f1, f2);
//}, __doc_igl_collapse_edge,
//py::arg("cost_and_placement"), py::arg("V"), py::arg("F"), py::arg("E"), py::arg("EMAP"), py::arg("EF"), py::arg("EI"), py::arg("Q"), py::arg("Qit"), py::arg("C"), py::arg("e"), py::arg("e1"), py::arg("e2"), py::arg("f1"), py::arg("f2"));

