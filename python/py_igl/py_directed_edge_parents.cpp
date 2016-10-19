

m.def("directed_edge_parents", []
(
  const Eigen::MatrixXi& E,
  Eigen::MatrixXi& P
)
{
  Eigen::VectorXi Pv;
  igl::directed_edge_parents(E, Pv);
  P = Pv;
}, __doc_igl_directed_edge_parents,
py::arg("E"), py::arg("P"));

