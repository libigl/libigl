m.def("edge_lengths", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXd& L
)
{
  return igl::edge_lengths(V, F, L);
}, __doc_igl_edge_lengths,
py::arg("V"), py::arg("F"), py::arg("L"));

