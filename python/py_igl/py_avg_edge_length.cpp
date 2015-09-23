m.def("avg_edge_length", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F
)
{
  return igl::avg_edge_length(V,F);
}, __doc_igl_avg_edge_length,
py::arg("V"), py::arg("F"));
