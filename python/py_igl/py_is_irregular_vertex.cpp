

m.def("is_irregular_vertex", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F
)
{
  return igl::is_irregular_vertex(V, F);
}, __doc_igl_is_irregular_vertex,
py::arg("V"), py::arg("F"));

