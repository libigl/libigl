

m.def("barycentric_to_global", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& bc
)
{
  return igl::barycentric_to_global(V, F, bc);
}, __doc_igl_barycentric_to_global,
py::arg("V"), py::arg("F"), py::arg("bc"));
