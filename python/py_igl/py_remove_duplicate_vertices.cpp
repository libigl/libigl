m.def("remove_duplicate_vertices", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const double epsilon,
  Eigen::MatrixXd& SV,
  Eigen::MatrixXi& SVI,
  Eigen::MatrixXi& SVJ,
  Eigen::MatrixXi& SF
)
{
  return igl::remove_duplicate_vertices(V, F, epsilon, SV, SVI, SVJ, SF);
}, __doc_igl_remove_duplicate_vertices,
py::arg("V"), py::arg("F"), py::arg("epsilon"), py::arg("SV"), py::arg("SVI"), py::arg("SVJ"), py::arg("SF"));
