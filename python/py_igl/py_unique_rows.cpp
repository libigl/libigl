m.def("unique_rows", []
(
  const Eigen::MatrixXd& A,
  Eigen::MatrixXd& C,
  Eigen::MatrixXi& IA,
  Eigen::MatrixXi& IC
)
{
  return igl::unique_rows(A,C,IA,IC);
}, __doc_igl_unique,
py::arg("A"), py::arg("C"), py::arg("IA"), py::arg("IC"));

m.def("unique_rows", []
(
  const Eigen::MatrixXi& A,
  Eigen::MatrixXi& C,
  Eigen::MatrixXi& IA,
  Eigen::MatrixXi& IC
)
{
  return igl::unique_rows(A,C,IA,IC);
}, __doc_igl_unique,
py::arg("A"), py::arg("C"), py::arg("IA"), py::arg("IC"));

