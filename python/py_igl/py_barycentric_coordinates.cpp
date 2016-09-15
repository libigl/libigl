

m.def("barycentric_coordinates", []
(
  const Eigen::MatrixXd& P,
  const Eigen::MatrixXd& A,
  const Eigen::MatrixXd& B,
  const Eigen::MatrixXd& C,
  const Eigen::MatrixXd& D,
  Eigen::MatrixXd& L
)
{
  return igl::barycentric_coordinates(P, A, B, C, D, L);
}, __doc_igl_barycentric_coordinates,
py::arg("P"), py::arg("A"), py::arg("B"), py::arg("C"), py::arg("D"), py::arg("L"));

m.def("barycentric_coordinates", []
(
  const Eigen::MatrixXd& P,
  const Eigen::MatrixXd& A,
  const Eigen::MatrixXd& B,
  const Eigen::MatrixXd& C,
  Eigen::MatrixXd& L
)
{
  return igl::barycentric_coordinates(P, A, B, C, L);
}, __doc_igl_barycentric_coordinates,
py::arg("P"), py::arg("A"), py::arg("B"), py::arg("C"), py::arg("L"));

