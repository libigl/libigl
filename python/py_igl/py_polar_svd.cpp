m.def("polar_svd", []
(
  const Eigen::MatrixXd& A,
  Eigen::MatrixXd& R,
  Eigen::MatrixXd& T,
  Eigen::MatrixXd& U,
  Eigen::MatrixXd& S,
  Eigen::MatrixXd& V
)
{
  return igl::polar_svd(A, R, T, U, S, V);
}, __doc_igl_polar_svd,
py::arg("A"), py::arg("R"), py::arg("T"), py::arg("U"), py::arg("S"), py::arg("V"));


m.def("polar_svd", []
(
  const Eigen::MatrixXd& A,
  Eigen::MatrixXd& R,
  Eigen::MatrixXd& T
)
{
  return igl::polar_svd(A, R, T);
}, __doc_igl_polar_svd,
py::arg("A"), py::arg("R"), py::arg("T"));

