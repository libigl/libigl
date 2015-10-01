m.def("nrosy", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXi& b,
  const Eigen::MatrixXd& bc,
  const Eigen::MatrixXi& b_soft,
  const Eigen::MatrixXd& w_soft,
  const Eigen::MatrixXd& bc_soft,
  const int N,
  const double soft,
  Eigen::MatrixXd& R,
  Eigen::MatrixXd& S
)
{
  assert_is_VectorX("b",b);
  assert_is_VectorX("b_soft",b_soft);
  assert_is_VectorX("w_soft",w_soft);
  Eigen::VectorXd St;
  igl::comiso::nrosy(V,F,b,bc,b_soft,w_soft,bc_soft,N,soft,R,St);
  S = St;
}, __doc_igl_comiso_nrosy,
py::arg("V"), py::arg("F"), py::arg("b"), py::arg("bc"), py::arg("b_soft"), py::arg("w_soft"), py::arg("bc_soft"), py::arg("N"), py::arg("soft"), py::arg("R"), py::arg("S"));

m.def("nrosy", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXi& b,
  const Eigen::MatrixXd& bc,
  const int N,
  Eigen::MatrixXd& R,
  Eigen::MatrixXd& S
)
{
  assert_is_VectorX("b",b);
  Eigen::VectorXd St;
  igl::comiso::nrosy(V,F,b,bc,N,R,St);
  S = St;
}, __doc_igl_comiso_nrosy,
py::arg("V"), py::arg("F"), py::arg("b"), py::arg("bc"), py::arg("N"), py::arg("R"), py::arg("S"));
