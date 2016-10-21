// COMPLETE BINDINGS ========================


m.def("lbs_matrix", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXd& W,
  Eigen::MatrixXd& M
)
{
  return igl::lbs_matrix(V, W, M);
}, __doc_igl_lbs_matrix,
py::arg("V"), py::arg("W"), py::arg("M"));

m.def("lbs_matrix_column", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXd& W,
  Eigen::MatrixXd& M
)
{
  return igl::lbs_matrix_column(V, W, M);
}, __doc_igl_lbs_matrix_column,
py::arg("V"), py::arg("W"), py::arg("M"));

m.def("lbs_matrix_column", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXd& W,
  const Eigen::MatrixXi& WI,
  Eigen::MatrixXd& M
)
{
  return igl::lbs_matrix_column(V, W, WI, M);
}, __doc_igl_lbs_matrix_column,
py::arg("V"), py::arg("W"), py::arg("WI"), py::arg("M"));





// INCOMPLETE BINDINGS ========================


m.def("lbs_matrix_column", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXd& W,
  Eigen::SparseMatrix<double>& M
)
{
  return igl::lbs_matrix_column(V, W, M);
}, __doc_igl_lbs_matrix_column,
py::arg("V"), py::arg("W"), py::arg("M"));

m.def("lbs_matrix_column", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXd& W,
  const Eigen::MatrixXi& WI,
  Eigen::SparseMatrix<double>& M
)
{
  return igl::lbs_matrix_column(V, W, WI, M);
}, __doc_igl_lbs_matrix_column,
py::arg("V"), py::arg("W"), py::arg("WI"), py::arg("M"));

