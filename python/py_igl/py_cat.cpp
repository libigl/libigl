m.def("cat", []
(
  const int dim,
  const Eigen::MatrixXd& A,
  const Eigen::MatrixXd& B,
  Eigen::MatrixXd& C
)
{
  return igl::cat(dim, A, B, C);
}, __doc_igl_cat,
py::arg("dim"), py::arg("A"), py::arg("B"), py::arg("C"));

m.def("cat", []
(
  const int dim,
  Eigen::MatrixXd& A,
  Eigen::MatrixXd& B
)
{
  return igl::cat(dim, A, B);
}, __doc_igl_cat,
py::arg("dim"), py::arg("A"), py::arg("B"));

m.def("cat", []
(
  const int dim,
  Eigen::MatrixXi& A,
  Eigen::MatrixXi& B
)
{
  return igl::cat(dim, A, B);
}, __doc_igl_cat,
py::arg("dim"), py::arg("A"), py::arg("B"));

//m.def("cat", []
//(
//  const std::vector<std::vector<Eigen::MatrixXd > > & A, 
//  Eigen::MatrixXd & C
//)
//{
//  return igl::cat(A, C);
//}, __doc_igl_cat,
//py::arg("A"), py::arg("C"));

m.def("cat", []
(
  const int dim,
  const Eigen::SparseMatrix<double>& A,
  const Eigen::SparseMatrix<double>& B,
  Eigen::SparseMatrix<double>& C
)
{
  return igl::cat(dim, A, B, C);
}, __doc_igl_cat,
py::arg("dim"), py::arg("A"), py::arg("B"), py::arg("C"));

