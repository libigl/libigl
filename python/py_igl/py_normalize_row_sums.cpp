

m.def("normalize_row_sums", []
(
  const Eigen::MatrixXd& A,
  Eigen::MatrixXd& B
)
{
  return igl::normalize_row_sums(A, B);
}, __doc_igl_normalize_row_sums,
py::arg("A"), py::arg("B"));

