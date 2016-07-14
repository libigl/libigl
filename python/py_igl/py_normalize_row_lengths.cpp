

m.def("normalize_row_lengths", []
(
  const Eigen::MatrixXd& A,
  Eigen::MatrixXd& B
)
{
  return igl::normalize_row_lengths(A, B);
}, __doc_igl_normalize_row_lengths,
py::arg("A"), py::arg("B"));

