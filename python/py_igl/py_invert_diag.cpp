m.def("invert_diag", []
(
  const Eigen::SparseMatrix<double>& X,
  Eigen::SparseMatrix<double>& Y
)
{
  return igl::invert_diag(X,Y);
}, __doc_igl_invert_diag,
py::arg("X"), py::arg("Y"));
