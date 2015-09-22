m.def("cotmatrix", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::SparseMatrix<double>& L
)
{
  return igl::cotmatrix(V,F,L);
}, __doc_igl_cotmatrix,
py::arg("V"), py::arg("F"), py::arg("L"));
