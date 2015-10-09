m.def("local_basis", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXd& B1,
  Eigen::MatrixXd& B2,
  Eigen::MatrixXd& B3
)
{
  return igl::local_basis(V,F,B1,B2,B3);
}, __doc_igl_local_basis,
py::arg("V"), py::arg("F"), py::arg("B1"), py::arg("B2"), py::arg("B3"));
