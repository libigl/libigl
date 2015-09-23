m.def("setdiff", []
(
  const Eigen::MatrixXi& A,
  const Eigen::MatrixXi& B,
  Eigen::MatrixXi& C,
  Eigen::MatrixXi& IA
)
{
  return igl::setdiff(A,B,C,IA);
}, __doc_igl_setdiff,
py::arg("A"), py::arg("B"), py::arg("C"), py::arg("IA"));
