m.def("harmonic", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXi& b,
  const Eigen::MatrixXd& bc,
  const int k,
  Eigen::MatrixXd& W
)
{
  assert_is_VectorX("b",b);
  return igl::harmonic(V,F,b,bc,k,W);
}, __doc_igl_harmonic,
py::arg("V"), py::arg("F"), py::arg("b"), py::arg("bc"), py::arg("k"), py::arg("W"));
