m.def("lscm", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXi& b,
  const Eigen::MatrixXd& bc,
  Eigen::MatrixXd& V_uv
)
{
  assert_is_VectorX("b",b);
  return igl::lscm(V,F,b,bc,V_uv);
}, __doc_igl_lscm,
py::arg("V"), py::arg("F"), py::arg("b"), py::arg("bc"), py::arg("V_uv"));
