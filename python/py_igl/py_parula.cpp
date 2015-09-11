m.def("parula", []
(
  const Eigen::MatrixXd& Z,
  const bool normalize,
  Eigen::MatrixXd& C
)
{
  assert_is_VectorX("Z",Z);
  return igl::parula(Z,normalize,C);
}, __doc_igl_parula,
py::arg("Z"), py::arg("normalize"), py::arg("C"));

m.def("parula", []
(
  const Eigen::MatrixXd& Z,
  const double min_Z,
  const double max_Z,
  Eigen::MatrixXd& C
)
{
  assert_is_VectorX("Z",Z);
  return igl::parula(Z,min_Z,max_Z,C);
}, __doc_igl_parula,
py::arg("Z"), py::arg("min_Z"), py::arg("max_Z"), py::arg("C"));
