m.def("writePLY", []
(
  const std::string str,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& N,
  const Eigen::MatrixXd& UV
)
{
  return igl::writePLY(str,V,F,N,UV);
}, __doc_igl_writePLY,
py::arg("str"), py::arg("V"), py::arg("F"), py::arg("N"), py::arg("UV"));

m.def("writePLY", []
(
  const std::string str,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F
)
{
  return igl::writePLY(str,V,F);
}, __doc_igl_writePLY,
py::arg("str"), py::arg("V"), py::arg("F"));
