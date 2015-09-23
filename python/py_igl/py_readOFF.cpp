m.def("readOFF", []
(
  const std::string str,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F
)
{
  return igl::readOFF(str,V,F);
}, __doc_igl_readOFF,
py::arg("str"), py::arg("V"), py::arg("F"));

m.def("readOFF", []
(
  const std::string str,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F,
  Eigen::MatrixXd& N
)
{
  return igl::readOFF(str,V,F,N);
}, __doc_igl_readOFF,
py::arg("str"), py::arg("V"), py::arg("F"), py::arg("N"));
