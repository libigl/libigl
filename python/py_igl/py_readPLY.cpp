m.def("readPLY", []
(
  const std::string str,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F,
  Eigen::MatrixXd& N, 
  Eigen::MatrixXd& UV
)
{
  return igl::readPLY(str,V,F,N,UV);
}, __doc_igl_readPLY, 
py::arg("str"), py::arg("V"), py::arg("F"), py::arg("N"), py::arg("UV"));
