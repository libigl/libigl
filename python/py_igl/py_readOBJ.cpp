m.def("readOBJ", []
(
  const std::string str,
  Eigen::MatrixXd& V,
  Eigen::MatrixXd& TC,
  Eigen::MatrixXd& CN,
  Eigen::MatrixXi& F,
  Eigen::MatrixXi& FTC,
  Eigen::MatrixXi& FN
)
{
  return igl::readOBJ(str,V,TC,CN,F,FTC,FN);
}, __doc_igl_readOBJ,
py::arg("str"), py::arg("V"), py::arg("TC"), py::arg("CN"), py::arg("F"), py::arg("FTC"), py::arg("FN"));

m.def("readOBJ", []
(
  const std::string str,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F
)
{
  return igl::readOBJ(str,V,F);
}, __doc_igl_readOBJ,
py::arg("str"), py::arg("V"), py::arg("F"));
