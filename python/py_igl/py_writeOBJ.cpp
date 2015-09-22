m.def("writeOBJ", []
(
  const std::string str,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& CN,
  const Eigen::MatrixXi& FN,
  const Eigen::MatrixXd& TC,
  const Eigen::MatrixXi& FTC
)
{
  return igl::writeOBJ(str,V,F,CN,FN,TC,FTC);
}, __doc_igl_writeOBJ,
py::arg("str"), py::arg("V"), py::arg("F"), py::arg("CN"), py::arg("FN"), py::arg("TC"), py::arg("FTC"));

m.def("writeOBJ", []
(
  const std::string str,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F
)
{
  return igl::writeOBJ(str,V,F);
}, __doc_igl_writeOBJ,
py::arg("str"), py::arg("V"), py::arg("F"));
