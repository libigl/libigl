m.def("readMESH", []
(
  const std::string mesh_file_name,
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& T,
  Eigen::MatrixXi& F
)
{
  return igl::readMESH(mesh_file_name, V, T, F);
}, __doc_igl_readMESH,
py::arg("mesh_file_name"), py::arg("V"), py::arg("T"), py::arg("F"));


m.def("readMESH", []
(
  const std::string mesh_file_name,
  std::vector<std::vector<double> > & V,
  std::vector<std::vector<int> > & T,
  std::vector<std::vector<int> > & F
)
{
  return igl::readMESH(mesh_file_name, V, T, F);
}, __doc_igl_readMESH,
py::arg("mesh_file_name"), py::arg("V"), py::arg("T"), py::arg("F"));



