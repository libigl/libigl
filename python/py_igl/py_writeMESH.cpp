m.def("writeMESH", []
(
  const std::string mesh_file_name,
  const std::vector<std::vector<double> > & V,
  const std::vector<std::vector<int> > & T,
  const std::vector<std::vector<int> > & F
)
{
  return igl::writeMESH(mesh_file_name, V, T, F);
}, __doc_igl_writeMESH,
py::arg("mesh_file_name"), py::arg("V"), py::arg("T"), py::arg("F"));

m.def("writeMESH", []
(
  const std::string str,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXd& T,
  const Eigen::MatrixXi& F
)
{
  return igl::writeMESH(str, V, T, F);
}, __doc_igl_writeMESH,
py::arg("str"), py::arg("V"), py::arg("T"), py::arg("F"));

