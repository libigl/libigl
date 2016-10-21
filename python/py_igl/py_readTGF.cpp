// COMPLETE BINDINGS ========================


m.def("readTGF", []
(
  const std::string tgf_filename,
  Eigen::MatrixXd& C,
  Eigen::MatrixXi& E,
  Eigen::MatrixXi& P,
  Eigen::MatrixXi& BE,
  Eigen::MatrixXi& CE,
  Eigen::MatrixXi& PE
)
{
  Eigen::VectorXi Pv;
  bool ret = igl::readTGF(tgf_filename, C, E, Pv, BE, CE, PE);
  P = Pv;
  return ret;
}, __doc_igl_readTGF,
py::arg("tgf_filename"), py::arg("C"), py::arg("E"), py::arg("P"), py::arg("BE"), py::arg("CE"), py::arg("PE"));

m.def("readTGF", []
(
  const std::string tgf_filename,
  Eigen::MatrixXd& C,
  Eigen::MatrixXi& E
)
{
  return igl::readTGF(tgf_filename, C, E);
}, __doc_igl_readTGF,
py::arg("tgf_filename"), py::arg("C"), py::arg("E"));





// INCOMPLETE BINDINGS ========================


//m.def("readTGF", []
//(
//  const std::string tgf_filename,
//  std::vector<std::vector<double> > & C,
//  std::vector<std::vector<int> > & E,
//  std::vector<int> & P,
//  std::vector<std::vector<int> > & BE,
//  std::vector<std::vector<int> > & CE,
//  std::vector<std::vector<int> > & PE
//)
//{
//  return igl::readTGF(tgf_filename, C, E, P, BE, CE, PE);
//}, __doc_igl_readTGF,
//py::arg("tgf_filename"), py::arg("C"), py::arg("E"), py::arg("P"), py::arg("BE"), py::arg("CE"), py::arg("PE"));

