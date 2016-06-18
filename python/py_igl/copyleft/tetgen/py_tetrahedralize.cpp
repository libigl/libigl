
m.def("tetrahedralize", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const std::string switches,
  Eigen::MatrixXd& TV,
  Eigen::MatrixXi& TT,
  Eigen::MatrixXi& TF
)
{
  return igl::copyleft::tetgen::tetrahedralize(V, F, switches, TV, TT, TF);
}, __doc_igl_copyleft_tetgen_tetrahedralize,
py::arg("V"), py::arg("F"), py::arg("switches"), py::arg("TV"), py::arg("TT"), py::arg("TF"));


