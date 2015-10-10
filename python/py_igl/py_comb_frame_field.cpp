m.def("comb_frame_field", []
(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  const Eigen::MatrixXd &PD1,
  const Eigen::MatrixXd &PD2,
  const Eigen::MatrixXd &BIS1_combed,
  const Eigen::MatrixXd &BIS2_combed,
  Eigen::MatrixXd &PD1_combed,
  Eigen::MatrixXd &PD2_combed
)
{
  return igl::comb_frame_field(V,F,PD1,PD2,BIS1_combed,BIS2_combed,PD1_combed,PD2_combed);
}, __doc_igl_comb_frame_field,
py::arg("V"), py::arg("F"), py::arg("PD1"), py::arg("PD2"), py::arg("BIS1_combed"), py::arg("BIS2_combed"), py::arg("PD1_combed"), py::arg("PD2_combed"));
