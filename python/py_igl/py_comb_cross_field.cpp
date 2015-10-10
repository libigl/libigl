m.def("comb_cross_field", []
(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  const Eigen::MatrixXd &PD1in,
  const Eigen::MatrixXd &PD2in,
  Eigen::MatrixXd &PD1out,
  Eigen::MatrixXd &PD2out
)
{
  return igl::comb_cross_field(V,F,PD1in,PD2in,PD1out,PD2out);
}, __doc_igl_comb_cross_field,
py::arg("V"), py::arg("F"), py::arg("PD1in"), py::arg("PD2in"), py::arg("PD1out"), py::arg("PD2out"));
