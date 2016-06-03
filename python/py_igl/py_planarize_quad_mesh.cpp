
m.def("planarize_quad_mesh", []
(
  const Eigen::MatrixXd& Vin,
  const Eigen::MatrixXi& F,
  const int maxIter,
  const double & threshold,
  Eigen::MatrixXd& Vout
)
{
  return igl::planarize_quad_mesh(Vin, F, maxIter, threshold, Vout);
}, __doc_igl_planarize_quad_mesh,
py::arg("Vin"), py::arg("F"), py::arg("maxIter"), py::arg("threshold"), py::arg("Vout"));

