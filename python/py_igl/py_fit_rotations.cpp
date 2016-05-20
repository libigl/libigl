m.def("fit_rotations", []
(
  const Eigen::MatrixXd& S,
  const bool single_precision,
  Eigen::MatrixXd& R
)
{
  return igl::fit_rotations(S, single_precision, R);
}, __doc_igl_fit_rotations,
py::arg("S"), py::arg("single_precision"), py::arg("R"));


m.def("fit_rotations_planar", []
(
  const Eigen::MatrixXd& S,
  Eigen::MatrixXd& R
)
{
  return igl::fit_rotations_planar(S, R);
}, __doc_igl_fit_rotations_planar,
py::arg("S"), py::arg("R"));

