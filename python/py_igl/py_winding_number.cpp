// COMPLETE BINDINGS ========================


m.def("winding_number", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& O,
  Eigen::MatrixXd& W
)
{
  Eigen::VectorXd Wv;
  igl::winding_number(V, F, O, Wv);
  W = Wv;
}, __doc_igl_winding_number,
py::arg("V"), py::arg("F"), py::arg("O"), py::arg("W"));


