

m.def("boundary_conditions", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& Ele,
  const Eigen::MatrixXd& C,
  const Eigen::MatrixXi& P,
  const Eigen::MatrixXi& BE,
  const Eigen::MatrixXi& CE,
  Eigen::MatrixXi& b,
  Eigen::MatrixXd& bc
)
{
  assert_is_VectorX("P", P);
  Eigen::VectorXi Pv;
  if (P.size() != 0)
    Pv = P;
  Eigen::VectorXi bv;
  igl::boundary_conditions(V, Ele, C, Pv, BE, CE, bv, bc);
  b = bv;
}, __doc_igl_boundary_conditions,
py::arg("V"), py::arg("Ele"), py::arg("C"), py::arg("P"), py::arg("BE"), py::arg("CE"), py::arg("b"), py::arg("bc"));

