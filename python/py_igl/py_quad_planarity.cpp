

m.def("quad_planarity", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXd& P
)
{
  Eigen::VectorXd Pv;
  igl::quad_planarity(V, F, Pv);
  P = Pv;
}, __doc_igl_quad_planarity,
py::arg("V"), py::arg("F"), py::arg("P"));

