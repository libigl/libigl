

m.def("internal_angles", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXd& K
)
{
  return igl::internal_angles(V, F, K);
}, __doc_igl_internal_angles,
py::arg("V"), py::arg("F"), py::arg("K"));

m.def("internal_angles_using_squared_edge_lengths", []
(
  const Eigen::MatrixXd& L_sq,
  Eigen::MatrixXd& K
)
{
  return igl::internal_angles_using_squared_edge_lengths(L_sq, K);
}, __doc_igl_internal_angles,
py::arg("L_sq"), py::arg("K"));

m.def("internal_angles_using_edge_lengths", []
(
  const Eigen::MatrixXd& L,
  Eigen::MatrixXd& K
)
{
  return igl::internal_angles_using_edge_lengths(L, K);
}, __doc_igl_internal_angles,
py::arg("L"), py::arg("K"));

