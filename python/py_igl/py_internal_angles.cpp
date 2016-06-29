

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

m.def("internal_angles", []
(
  const Eigen::MatrixXd& L,
  Eigen::MatrixXd& K
)
{
  return igl::internal_angles(L, K);
}, __doc_igl_internal_angles,
py::arg("L"), py::arg("K"));

