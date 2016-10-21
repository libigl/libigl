m.def("column_to_quats", []
(
  const Eigen::MatrixXd& Q,
  RotationList& vQ
)
{
  return igl::column_to_quats(Q, vQ);
}, __doc_igl_column_to_quats,
py::arg("Q"), py::arg("vQ"));

