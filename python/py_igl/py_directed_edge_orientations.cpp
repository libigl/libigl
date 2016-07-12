
m.def("directed_edge_orientations", []
(
  const Eigen::MatrixXd& C,
  const Eigen::MatrixXi& E,
  RotationList& Q
)
{
  return igl::directed_edge_orientations(C, E, Q);
}, __doc_igl_directed_edge_orientations,
py::arg("C"), py::arg("E"), py::arg("Q"));

