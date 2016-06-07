

m.def("triangulate", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& E,
  const Eigen::MatrixXd& H,
  const std::string flags,
  Eigen::MatrixXd& V2,
  Eigen::MatrixXi& F2
)
{
  return igl::triangle::triangulate(V, E, H, flags, V2, F2);
}, __doc_igl_triangle_triangulate,
py::arg("V"), py::arg("E"), py::arg("H"), py::arg("flags"), py::arg("V2"), py::arg("F2"));

