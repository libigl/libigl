

m.def("line_mesh_intersection", []
(
  const Eigen::MatrixXd& V_source,
  const Eigen::MatrixXd& N_source,
  const Eigen::MatrixXd& V_target,
  const Eigen::MatrixXi& F_target
)
{
  return igl::embree::line_mesh_intersection(V_source, N_source, V_target, F_target);
}, __doc_igl_embree_line_mesh_intersection,
py::arg("V_source"), py::arg("N_source"), py::arg("V_target"), py::arg("F_target"));
