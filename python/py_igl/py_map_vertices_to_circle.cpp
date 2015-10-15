m.def("map_vertices_to_circle", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& bnd,
  Eigen::MatrixXd& UV
)
{
  assert_is_VectorX("bnd",bnd);
  return igl::map_vertices_to_circle(V,bnd,UV);
}, __doc_igl_map_vertices_to_circle,
py::arg("V"), py::arg("bnd"), py::arg("UV"));
