m.def("point_mesh_squared_distance", []
(
  const Eigen::MatrixXd& P,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& Ele,
  Eigen::MatrixXd& sqrD,
  Eigen::MatrixXi& I,
  Eigen::MatrixXd& C
)
{
//  assert_is_VectorX("I",I);
  return igl::point_mesh_squared_distance(P, V, Ele, sqrD, I, C);
}, __doc_igl_point_mesh_squared_distance,
py::arg("P"), py::arg("V"), py::arg("Ele"), py::arg("sqrD"), py::arg("I"), py::arg("C"));

