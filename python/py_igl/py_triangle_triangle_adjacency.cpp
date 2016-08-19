
m.def("triangle_triangle_adjacency", []
(
  const Eigen::MatrixXi& F,
  Eigen::MatrixXi& TT,
  Eigen::MatrixXi& TTi
)
{
  return igl::triangle_triangle_adjacency(F, TT, TTi);
}, __doc_igl_triangle_triangle_adjacency,
py::arg("F"), py::arg("TT"), py::arg("TTi"));

m.def("triangle_triangle_adjacency", []
(
  const Eigen::MatrixXi& F,
  Eigen::MatrixXi& TT
)
{
  return igl::triangle_triangle_adjacency(F, TT);
}, __doc_igl_triangle_triangle_adjacency,
py::arg("F"), py::arg("TT"));
