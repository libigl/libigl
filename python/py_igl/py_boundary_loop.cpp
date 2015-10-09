m.def("boundary_loop", []
(
  const Eigen::MatrixXi& F,
  std::vector<std::vector<int> >& L
)
{
  return igl::boundary_loop(F,L);
}, __doc_igl_boundary_loop,
py::arg("F"), py::arg("L"));

m.def("boundary_loop", []
(
  const Eigen::MatrixXi& F,
  std::vector<int>& L
)
{
  return igl::boundary_loop(F,L);
}, __doc_igl_boundary_loop,
py::arg("F"), py::arg("L"));

m.def("boundary_loop", []
(
  const Eigen::MatrixXi& F,
  Eigen::MatrixXi& L
)
{
  Eigen::VectorXi T;
  igl::boundary_loop(F,T);
  L = T;
}, __doc_igl_boundary_loop,
py::arg("F"), py::arg("L"));
