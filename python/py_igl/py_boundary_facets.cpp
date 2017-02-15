m.def("boundary_facets", []
(
  const Eigen::MatrixXi& T,
  Eigen::MatrixXi& F
)
{
  return igl::boundary_facets(T,F);
}, __doc_igl_boundary_facets,
py::arg("T"), py::arg("F"));

m.def("boundary_facets", []
(
  const Eigen::MatrixXi& T
)
{
  Eigen::MatrixXi F;
  igl::boundary_facets(T,F);
  return F;
}, __doc_igl_boundary_facets,
py::arg("T"));

m.def("boundary_facets", []
(
  const std::vector<std::vector<int> > & T,
  std::vector<std::vector<int> > & F
)
{
  return igl::boundary_facets(T,F);
}, __doc_igl_boundary_facets,
py::arg("T"), py::arg("F"));
