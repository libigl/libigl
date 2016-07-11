
m.def("directed_edge_orientations", []
(
  const Eigen::MatrixXd& C,
  const Eigen::MatrixXi& E,
  py::list Q
)
{
  std::vector<Eigen::Quaterniond, Eigen::aligned_allocator<Eigen::Quaterniond> > Ql;
  igl::directed_edge_orientations(C, E, Ql);
  for (auto item : Ql) {
    py::object obj = py::cast(item);
    Q.append(obj);
  }
}, __doc_igl_directed_edge_orientations,
py::arg("C"), py::arg("E"), py::arg("Q"));

