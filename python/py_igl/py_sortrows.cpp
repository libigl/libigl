m.def("sortrows", []
(
  const Eigen::MatrixXd& X,
  const bool ascending,
  Eigen::MatrixXd& Y,
  Eigen::MatrixXi& I
)
{
  return igl::sortrows(X,ascending,Y,I);
}, __doc_igl_sortrows,
py::arg("X"), py::arg("ascending"), py::arg("Y"), py::arg("I"));

m.def("sortrows", []
(
  const Eigen::MatrixXi& X,
  const bool ascending,
  Eigen::MatrixXi& Y,
  Eigen::MatrixXi& I
)
{
  return igl::sortrows(X,ascending,Y,I);
}, __doc_igl_sortrows,
py::arg("X"), py::arg("ascending"), py::arg("Y"), py::arg("I"));
