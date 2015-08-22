m.def("per_corner_normals", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const double corner_threshold,
  Eigen::MatrixXd& CN
)
{
  return igl::per_corner_normals(V,F,corner_threshold,CN);
}, __doc_igl_per_corner_normals,
py::arg("V"), py::arg("F"), py::arg("corner_threshold"), py::arg("CN"));

m.def("per_corner_normals", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& FN,
  const double corner_threshold,
  Eigen::MatrixXd& CN
)
{
  return igl::per_corner_normals(V,F,FN,corner_threshold,CN);
}, __doc_igl_per_corner_normals,
py::arg("V"), py::arg("F"), py::arg("FN"), py::arg("corner_threshold"), py::arg("CN"));

m.def("per_corner_normals", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& FN,
  const double corner_threshold,
  const std::vector<std::vector<int> >& VF,
  Eigen::MatrixXd& CN
)
{
  return igl::per_corner_normals(V,F,FN,VF,corner_threshold,CN);
}, __doc_igl_per_corner_normals,
py::arg("V"), py::arg("F"), py::arg("FN"), py::arg("corner_threshold"), py::arg("VF"), py::arg("CN"));
