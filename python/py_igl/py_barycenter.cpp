m.def("barycenter", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXd& BC
)
{
  return igl::barycenter(V,F,BC);
}, __doc_igl_barycenter,
py::arg("V"), py::arg("F"), py::arg("BC"));
