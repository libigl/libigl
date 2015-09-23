m.def("gaussian_curvature", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXd& K
)
{
  return igl::gaussian_curvature(V,F,K);
}, __doc_igl_gaussian_curvature,
py::arg("V"), py::arg("F"), py::arg("K"));
