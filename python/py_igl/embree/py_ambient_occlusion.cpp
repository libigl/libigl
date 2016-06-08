

m.def("ambient_occlusion", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& P,
  const Eigen::MatrixXd& N,
  const int num_samples,
  Eigen::MatrixXd& S
)
{
  return igl::embree::ambient_occlusion(V, F, P, N, num_samples, S);
}, __doc_igl_embree_ambient_occlusion,
py::arg("V"), py::arg("F"), py::arg("P"), py::arg("N"), py::arg("num_samples"), py::arg("S"));

