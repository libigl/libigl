m.def("upsample", []
(
  Eigen::MatrixXd& V,
  Eigen::MatrixXi& F
)
{
  return igl::upsample(V, F);
}, __doc_igl_upsample,
py::arg("V"), py::arg("F"));

m.def("upsample", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXd& NV,
  Eigen::MatrixXi& NF
)
{
  return igl::upsample(V, F, NV, NF);
}, __doc_igl_upsample,
py::arg("V"), py::arg("F"), py::arg("NV"), py::arg("NF"));

