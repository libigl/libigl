m.def("randperm", []
(
  int n,
  Eigen::MatrixXi& I
)
{
  return igl::randperm(n, I);
}, __doc_igl_randperm,
py::arg("n"), py::arg("I"));

