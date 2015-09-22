
m.def("readDMAT", []
(
  const std::string str,
  Eigen::MatrixXd& W
)
{
  return igl::readDMAT(str,W);
}, __doc_igl_readDMAT,
py::arg("str"), py::arg("W"));
