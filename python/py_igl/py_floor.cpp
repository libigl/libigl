m.def("floor", []
(
  const Eigen::MatrixXd& X,
  Eigen::MatrixXi& Y
)
{
  return igl::floor(X,Y);
}, __doc_igl_floor,
py::arg("X"), py::arg("Y"));
