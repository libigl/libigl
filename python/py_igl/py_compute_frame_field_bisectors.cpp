m.def("compute_frame_field_bisectors", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& B1,
  const Eigen::MatrixXd& B2,
  const Eigen::MatrixXd& PD1,
  const Eigen::MatrixXd& PD2,
  Eigen::MatrixXd& BIS1,
  Eigen::MatrixXd& BIS2
)
{
  return igl::compute_frame_field_bisectors(V,F,B1,B2,PD1,PD2,BIS1,BIS2);
}, __doc_igl_compute_frame_field_bisectors,
py::arg("V"), py::arg("F"), py::arg("B1"), py::arg("B2"), py::arg("PD1"), py::arg("PD2"), py::arg("BIS1"), py::arg("BIS2"));

m.def("compute_frame_field_bisectors", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& PD1,
  const Eigen::MatrixXd& PD2,
  Eigen::MatrixXd& BIS1,
  Eigen::MatrixXd& BIS2
)
{
  return igl::compute_frame_field_bisectors(V,F,PD1,PD2,BIS1,BIS2);
}, __doc_igl_compute_frame_field_bisectors,
py::arg("V"), py::arg("F"), py::arg("PD1"), py::arg("PD2"), py::arg("BIS1"), py::arg("BIS2"));
