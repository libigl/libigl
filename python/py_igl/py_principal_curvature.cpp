m.def("principal_curvature", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXd& PD1,
  Eigen::MatrixXd& PD2,
  Eigen::MatrixXd& PV1,
  Eigen::MatrixXd& PV2,
  unsigned radius,
  bool useKring
)
{
  return igl::principal_curvature(V,F,PD1,PD2,PV1,PV2,radius,useKring);
}, __doc_igl_principal_curvature,
py::arg("V"), py::arg("F"), py::arg("PD1"), py::arg("PD2"), py::arg("PV1"), py::arg("PV2"), py::arg("radius") = 5, py::arg("useKring") = true);
