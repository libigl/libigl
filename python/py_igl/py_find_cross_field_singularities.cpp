m.def("find_cross_field_singularities", []
(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  const Eigen::MatrixXi &Handle_MMatch,
  Eigen::MatrixXi &isSingularity,
  Eigen::MatrixXi &singularityIndex
)
{
  return igl::find_cross_field_singularities(V,F,Handle_MMatch,isSingularity,singularityIndex);
}, __doc_igl_find_cross_field_singularities,
py::arg("V"), py::arg("F"), py::arg("Handle_MMatch"), py::arg("isSingularity"), py::arg("singularityIndex"));

m.def("find_cross_field_singularities", []
(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  const Eigen::MatrixXd &PD1,
  const Eigen::MatrixXd &PD2,
  Eigen::MatrixXi &isSingularity,
  Eigen::MatrixXi &singularityIndex,
  bool isCombed
)
{
  return igl::find_cross_field_singularities(V,F,PD1,PD2,isSingularity,singularityIndex,isCombed);
}, __doc_igl_find_cross_field_singularities,
py::arg("V"), py::arg("F"), py::arg("PD1"), py::arg("PD2"), py::arg("isSingularity"), py::arg("singularityIndex"),  py::arg("isCombed") = false);
