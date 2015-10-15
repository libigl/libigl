m.def("cut_mesh_from_singularities", []
(
  const Eigen::MatrixXd &V,
  const Eigen::MatrixXi &F,
  const Eigen::MatrixXi &MMatch,
  Eigen::MatrixXi &seams
)
{
  return igl::cut_mesh_from_singularities(V,F,MMatch,seams);
}, __doc_igl_cut_mesh_from_singularities,
py::arg("V"), py::arg("F"), py::arg("MMatch"), py::arg("seams"));
