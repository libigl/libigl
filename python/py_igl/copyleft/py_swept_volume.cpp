
m.def("swept_volume", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const std::function<Eigen::Affine3d (const double)> & transform,
  const size_t steps,
  const size_t grid_res,
  const size_t isolevel,
  Eigen::MatrixXd& SV,
  Eigen::MatrixXi& SF
)
{
  return igl::copyleft::swept_volume(V, F, transform, steps, grid_res, isolevel, SV, SF);
}, __doc_igl_copyleft_swept_volume,
py::arg("V"), py::arg("F"), py::arg("transform"), py::arg("steps"), py::arg("grid_res"), py::arg("isolevel"), py::arg("SV"), py::arg("SF"));



