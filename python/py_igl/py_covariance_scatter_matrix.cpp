m.def("covariance_scatter_matrix", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const igl::ARAPEnergyType energy,
  Eigen::SparseMatrix<double>& CSM
)
{
  return igl::covariance_scatter_matrix(V,F,energy,CSM);
}, __doc_igl_covariance_scatter_matrix,
py::arg("V"), py::arg("F"), py::arg("energy"), py::arg("CSM"));
