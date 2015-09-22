py::enum_<igl::MassMatrixType>(m, "MassMatrixType")
    .value("MASSMATRIX_TYPE_BARYCENTRIC", igl::MASSMATRIX_TYPE_BARYCENTRIC)
    .value("MASSMATRIX_TYPE_VORONOI", igl::MASSMATRIX_TYPE_VORONOI)
    .value("MASSMATRIX_TYPE_FULL", igl::MASSMATRIX_TYPE_FULL)
    .value("MASSMATRIX_TYPE_DEFAULT", igl::MASSMATRIX_TYPE_DEFAULT)
    .value("NUM_MASSMATRIX_TYPE", igl::NUM_MASSMATRIX_TYPE)
    .export_values();

m.def("massmatrix", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const igl::MassMatrixType type,
  Eigen::SparseMatrix<double>& M
)
{
  return igl::massmatrix(V,F,type,M);
}, __doc_igl_massmatrix,
py::arg("V"), py::arg("F"), py::arg("type"), py::arg("M"));
