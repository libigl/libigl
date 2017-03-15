py::enum_<igl::EigsType>(m, "EigsType")
    .value("EIGS_TYPE_SM", igl::EIGS_TYPE_SM)
    .value("EIGS_TYPE_LM", igl::EIGS_TYPE_LM)
    .value("NUM_EIGS_TYPES", igl::NUM_EIGS_TYPES)
    .export_values();

m.def("eigs", []
(
  const Eigen::SparseMatrix<double>& A,
  const Eigen::SparseMatrix<double>& B,
  const size_t k,
  const igl::EigsType type,
  Eigen::MatrixXd& sU,
  Eigen::MatrixXd& sS
)
{
  Eigen::VectorXd sSt;
  bool ret = igl::eigs(A,B,k,type,sU,sSt);
  sS = sSt;
  return ret;
}, __doc_igl_eigs,
py::arg("A"), py::arg("B"), py::arg("k"), py::arg("type"), py::arg("sU"), py::arg("sS"));
