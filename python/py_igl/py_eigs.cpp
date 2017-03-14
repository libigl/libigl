py::enum_<igl::EigsType>(m, "EigsType")
    .value("EIGS_TYPE_SM", igl::EIGS_TYPE_SM)
    .value("EIGS_TYPE_LM", igl::EIGS_TYPE_LM)
    .value("NUM_EIGS_TYPES", igl::NUM_EIGS_TYPES)
    .export_values();

m.def("eigs", []
(
  const Eigen::SparseMatrix<double>& A,
  const Eigen::SparseMatrix<double>& B,
  const igl::EigsType type,
  Eigen::MatrixXd& sU,
  Eigen::MatrixXd& sS,
  const size_t k,
  unsigned int max_iter
)
{
  Eigen::VectorXd sSt;
  bool ret = igl::eigs(A,B,type,sU,sSt,k,max_iter);
  sS = sSt;
  return ret;
}, __doc_igl_eigs,
py::arg("A"), py::arg("B"), py::arg("type"), py::arg("sU"), py::arg("sS"), py::arg("k"), py::arg("max_iter"));
