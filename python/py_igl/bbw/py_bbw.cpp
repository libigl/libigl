py::enum_<igl::bbw::QPSolver>(m, "QPSolver")
    .value("QP_SOLVER_IGL_ACTIVE_SET", igl::bbw::QP_SOLVER_IGL_ACTIVE_SET)
    .value("QP_SOLVER_MOSEK", igl::bbw::QP_SOLVER_MOSEK)
    .value("NUM_QP_SOLVERS", igl::bbw::NUM_QP_SOLVERS)
    .export_values();

// Wrap the BBWData class
py::class_<igl::bbw::BBWData > BBWData(m, "BBWData");

BBWData
.def(py::init<>())
.def_readwrite("partition_unity", &igl::bbw::BBWData::partition_unity)
.def_readwrite("W0", &igl::bbw::BBWData::W0)
.def_readwrite("active_set_params", &igl::bbw::BBWData::active_set_params)
.def_readwrite("qp_solver", &igl::bbw::BBWData::qp_solver)
.def_readwrite("verbosity", &igl::bbw::BBWData::verbosity)
#ifndef IGL_NO_MOSEK
.def_readwrite("mosek_data", &igl::bbw::BBWData::mosek_data)
#endif
.def("print", [](igl::bbw::BBWData& data)
{
    return data.print();
})
;

m.def("bbw", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& Ele,
  const Eigen::MatrixXi& b,
  const Eigen::MatrixXd& bc,
  igl::bbw::BBWData& data,
  Eigen::MatrixXd& W
)
{
  assert_is_VectorX("b",b);
  Eigen::VectorXi bv;
  if (b.size() != 0)
    bv = b;
  return igl::bbw::bbw(V, Ele, bv, bc, data, W);
}, __doc_igl_bbw_bbw,
py::arg("V"), py::arg("Ele"), py::arg("b"), py::arg("bc"), py::arg("data"), py::arg("W"));
