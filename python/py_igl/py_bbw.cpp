// Wrap the BBWData class
py::class_<iglData > BBWData(m, "BBWData");

BBWData
.def(py::init<>())
.def_readwrite("partition_unity", &iglData::partition_unity)
.def_readwrite("W0", &iglData::W0)
.def_readwrite("active_set_params", &iglData::active_set_params)
.def_readwrite("verbosity", &iglData::verbosity)
.def("print", [](iglData& data)
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
  iglData& data,
  Eigen::MatrixXd& W
)
{
  assert_is_VectorX("b",b);
  Eigen::VectorXi bv;
  if (b.size() != 0)
    bv = b;
  return igl(V, Ele, bv, bc, data, W);
}, __doc_igl_bbw,
py::arg("V"), py::arg("Ele"), py::arg("b"), py::arg("bc"), py::arg("data"), py::arg("W"));
