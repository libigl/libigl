// Wrap the BBWData class
py::class_<igl::BBWData > BBWData(m, "BBWData");

BBWData
.def(py::init<>())
.def_readwrite("partition_unity", &igl::BBWData::partition_unity)
.def_readwrite("W0", &igl::BBWData::W0)
.def_readwrite("active_set_params", &igl::BBWData::active_set_params)
.def_readwrite("verbosity", &igl::BBWData::verbosity)

.def("print", [](igl::BBWData& data)
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
  igl::BBWData& data,
  Eigen::MatrixXd& W
)
{
  assert_is_VectorX("b",b);
  Eigen::VectorXi bv;
  if (b.size() != 0)
    bv = b;
  return igl::bbw(V, Ele, bv, bc, data, W);
}, __doc_igl_bbw,
py::arg("V"), py::arg("Ele"), py::arg("b"), py::arg("bc"), py::arg("data"), py::arg("W"));
