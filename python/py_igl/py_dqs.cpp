

m.def("dqs", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXd& W,
  const RotationList& vQ,
  const std::vector<Eigen::MatrixXd> & vT,
  Eigen::MatrixXd& U
)
{
  std::vector<Eigen::Vector3d> vTv;
  for (auto item : vT) {
    assert_is_Vector3("item", item);
    Eigen::Vector3d obj = Eigen::Vector3d(item);
    vTv.push_back(obj);
  }
  return igl::dqs(V, W, vQ, vTv, U);
}, __doc_igl_dqs,
py::arg("V"), py::arg("W"), py::arg("vQ"), py::arg("vT"), py::arg("U"));

