m.def("slice_tets", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& T,
  const Eigen::MatrixXd& plane,
  Eigen::MatrixXd& U,
  Eigen::MatrixXi& G,
  Eigen::MatrixXi& J,
  Eigen::SparseMatrix<double> & BC
)
{
  assert_is_VectorX("plane", plane);
  Eigen::VectorXd pl;
  if (plane.size() != 0)
    pl = plane;
  assert_is_VectorX("J", J);
  Eigen::VectorXi Jv;
  if (J.size() != 0)
    Jv = J;
  return igl::slice_tets(V, T, pl, U, G, Jv, BC);
}, __doc_igl_slice_tets,
py::arg("V"), py::arg("T"), py::arg("plane"), py::arg("U"), py::arg("G"), py::arg("J"), py::arg("BC"));

