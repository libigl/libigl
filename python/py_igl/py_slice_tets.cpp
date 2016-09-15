m.def("slice_tets", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& T,
  const Eigen::MatrixXd& plane,
  Eigen::MatrixXd& U,
  Eigen::MatrixXi& G,
  Eigen::MatrixXi& J,
  Eigen::SparseMatrix<double>& BC
)
{
  assert_is_VectorX("plane", plane);
  Eigen::VectorXd planev;
  if (plane.size() != 0)
    planev = plane;
  Eigen::VectorXi Jv;
  igl::slice_tets(V, T, planev, U, G, Jv, BC);
  J = Jv;
}, __doc_igl_slice_tets,
py::arg("V"), py::arg("T"), py::arg("plane"), py::arg("U"), py::arg("G"), py::arg("J"), py::arg("BC"));

