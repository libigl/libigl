m.def("n_polyvector", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXi& b,
  const Eigen::MatrixXd& bc,
  Eigen::MatrixXd &output
)
{
  assert_is_VectorX("b",b);

  Eigen::VectorXi bt;
  if (b.size() != 0)
    bt = b;

  igl::n_polyvector(V,F,bt,bc,output);

}, __doc_igl_n_polyvector,
py::arg("V"), py::arg("F"), py::arg("b"), py::arg("bc"), py::arg("output"));
