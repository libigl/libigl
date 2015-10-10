m.def("rotate_vectors", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXd& A,
  const Eigen::MatrixXd& B1,
  const Eigen::MatrixXd& B2
)
{
  assert_is_VectorX("A",A);
  return igl::rotate_vectors(V,A,B1,B2);
}, __doc_igl_rotate_vectors,
py::arg("V"), py::arg("A"), py::arg("B1"), py::arg("B2"));
