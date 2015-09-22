m.def("per_face_normals", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& Z,
  Eigen::MatrixXd& N
)
{
  assert_is_VectorX("Z",Z);
  return igl::per_face_normals(V,F,Z,N);
}, __doc_igl_per_face_normals,
py::arg("V"), py::arg("F"), py::arg("Z"), py::arg("N"));

m.def("per_face_normals", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXd& N
)
{
  return igl::per_face_normals(V,F,N);
}, __doc_igl_per_face_normals,
py::arg("V"), py::arg("F"), py::arg("N"));

m.def("per_face_normals_stable", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXd& N
)
{
  return igl::per_face_normals_stable(V,F,N);
}, __doc_igl_per_face_normals,
py::arg("V"), py::arg("F"), py::arg("N"));
