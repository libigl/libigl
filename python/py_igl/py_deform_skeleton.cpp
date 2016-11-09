// COMPLETE BINDINGS ========================


m.def("deform_skeleton", []
(
  const Eigen::MatrixXd& C,
  const Eigen::MatrixXi& BE,
  const Eigen::MatrixXd& T,
  Eigen::MatrixXd& CT,
  Eigen::MatrixXi& BET
)
{
  return igl::deform_skeleton(C, BE, T, CT, BET);
}, __doc_igl_deform_skeleton,
py::arg("C"), py::arg("BE"), py::arg("T"), py::arg("CT"), py::arg("BET"));





// INCOMPLETE BINDINGS ========================


//m.def("deform_skeleton", []
//(
//  const Eigen::MatrixXd& C,
//  const Eigen::MatrixXi& BE,
//  std::vector<Eigen::Affine3d, Eigen::aligned_allocator<Eigen::Affine3d> > & vA,
//  Eigen::MatrixXd& CT,
//  Eigen::MatrixXi& BET
//)
//{
//  return igl::deform_skeleton(C, BE, vA, CT, BET);
//}, __doc_igl_deform_skeleton,
//py::arg("C"), py::arg("BE"), py::arg("vA"), py::arg("CT"), py::arg("BET"));

