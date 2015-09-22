m.def("unique", []
(
  const std::vector<double> & A,
  std::vector<double> & C,
  std::vector<size_t> & IA,
  std::vector<size_t> & IC
)
{
  return igl::unique(A,C,IA,IC);
}, __doc_igl_unique,
py::arg("A"), py::arg("C"), py::arg("IA"), py::arg("IC"));

m.def("unique", []
(
  const std::vector<double> & A,
  std::vector<double> & C
)
{
  return igl::unique(A,C);
}, __doc_igl_unique,
py::arg("A"), py::arg("C"));

m.def("unique", []
(
  const Eigen::MatrixXd& A,
  Eigen::MatrixXd& C,
  Eigen::MatrixXi& IA,
  Eigen::MatrixXi& IC
)
{
  return igl::unique(A,C,IA,IC);
}, __doc_igl_unique,
py::arg("A"), py::arg("C"), py::arg("IA"), py::arg("IC"));

m.def("unique", []
(
  const Eigen::MatrixXd& A,
  Eigen::MatrixXd& C
)
{
  return igl::unique(A,C);
}, __doc_igl_unique,
py::arg("A"), py::arg("C"));

m.def("unique_rows", []
(
  const Eigen::MatrixXd& A,
  Eigen::MatrixXd& C,
  Eigen::MatrixXi& IA,
  Eigen::MatrixXi& IC
)
{
  return igl::unique_rows(A,C,IA,IC);
}, __doc_igl_unique,
py::arg("A"), py::arg("C"), py::arg("IA"), py::arg("IC"));



// int


m.def("unique", []
(
  const std::vector<int> & A,
  std::vector<int> & C,
  std::vector<size_t> & IA,
  std::vector<size_t> & IC
)
{
  return igl::unique(A,C,IA,IC);
}, __doc_igl_unique,
py::arg("A"), py::arg("C"), py::arg("IA"), py::arg("IC"));

m.def("unique", []
(
  const std::vector<int> & A,
  std::vector<int> & C
)
{
  return igl::unique(A,C);
}, __doc_igl_unique,
py::arg("A"), py::arg("C"));

m.def("unique", []
(
  const Eigen::MatrixXi& A,
  Eigen::MatrixXi& C,
  Eigen::MatrixXi& IA,
  Eigen::MatrixXi& IC
)
{
  return igl::unique(A,C,IA,IC);
}, __doc_igl_unique,
py::arg("A"), py::arg("C"), py::arg("IA"), py::arg("IC"));

m.def("unique", []
(
  const Eigen::MatrixXi& A,
  Eigen::MatrixXi& C
)
{
  return igl::unique(A,C);
}, __doc_igl_unique,
py::arg("A"), py::arg("C"));

m.def("unique_rows", []
(
  const Eigen::MatrixXi& A,
  Eigen::MatrixXi& C,
  Eigen::MatrixXi& IA,
  Eigen::MatrixXi& IC
)
{
  return igl::unique_rows(A,C,IA,IC);
}, __doc_igl_unique,
py::arg("A"), py::arg("C"), py::arg("IA"), py::arg("IC"));
