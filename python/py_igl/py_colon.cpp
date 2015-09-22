m.def("colon", []
(
  const double low,
  const double step,
  const double high,
  Eigen::MatrixXd& I
)
{
  Eigen::Matrix<double,Eigen::Dynamic,1> temp;
  igl::colon<double>(low,step,high,temp);
  I = temp;
}, __doc_igl_colon,
py::arg("low"), py::arg("step"), py::arg("high"), py::arg("I"));

m.def("colon", []
(
  const double low,
  const double high,
  Eigen::MatrixXd& I
)
{
  Eigen::Matrix<double,Eigen::Dynamic,1> temp;
  igl::colon<double>(low,high,temp);
  I = temp;
}, __doc_igl_colon,
py::arg("low"), py::arg("high"), py::arg("I"));

m.def("colon", []
(
  const double& low,
  const double& high
)
{
  return Eigen::MatrixXd(igl::colon<double>(low,high));
}, __doc_igl_colon,
py::arg("low"), py::arg("high"));


m.def("coloni", []
(
  const int low,
  const int step,
  const int high,
  Eigen::MatrixXi& I
)
{
  Eigen::Matrix<int,Eigen::Dynamic,1> temp;
  igl::colon<int>(low,step,high,temp);
  I = temp;
}, __doc_igl_colon,
py::arg("low"), py::arg("step"), py::arg("high"), py::arg("I"));

m.def("coloni", []
(
  const int low,
  const int high,
  Eigen::MatrixXi& I
)
{
  Eigen::Matrix<int,Eigen::Dynamic,1> temp;
  igl::colon<int>(low,high,temp);
  I = temp;
}, __doc_igl_colon,
py::arg("low"), py::arg("high"), py::arg("I"));

m.def("coloni", []
(
  const int& low,
  const int& high
)
{
  return Eigen::MatrixXi(igl::colon<int>(low,high));
}, __doc_igl_colon,
py::arg("low"), py::arg("high"));
