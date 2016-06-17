

m.def("writePNG", []
(
  const Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> & R,
  const Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> & G,
  const Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> & B,
  const Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> & A,
  const std::string png_file
)
{
  return igl::png::writePNG(R, G, B, A, png_file);
}, __doc_igl_png_writePNG,
py::arg("R"), py::arg("G"), py::arg("B"), py::arg("A"), py::arg("png_file"));

