
m.def("readPNG", []
(
  const std::string png_file,
  Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> & R,
  Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> & G,
  Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> & B,
  Eigen::Matrix<unsigned char, Eigen::Dynamic, Eigen::Dynamic> & A
)
{
  return igl::png::readPNG(png_file, R, G, B, A);
}, __doc_igl_png_readPNG,
py::arg("png_file"), py::arg("R"), py::arg("G"), py::arg("B"), py::arg("A"));

