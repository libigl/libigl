// COMPLETE BINDINGS ========================


m.def("hsv_to_rgb", []
(
  const Eigen::MatrixXd& H,
  Eigen::MatrixXd& R
)
{
  return igl::hsv_to_rgb(H, R);
}, __doc_igl_hsv_to_rgb,
py::arg("H"), py::arg("R"));





// INCOMPLETE BINDINGS ========================


//m.def("hsv_to_rgb", []
//(
//  T * hsv,
//  T * rgb
//)
//{
//  return igl::hsv_to_rgb(hsv, rgb);
//}, __doc_igl_hsv_to_rgb,
//py::arg("hsv"), py::arg("rgb"));

//m.def("hsv_to_rgb", []
//(
//  T & h,
//  T & s,
//  T & v,
//  T & r,
//  T & g,
//  T & b
//)
//{
//  return igl::hsv_to_rgb(h, s, v, r, g, b);
//}, __doc_igl_hsv_to_rgb,
//py::arg("h"), py::arg("s"), py::arg("v"), py::arg("r"), py::arg("g"), py::arg("b"));

