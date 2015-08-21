#include <Eigen/Dense>

#include "python.h"
#include <igl/readOFF.h>
#include <igl/writeOBJ.h>
#include <igl/per_face_normals.h>

void python_export_igl(py::module &m) {

// readOFF.h

  m.def("readOFF", []
  (
    const std::string str,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F
  )
  {
    return igl::readOFF(str,V,F);
  }, __doc_readOFF,
  py::arg("str"), py::arg("V"), py::arg("F"));

  m.def("readOFF", []
  (
    const std::string str,
    Eigen::MatrixXd& V,
    Eigen::MatrixXi& F,
    Eigen::MatrixXd& N
  )
  {
    return igl::readOFF(str,V,F,N);
  }, __doc_readOFF,
  py::arg("str"), py::arg("V"), py::arg("F"), py::arg("N"));

// writeOBJ.h
m.def("writeOBJ", []
(
  const std::string str,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::MatrixXd& CN,
  const Eigen::MatrixXi& FN,
  const Eigen::MatrixXd& TC,
  const Eigen::MatrixXi& FTC
)
{
  return igl::writeOBJ(str,V,F,CN,FN,TC,FTC);
}, __doc_writeOBJ,
py::arg("str"), py::arg("V"), py::arg("F"), py::arg("CN"), py::arg("FN"), py::arg("TC"), py::arg("FTC"));

m.def("writeOBJ", []
(
  const std::string str,
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F
)
{
  return igl::writeOBJ(str,V,F);
}, __doc_writeOBJ,
py::arg("str"), py::arg("V"), py::arg("F"));

// per_face_normals

m.def("per_face_normals", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  const Eigen::VectorXd& Z,
  Eigen::MatrixXd& N
)
{
  return igl::per_face_normals(V,F,Z,N);
}, __doc_per_face_normals,
py::arg("V"), py::arg("F"), py::arg("Z"), py::arg("N"));

m.def("per_face_normals", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXd& N
)
{
  return igl::per_face_normals(V,F,N);
}, __doc_per_face_normals,
py::arg("V"), py::arg("F"), py::arg("N"));

m.def("per_face_normals_stable", []
(
  const Eigen::MatrixXd& V,
  const Eigen::MatrixXi& F,
  Eigen::MatrixXd& N
)
{
  return igl::per_face_normals_stable(V,F,N);
}, __doc_per_face_normals,
py::arg("V"), py::arg("F"), py::arg("N"));

}
