#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "python.h"
#include <igl/viewer/Viewer.h>

void python_export_igl_viewer(py::module &m)
{
  py::module me = m.def_submodule(
    "viewer", "Mesh viewer");

    py::class_<igl::viewer::ViewerData> viewerdata_class(me, "ViewerData");
    viewerdata_class
    .def(py::init<>())
    .def("set_mesh", &igl::viewer::ViewerData::set_mesh)
    ;

    py::class_<igl::viewer::ViewerCore> viewercore_class(me, "ViewerCore");
    viewercore_class
    .def(py::init<>())
    .def("align_camera_center", [](igl::viewer::ViewerCore& core, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F){return core.align_camera_center(V,F);})
    //.def("align_camera_center", &igl::viewer::ViewerCore::align_camera_center)
    .def("init", &igl::viewer::ViewerCore::init)
    ;

    py::class_<igl::viewer::Viewer> viewer_class(me, "Viewer");
    viewer_class
    .def(py::init<>())
    .def_readwrite("data", &igl::viewer::Viewer::data)
    .def_readwrite("core", &igl::viewer::Viewer::core)
    .def("launch", &igl::viewer::Viewer::launch, py::arg("resizable") = true, py::arg("fullscreen") = false)
    ;
}
