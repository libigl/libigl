#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "python.h"
#include <igl/viewer/Viewer.h>
#include <igl/viewer/ViewerCore.h>

void python_export_igl_viewer(py::module &m)
{

  py::module me = m.def_submodule(
    "viewer", "Mesh viewer");

    py::class_<igl::viewer::ViewerData> viewerdata_class(me, "ViewerData");
    viewerdata_class
    .def(py::init<>())
    .def("set_mesh", &igl::viewer::ViewerData::set_mesh)
    .def("clear", &igl::viewer::ViewerData::clear)
    ;

    py::class_<igl::viewer::ViewerCore> viewercore_class(me, "ViewerCore");
    viewercore_class
    .def(py::init<>())
    //.def("align_camera_center", [](igl::viewer::ViewerCore& core, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F){return core.align_camera_center(V,F);})
    .def("align_camera_center",
       (void (igl::viewer::ViewerCore::*) (const Eigen::MatrixXd &, const Eigen::MatrixXi &)) &igl::viewer::ViewerCore::align_camera_center
    )
    .def("init", &igl::viewer::ViewerCore::init)
    ;

    py::class_<igl::viewer::Viewer>(me, "Viewer")
    .def(py::init<>())
    .def_readwrite("data", &igl::viewer::Viewer::data)
    .def_readwrite("core", &igl::viewer::Viewer::core)
    .def("launch", &igl::viewer::Viewer::launch, py::arg("resizable") = true, py::arg("fullscreen") = false)

    // Callbacks
    .def_readwrite("callback_init", &igl::viewer::Viewer::callback_init)
    .def_readwrite("callback_pre_draw", &igl::viewer::Viewer::callback_pre_draw)
    .def_readwrite("callback_post_draw", &igl::viewer::Viewer::callback_post_draw)
    .def_readwrite("callback_mouse_down", &igl::viewer::Viewer::callback_mouse_down)
    .def_readwrite("callback_mouse_up", &igl::viewer::Viewer::callback_mouse_up)
    .def_readwrite("callback_mouse_move", &igl::viewer::Viewer::callback_mouse_move)
    .def_readwrite("callback_mouse_scroll", &igl::viewer::Viewer::callback_mouse_scroll)
    .def_readwrite("callback_key_pressed", &igl::viewer::Viewer::callback_key_pressed)
    .def_readwrite("callback_key_down", &igl::viewer::Viewer::callback_key_down)
    .def_readwrite("callback_key_up", &igl::viewer::Viewer::callback_key_up)
    ;
}
