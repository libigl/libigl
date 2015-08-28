#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "python.h"
#include <igl/viewer/Viewer.h>
#include <igl/viewer/ViewerCore.h>

void python_export_igl_viewer(py::module &m)
{

  py::module me = m.def_submodule(
    "viewer", "Mesh viewer");

/////////////////////// DATA

    py::class_<igl::viewer::ViewerData> viewerdata_class(me, "ViewerData");
    viewerdata_class
    .def(py::init<>())
    .def("set_mesh", &igl::viewer::ViewerData::set_mesh)
    .def("set_colors", &igl::viewer::ViewerData::set_colors)
    .def("clear", &igl::viewer::ViewerData::clear)
    ;


//////////////////////// CORE

    py::class_<igl::viewer::ViewerCore> viewercore_class(me, "ViewerCore");
    viewercore_class
    .def(py::init<>())
    //.def("align_camera_center", [](igl::viewer::ViewerCore& core, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F){return core.align_camera_center(V,F);})
    .def("init", &igl::viewer::ViewerCore::init)
    .def("shut", &igl::viewer::ViewerCore::shut)
    //.def("InitSerialization", &igl::viewer::ViewerCore::InitSerialization)
    .def("align_camera_center",
       (void (igl::viewer::ViewerCore::*) (const Eigen::MatrixXd &, const Eigen::MatrixXi &)) &igl::viewer::ViewerCore::align_camera_center
    )

    .def("align_camera_center",
       (void (igl::viewer::ViewerCore::*) (const Eigen::MatrixXd &)) &igl::viewer::ViewerCore::align_camera_center
    )

    .def("clear_framebuffers",&igl::viewer::ViewerCore::clear_framebuffers)
    .def("draw",&igl::viewer::ViewerCore::draw)
    .def("draw_buffer",&igl::viewer::ViewerCore::draw_buffer)

    .def_readwrite("textrenderer",&igl::viewer::ViewerCore::textrenderer)
    .def_readwrite("shininess",&igl::viewer::ViewerCore::shininess)

    .def_property("background_color",
    [](const igl::viewer::ViewerCore& core) {return Eigen::MatrixXd(core.background_color);},
    [](igl::viewer::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Vector3d("background_color",v);
      core.background_color = Vector3f(v.cast<float>());
    })

    // // Colors
    // Eigen::Vector3f background_color;
    // Eigen::Vector3f line_color;
    //
    // // Lighting
    // Eigen::Vector3f light_position;
    // float lighting_factor;
    //
    // // Trackball angle (quaternion)
    // enum RotationType
    // {
    //   ROTATION_TYPE_TRACKBALL = 0,
    //   ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP = 1,
    //   NUM_ROTATION_TYPES = 2
    // } rotation_type;
    // Eigen::Quaternionf trackball_angle;
    //
    // // Model viewing parameters
    // float model_zoom;
    // Eigen::Vector3f model_translation;
    //
    // // Model viewing paramters (uv coordinates)
    // float model_zoom_uv;
    // Eigen::Vector3f model_translation_uv;
    //
    // // Camera parameters
    // float camera_zoom;
    // bool orthographic;
    // Eigen::Vector3f camera_eye;
    // Eigen::Vector3f camera_up;
    // Eigen::Vector3f camera_center;
    // float camera_view_angle;
    // float camera_dnear;
    // float camera_dfar;
    //
    // // Visualization options
    // bool show_overlay;
    // bool show_overlay_depth;
    // bool show_texture;
    // bool show_faces;
    // bool show_lines;
    // bool show_vertid;
    // bool show_faceid;
    // bool invert_normals;
    // bool depth_test;
    //
    // // Point size / line width
    // float point_size;
    // float line_width;
    //
    // // Animation
    // bool is_animating;
    // double animation_max_fps;
    //
    // // Caches the two-norm between the min/max point of the bounding box
    // float object_scale;
    //
    // // Viewport size
    // Eigen::Vector4f viewport;
    //
    // // Save the OpenGL transformation matrices used for the previous rendering pass
    // Eigen::Matrix4f view;
    // Eigen::Matrix4f model;
    // Eigen::Matrix4f proj;

    ///
    ;

///////////////////////// VIEWER

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
