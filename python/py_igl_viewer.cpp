#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "python.h"
#define ENABLE_SERIALIZATION
#include <igl/viewer/Viewer.h>
#include <igl/viewer/ViewerCore.h>
#include <igl/viewer/ViewerData.h>

void python_export_igl_viewer(py::module &m)
{

  py::module me = m.def_submodule(
    "viewer", "Mesh viewer");

/////////////////////// DATA

py::class_<igl::viewer::ViewerData> viewerdata_class(me, "ViewerData");

py::enum_<igl::viewer::ViewerData::DirtyFlags>(viewerdata_class, "DirtyFlags")
    .value("DIRTY_NONE", igl::viewer::ViewerData::DIRTY_NONE)
    .value("DIRTY_POSITION", igl::viewer::ViewerData::DIRTY_POSITION)
    .value("DIRTY_UV", igl::viewer::ViewerData::DIRTY_UV)
    .value("DIRTY_NORMAL", igl::viewer::ViewerData::DIRTY_NORMAL)
    .value("DIRTY_AMBIENT", igl::viewer::ViewerData::DIRTY_AMBIENT)
    .value("DIRTY_DIFFUSE", igl::viewer::ViewerData::DIRTY_DIFFUSE)
    .value("DIRTY_SPECULAR", igl::viewer::ViewerData::DIRTY_SPECULAR)
    .value("DIRTY_TEXTURE", igl::viewer::ViewerData::DIRTY_TEXTURE)
    .value("DIRTY_FACE", igl::viewer::ViewerData::DIRTY_FACE)
    .value("DIRTY_MESH", igl::viewer::ViewerData::DIRTY_MESH)
    .value("DIRTY_OVERLAY_LINES", igl::viewer::ViewerData::DIRTY_OVERLAY_LINES)
    .value("DIRTY_OVERLAY_POINTS", igl::viewer::ViewerData::DIRTY_OVERLAY_POINTS)
    .value("DIRTY_ALL", igl::viewer::ViewerData::DIRTY_ALL)
    .export_values();


    viewerdata_class
    .def(py::init<>())
    .def("set_mesh", &igl::viewer::ViewerData::set_mesh)
    .def("set_colors", &igl::viewer::ViewerData::set_colors)
    .def("clear", &igl::viewer::ViewerData::clear)
    .def("set_face_based", &igl::viewer::ViewerData::set_face_based)

    .def("set_vertices", &igl::viewer::ViewerData::set_vertices)
    .def("set_normals", &igl::viewer::ViewerData::set_normals)

    .def("set_uv",
       (void (igl::viewer::ViewerData::*) (const Eigen::MatrixXd &)) &igl::viewer::ViewerData::set_uv
    )

    .def("set_uv",
       (void (igl::viewer::ViewerData::*) (const Eigen::MatrixXd &, const Eigen::MatrixXi&)) &igl::viewer::ViewerData::set_uv
    )

    .def("set_texture", &igl::viewer::ViewerData::set_texture)
    .def("set_points", &igl::viewer::ViewerData::set_points)
    .def("add_points", &igl::viewer::ViewerData::add_points)
    .def("set_edges", &igl::viewer::ViewerData::set_edges)
    .def("add_edges", &igl::viewer::ViewerData::add_edges)

    .def("add_label", [] (igl::viewer::ViewerData& data, const Eigen::MatrixXd& P,  const std::string& str)
    {
      assert_is_VectorX("P",P);
      data.add_label(P,str);
    })

    .def("compute_normals", &igl::viewer::ViewerData::compute_normals)

    .def("uniform_colors", [] (igl::viewer::ViewerData& data, const Eigen::MatrixXd& ambient, const Eigen::MatrixXd& diffuse, const Eigen::MatrixXd& specular)
    {
      assert_is_Vector3("ambient",ambient);
      assert_is_Vector3("diffuse",diffuse);
      assert_is_Vector3("specular",specular);
      data.uniform_colors(ambient,diffuse, specular);
    })

    .def("grid_texture", &igl::viewer::ViewerData::grid_texture)

    .def_readwrite("V", &igl::viewer::ViewerData::V)
    .def_readwrite("F", &igl::viewer::ViewerData::F)

    .def_readwrite("F_normals", &igl::viewer::ViewerData::F_normals)
    .def_readwrite("F_material_ambient", &igl::viewer::ViewerData::F_material_ambient)
    .def_readwrite("F_material_diffuse", &igl::viewer::ViewerData::F_material_diffuse)
    .def_readwrite("F_material_specular", &igl::viewer::ViewerData::F_material_specular)

    .def_readwrite("V_normals", &igl::viewer::ViewerData::V_normals)
    .def_readwrite("V_material_ambient", &igl::viewer::ViewerData::V_material_ambient)
    .def_readwrite("V_material_diffuse", &igl::viewer::ViewerData::V_material_diffuse)
    .def_readwrite("V_material_specular", &igl::viewer::ViewerData::V_material_specular)

    .def_readwrite("V_uv", &igl::viewer::ViewerData::V_uv)
    .def_readwrite("F_uv", &igl::viewer::ViewerData::F_uv)

    .def_readwrite("texture_R", &igl::viewer::ViewerData::texture_R)
    .def_readwrite("texture_G", &igl::viewer::ViewerData::texture_G)
    .def_readwrite("texture_B", &igl::viewer::ViewerData::texture_B)

    .def_readwrite("lines", &igl::viewer::ViewerData::lines)
    .def_readwrite("points", &igl::viewer::ViewerData::points)
    .def_readwrite("labels_positions", &igl::viewer::ViewerData::labels_positions)
    .def_readwrite("labels_strings", &igl::viewer::ViewerData::labels_strings)
    .def_readwrite("dirty", &igl::viewer::ViewerData::dirty)
    .def_readwrite("face_based", &igl::viewer::ViewerData::face_based)
    .def("serialize", [](igl::viewer::ViewerData& data)
    {
      std::vector<char> a;
      igl::serialize(data,"Data",a);
      return a;
    })

    .def("deserialize", [](igl::viewer::ViewerData& data, const std::vector<char>& a)
    {
      igl::deserialize(data,"Data",a);
      return;
    })

    ;

//////////////////////// CORE

py::class_<igl::viewer::ViewerCore> viewercore_class(me, "ViewerCore");

    py::enum_<igl::viewer::ViewerCore::RotationType>(viewercore_class, "RotationType")
        .value("ROTATION_TYPE_TRACKBALL", igl::viewer::ViewerCore::ROTATION_TYPE_TRACKBALL)
        .value("ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP", igl::viewer::ViewerCore::ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP)
        .value("NUM_ROTATION_TYPES", igl::viewer::ViewerCore::NUM_ROTATION_TYPES)
        .export_values();

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
    [](const igl::viewer::ViewerCore& core) {return Eigen::MatrixXd(core.background_color.cast<double>());},
    [](igl::viewer::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Vector4("background_color",v);
      core.background_color << Eigen::Vector4f(v.cast<float>());
    })

    .def_property("line_color",
    [](const igl::viewer::ViewerCore& core) {return Eigen::MatrixXd(core.line_color.cast<double>());},
    [](igl::viewer::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Vector4("line_color",v);
      core.line_color = Eigen::Vector4f(v.cast<float>());
    })

    .def_property("light_position",
    [](const igl::viewer::ViewerCore& core) {return Eigen::MatrixXd(core.light_position.cast<double>());},
    [](igl::viewer::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Vector3("light_position",v);
      core.light_position = Eigen::Vector3f(v.cast<float>());
    })

    .def_readwrite("lighting_factor",&igl::viewer::ViewerCore::lighting_factor)

    .def_readwrite("model_zoom",&igl::viewer::ViewerCore::model_zoom)

    .def_property("model_translation",
    [](const igl::viewer::ViewerCore& core) {return Eigen::MatrixXd(core.model_translation.cast<double>());},
    [](igl::viewer::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Vector3("model_translation",v);
      core.model_translation = Eigen::Vector3f(v.cast<float>());
    })

    .def_readwrite("model_zoom_uv",&igl::viewer::ViewerCore::model_zoom_uv)

    .def_property("model_translation_uv",
    [](const igl::viewer::ViewerCore& core) {return Eigen::MatrixXd(core.model_translation_uv.cast<double>());},
    [](igl::viewer::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Vector3("model_translation_uv",v);
      core.model_translation_uv = Eigen::Vector3f(v.cast<float>());
    })

    .def_readwrite("camera_zoom",&igl::viewer::ViewerCore::camera_zoom)
    .def_readwrite("orthographic",&igl::viewer::ViewerCore::orthographic)

    .def_property("camera_eye",
    [](const igl::viewer::ViewerCore& core) {return Eigen::MatrixXd(core.camera_eye.cast<double>());},
    [](igl::viewer::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Vector3("camera_eye",v);
      core.camera_eye = Eigen::Vector3f(v.cast<float>());
    })

    .def_property("camera_up",
    [](const igl::viewer::ViewerCore& core) {return Eigen::MatrixXd(core.camera_up.cast<double>());},
    [](igl::viewer::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Vector3("camera_up",v);
      core.camera_up = Eigen::Vector3f(v.cast<float>());
    })

    .def_property("camera_center",
    [](const igl::viewer::ViewerCore& core) {return Eigen::MatrixXd(core.camera_center.cast<double>());},
    [](igl::viewer::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Vector3("camera_center",v);
      core.camera_center = Eigen::Vector3f(v.cast<float>());
    })

    .def_readwrite("camera_view_angle",&igl::viewer::ViewerCore::camera_view_angle)

    .def_readwrite("camera_dnear",&igl::viewer::ViewerCore::camera_dnear)
    .def_readwrite("camera_dfar",&igl::viewer::ViewerCore::camera_dfar)

    .def_readwrite("show_overlay",&igl::viewer::ViewerCore::show_overlay)
    .def_readwrite("show_overlay_depth",&igl::viewer::ViewerCore::show_overlay_depth)
    .def_readwrite("show_texture",&igl::viewer::ViewerCore::show_texture)
    .def_readwrite("show_faces",&igl::viewer::ViewerCore::show_faces)

    .def_readwrite("show_lines",&igl::viewer::ViewerCore::show_lines)
    .def_readwrite("show_vertid",&igl::viewer::ViewerCore::show_vertid)
    .def_readwrite("show_faceid",&igl::viewer::ViewerCore::show_faceid)
    .def_readwrite("invert_normals",&igl::viewer::ViewerCore::invert_normals)
    .def_readwrite("depth_test",&igl::viewer::ViewerCore::depth_test)

    .def_readwrite("point_size",&igl::viewer::ViewerCore::point_size)
    .def_readwrite("line_width",&igl::viewer::ViewerCore::line_width)

    .def_readwrite("is_animating",&igl::viewer::ViewerCore::is_animating)
    .def_readwrite("animation_max_fps",&igl::viewer::ViewerCore::animation_max_fps)

    .def_readwrite("object_scale",&igl::viewer::ViewerCore::object_scale)

    .def_property("viewport",
    [](const igl::viewer::ViewerCore& core) {return Eigen::MatrixXd(core.viewport.cast<double>());},
    [](igl::viewer::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Vector4("viewport",v);
      core.viewport = Eigen::Vector4f(v.cast<float>());
    })

    .def_property("view",
    [](const igl::viewer::ViewerCore& core) {return Eigen::MatrixXd(core.view.cast<double>());},
    [](igl::viewer::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Matrix4("view",v);
      core.view = Eigen::Matrix4f(v.cast<float>());
    })

    .def_property("model",
    [](const igl::viewer::ViewerCore& core) {return Eigen::MatrixXd(core.model.cast<double>());},
    [](igl::viewer::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Matrix4("model",v);
      core.model = Eigen::Matrix4f(v.cast<float>());
    })

    .def_property("proj",
    [](const igl::viewer::ViewerCore& core) {return Eigen::MatrixXd(core.proj.cast<double>());},
    [](igl::viewer::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Matrix4("proj",v);
      core.proj = Eigen::Matrix4f(v.cast<float>());
    })

    .def_readwrite("rotation_type",&igl::viewer::ViewerCore::rotation_type)

    .def("serialize", [](igl::viewer::ViewerCore& core)
    {
      std::vector<char> a;
      igl::serialize(core,"Core",a);
      return a;
    })

    .def("deserialize", [](igl::viewer::ViewerCore& core, const std::vector<char>& a)
    {
      igl::deserialize(core,"Core",a);
      return;
    })

    // TODO: wrap this!
    // Eigen::Quaternionf trackball_angle;
    ;

///////////////////////// VIEWER

// UI Enumerations
    py::class_<igl::viewer::Viewer> viewer_class(me, "Viewer");

    py::enum_<igl::viewer::Viewer::MouseButton>(viewer_class, "MouseButton")
        .value("Left", igl::viewer::Viewer::MouseButton::Left)
        .value("Middle", igl::viewer::Viewer::MouseButton::Middle)
        .value("Right", igl::viewer::Viewer::MouseButton::Right)
        .export_values();

    viewer_class
    .def(py::init<>())
    .def_readwrite("data", &igl::viewer::Viewer::data)
    .def_readwrite("core", &igl::viewer::Viewer::core)
    .def("launch", &igl::viewer::Viewer::launch, py::arg("resizable") = true, py::arg("fullscreen") = false)
    .def("launch_init", &igl::viewer::Viewer::launch_init, py::arg("resizable") = true, py::arg("fullscreen") = false)
    .def("launch_rendering", &igl::viewer::Viewer::launch_rendering, py::arg("loop") = true)
    .def("launch_shut", &igl::viewer::Viewer::launch_shut)
    .def("init", &igl::viewer::Viewer::init)
    .def("serialize", [](igl::viewer::Viewer& viewer)
    {
      std::vector<char> a;
      igl::serialize(viewer.core,"Core",a);
      igl::serialize(viewer.data,"Data",a);

      return a;
    })

    .def("deserialize", [](igl::viewer::Viewer& viewer, const std::vector<char>& a)
    {
      igl::deserialize(viewer.core,"Core",a);
      igl::deserialize(viewer.data,"Data",a);
      return;
    })

    // Scene IO
    .def("load_scene", [](igl::viewer::Viewer& viewer)
    {
      viewer.load_scene();
    })

    .def("load_scene", [](igl::viewer::Viewer& viewer, std::string str)
    {
      viewer.load_scene(str);
    })

    .def("save_scene", &igl::viewer::Viewer::save_scene)

    // Draw everything
    .def("draw", &igl::viewer::Viewer::draw)

    // OpenGL context resize
    .def("resize", &igl::viewer::Viewer::resize)

    // Helper functions
    .def("snap_to_canonical_quaternion", &igl::viewer::Viewer::snap_to_canonical_quaternion)
    .def("open_dialog_load_mesh", &igl::viewer::Viewer::open_dialog_load_mesh)
    .def("open_dialog_save_mesh", &igl::viewer::Viewer::open_dialog_save_mesh)

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
