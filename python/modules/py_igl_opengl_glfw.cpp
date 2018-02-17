// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2017 Sebastian Koch <s.koch@tu-berlin.de> and Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#include <Eigen/Dense>
#include <Eigen/Sparse>

#include "../python_shared.h"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/ViewerCore.h>
#include <igl/opengl/ViewerData.h>
#include <igl/opengl/MeshGL.h>
#include <igl/serialize.h>
#ifdef IGL_VIEWER_WITH_NANOGUI
#include "../../../external/nanogui/include/nanogui/formhelper.h"
#include "../../../external/nanogui/include/nanogui/screen.h"
#endif

void python_export_igl_glfw(py::module &m)
{

  py::module me = m.def_submodule(
    "glfw", "GLFW Mesh viewer");

/////////////////////// DATA

py::class_<igl::opengl::ViewerData> viewerdata_class(me, "ViewerData");

py::enum_<igl::opengl::MeshGL::DirtyFlags>(viewerdata_class, "DirtyFlags")
    .value("DIRTY_NONE", igl::opengl::MeshGL::DIRTY_NONE)
    .value("DIRTY_POSITION", igl::opengl::MeshGL::DIRTY_POSITION)
    .value("DIRTY_UV", igl::opengl::MeshGL::DIRTY_UV)
    .value("DIRTY_NORMAL", igl::opengl::MeshGL::DIRTY_NORMAL)
    .value("DIRTY_AMBIENT", igl::opengl::MeshGL::DIRTY_AMBIENT)
    .value("DIRTY_DIFFUSE", igl::opengl::MeshGL::DIRTY_DIFFUSE)
    .value("DIRTY_SPECULAR", igl::opengl::MeshGL::DIRTY_SPECULAR)
    .value("DIRTY_TEXTURE", igl::opengl::MeshGL::DIRTY_TEXTURE)
    .value("DIRTY_FACE", igl::opengl::MeshGL::DIRTY_FACE)
    .value("DIRTY_MESH", igl::opengl::MeshGL::DIRTY_MESH)
    .value("DIRTY_OVERLAY_LINES", igl::opengl::MeshGL::DIRTY_OVERLAY_LINES)
    .value("DIRTY_OVERLAY_POINTS", igl::opengl::MeshGL::DIRTY_OVERLAY_POINTS)
    .value("DIRTY_ALL", igl::opengl::MeshGL::DIRTY_ALL)
    .export_values();


    viewerdata_class
    .def(py::init<>())
    .def("set_mesh", &igl::opengl::ViewerData::set_mesh)
    .def("set_colors", &igl::opengl::ViewerData::set_colors)
    .def("clear", &igl::opengl::ViewerData::clear)
    .def("set_face_based", &igl::opengl::ViewerData::set_face_based)

    .def("set_vertices", &igl::opengl::ViewerData::set_vertices)
    .def("set_normals", &igl::opengl::ViewerData::set_normals)

    .def("set_uv",
       (void (igl::opengl::ViewerData::*) (const Eigen::MatrixXd &)) &igl::opengl::ViewerData::set_uv
    )

    .def("set_uv",
       (void (igl::opengl::ViewerData::*) (const Eigen::MatrixXd &, const Eigen::MatrixXi&)) &igl::opengl::ViewerData::set_uv
    )

    .def("set_texture",
       (void (igl::opengl::ViewerData::*) (
         const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>&,
         const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>&,
         const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>&)
       ) &igl::opengl::ViewerData::set_texture
    )

    .def("set_texture",
       (void (igl::opengl::ViewerData::*) (
         const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>&,
         const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>&,
         const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>&,
         const Eigen::Matrix<unsigned char,Eigen::Dynamic,Eigen::Dynamic>&)
       ) &igl::opengl::ViewerData::set_texture
    )

    .def("set_points", &igl::opengl::ViewerData::set_points)
    .def("add_points", &igl::opengl::ViewerData::add_points)
    .def("set_edges", &igl::opengl::ViewerData::set_edges)
    .def("add_edges", &igl::opengl::ViewerData::add_edges)

    .def("add_label", [] (igl::opengl::ViewerData& data, const Eigen::MatrixXd& P,  const std::string& str)
    {
      assert_is_VectorX("P",P);
      data.add_label(P,str);
    })

    .def("compute_normals", &igl::opengl::ViewerData::compute_normals)

    .def("uniform_colors", [] (igl::opengl::ViewerData& data, const Eigen::MatrixXd& ambient, const Eigen::MatrixXd& diffuse, const Eigen::MatrixXd& specular)
    {
      if (ambient.cols() == 3)
      {
        assert_is_Vector3("ambient",ambient);
        assert_is_Vector3("diffuse",diffuse);
        assert_is_Vector3("specular",specular);
        Eigen::Vector3d vambient = ambient;
        Eigen::Vector3d vdiffuse = diffuse;
        Eigen::Vector3d vspecular = specular;
        data.uniform_colors(vambient,vdiffuse, vspecular);
      }

      if (ambient.cols() == 4)
      {
        assert_is_Vector4("ambient",ambient);
        assert_is_Vector4("diffuse",diffuse);
        assert_is_Vector4("specular",specular);
        Eigen::Vector4d vambient = ambient;
        Eigen::Vector4d vdiffuse = diffuse;
        Eigen::Vector4d vspecular = specular;
        data.uniform_colors(vambient,vdiffuse,vspecular);
      }

    })

    .def("grid_texture", &igl::opengl::ViewerData::grid_texture)

    .def_readwrite("V", &igl::opengl::ViewerData::V)
    .def_readwrite("F", &igl::opengl::ViewerData::F)

    .def_readwrite("F_normals", &igl::opengl::ViewerData::F_normals)
    .def_readwrite("F_material_ambient", &igl::opengl::ViewerData::F_material_ambient)
    .def_readwrite("F_material_diffuse", &igl::opengl::ViewerData::F_material_diffuse)
    .def_readwrite("F_material_specular", &igl::opengl::ViewerData::F_material_specular)

    .def_readwrite("V_normals", &igl::opengl::ViewerData::V_normals)
    .def_readwrite("V_material_ambient", &igl::opengl::ViewerData::V_material_ambient)
    .def_readwrite("V_material_diffuse", &igl::opengl::ViewerData::V_material_diffuse)
    .def_readwrite("V_material_specular", &igl::opengl::ViewerData::V_material_specular)

    .def_readwrite("V_uv", &igl::opengl::ViewerData::V_uv)
    .def_readwrite("F_uv", &igl::opengl::ViewerData::F_uv)

    .def_readwrite("texture_R", &igl::opengl::ViewerData::texture_R)
    .def_readwrite("texture_G", &igl::opengl::ViewerData::texture_G)
    .def_readwrite("texture_B", &igl::opengl::ViewerData::texture_B)

    .def_readwrite("lines", &igl::opengl::ViewerData::lines)
    .def_readwrite("points", &igl::opengl::ViewerData::points)
    .def_readwrite("labels_positions", &igl::opengl::ViewerData::labels_positions)
    .def_readwrite("labels_strings", &igl::opengl::ViewerData::labels_strings)
    // .def_readwrite("dirty", &igl::opengl::MeshGL::dirty)
    .def_readwrite("face_based", &igl::opengl::ViewerData::face_based)
    .def("serialize", [](igl::opengl::ViewerData& data)
    {
      std::vector<char> a;
      igl::serialize(data,"Data",a);
      return a;
    })

    .def("deserialize", [](igl::opengl::ViewerData& data, const std::vector<char>& a)
    {
      igl::deserialize(data,"Data",a);
      return;
    })

    .def_readwrite("shininess",&igl::opengl::ViewerData::shininess)

    .def_property("line_color",
    [](const igl::opengl::ViewerData& data) {return Eigen::MatrixXd(data.line_color.cast<double>());},
    [](igl::opengl::ViewerData& data, const Eigen::MatrixXd& v)
    {
      assert_is_Vector4("line_color",v);
      data.line_color = Eigen::Vector4f(v.cast<float>());
    })

    .def_readwrite("show_overlay",&igl::opengl::ViewerData::show_overlay)
    .def_readwrite("show_overlay_depth",&igl::opengl::ViewerData::show_overlay_depth)
    .def_readwrite("show_texture",&igl::opengl::ViewerData::show_texture)
    .def_readwrite("show_faces",&igl::opengl::ViewerData::show_faces)

    .def_readwrite("show_lines",&igl::opengl::ViewerData::show_lines)
    .def_readwrite("show_vertid",&igl::opengl::ViewerData::show_vertid)
    .def_readwrite("show_faceid",&igl::opengl::ViewerData::show_faceid)
    .def_readwrite("invert_normals",&igl::opengl::ViewerData::invert_normals)

    .def_readwrite("point_size",&igl::opengl::ViewerData::point_size)
    .def_readwrite("line_width",&igl::opengl::ViewerData::line_width)

    ;

//////////////////////// OPENGL_State

// py::class_<igl::opengl::State> opengl_state_class(me, "OpenGL_state");

//     opengl_state_class
//     .def(py::init<>())
//     .def("init", &igl::opengl::State::init)

//     ;

//////////////////////// CORE

py::class_<igl::opengl::ViewerCore> viewercore_class(me, "ViewerCore");

    py::enum_<igl::opengl::ViewerCore::RotationType>(viewercore_class, "RotationType")
        .value("ROTATION_TYPE_TRACKBALL", igl::opengl::ViewerCore::ROTATION_TYPE_TRACKBALL)
        .value("ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP", igl::opengl::ViewerCore::ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP)
        .value("NUM_ROTATION_TYPES", igl::opengl::ViewerCore::NUM_ROTATION_TYPES)
        .export_values();

    viewercore_class
    .def(py::init<>())
    //.def("align_camera_center", [](igl::opengl::ViewerCore& core, const Eigen::MatrixXd& V, const Eigen::MatrixXi& F){return core.align_camera_center(V,F);})
    .def("init", &igl::opengl::ViewerCore::init)
    .def("shut", &igl::opengl::ViewerCore::shut)
    //.def("InitSerialization", &igl::opengl::ViewerCore::InitSerialization)
    .def("align_camera_center",
       (void (igl::opengl::ViewerCore::*) (const Eigen::MatrixXd &, const Eigen::MatrixXi &)) &igl::opengl::ViewerCore::align_camera_center
    )

    .def("align_camera_center",
       (void (igl::opengl::ViewerCore::*) (const Eigen::MatrixXd &)) &igl::opengl::ViewerCore::align_camera_center
    )

    .def("clear_framebuffers",&igl::opengl::ViewerCore::clear_framebuffers)
    .def("draw",&igl::opengl::ViewerCore::draw)
    .def("draw_buffer",&igl::opengl::ViewerCore::draw_buffer)

    .def_property("background_color",
    [](const igl::opengl::ViewerCore& core) {return Eigen::MatrixXd(core.background_color.cast<double>());},
    [](igl::opengl::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Vector4("background_color",v);
      core.background_color << Eigen::Vector4f(v.cast<float>());
    })

    .def_property("light_position",
    [](const igl::opengl::ViewerCore& core) {return Eigen::MatrixXd(core.light_position.cast<double>());},
    [](igl::opengl::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Vector3("light_position",v);
      core.light_position = Eigen::Vector3f(v.cast<float>());
    })

    .def_readwrite("lighting_factor",&igl::opengl::ViewerCore::lighting_factor)

    .def_property("trackball_angle",
    [](const igl::opengl::ViewerCore& core) {return Eigen::Quaterniond(core.trackball_angle.cast<double>());},
    [](igl::opengl::ViewerCore& core, const Eigen::Quaterniond& q)
    {
      core.trackball_angle = Eigen::Quaternionf(q.cast<float>());
    })

    .def_property("camera_base_translation",
    [](const igl::opengl::ViewerCore& core) {return Eigen::MatrixXd(core.camera_base_translation.cast<double>());},
    [](igl::opengl::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Vector3("camera_base_translation",v);
      core.camera_base_translation = Eigen::Vector3f(v.cast<float>());
    })

    .def_property("camera_translation",
    [](const igl::opengl::ViewerCore& core) {return Eigen::MatrixXd(core.camera_translation.cast<double>());},
    [](igl::opengl::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Vector3("camera_translation",v);
      core.camera_translation = Eigen::Vector3f(v.cast<float>());
    })

    .def_readwrite("camera_base_zoom",&igl::opengl::ViewerCore::camera_base_zoom)
    .def_readwrite("camera_zoom",&igl::opengl::ViewerCore::camera_zoom)
    .def_readwrite("orthographic",&igl::opengl::ViewerCore::orthographic)

    .def_property("camera_eye",
    [](const igl::opengl::ViewerCore& core) {return Eigen::MatrixXd(core.camera_eye.cast<double>());},
    [](igl::opengl::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Vector3("camera_eye",v);
      core.camera_eye = Eigen::Vector3f(v.cast<float>());
    })

    .def_property("camera_up",
    [](const igl::opengl::ViewerCore& core) {return Eigen::MatrixXd(core.camera_up.cast<double>());},
    [](igl::opengl::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Vector3("camera_up",v);
      core.camera_up = Eigen::Vector3f(v.cast<float>());
    })

    .def_property("camera_center",
    [](const igl::opengl::ViewerCore& core) {return Eigen::MatrixXd(core.camera_center.cast<double>());},
    [](igl::opengl::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Vector3("camera_center",v);
      core.camera_center = Eigen::Vector3f(v.cast<float>());
    })

    .def_readwrite("camera_view_angle",&igl::opengl::ViewerCore::camera_view_angle)

    .def_readwrite("camera_dnear",&igl::opengl::ViewerCore::camera_dnear)
    .def_readwrite("camera_dfar",&igl::opengl::ViewerCore::camera_dfar)

    .def_readwrite("depth_test",&igl::opengl::ViewerCore::depth_test)

    .def_readwrite("is_animating",&igl::opengl::ViewerCore::is_animating)
    .def_readwrite("animation_max_fps",&igl::opengl::ViewerCore::animation_max_fps)

    .def_readwrite("object_scale",&igl::opengl::ViewerCore::object_scale)

    .def_property("viewport",
    [](const igl::opengl::ViewerCore& core) {return Eigen::MatrixXd(core.viewport.cast<double>());},
    [](igl::opengl::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Vector4("viewport",v);
      core.viewport = Eigen::Vector4f(v.cast<float>());
    })

    .def_property("view",
    [](const igl::opengl::ViewerCore& core) {return Eigen::MatrixXd(core.view.cast<double>());},
    [](igl::opengl::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Matrix4("view",v);
      core.view = Eigen::Matrix4f(v.cast<float>());
    })

    .def_property("proj",
    [](const igl::opengl::ViewerCore& core) {return Eigen::MatrixXd(core.proj.cast<double>());},
    [](igl::opengl::ViewerCore& core, const Eigen::MatrixXd& v)
    {
      assert_is_Matrix4("proj",v);
      core.proj = Eigen::Matrix4f(v.cast<float>());
    })

    .def_readwrite("rotation_type",&igl::opengl::ViewerCore::rotation_type)

    .def("serialize", [](igl::opengl::ViewerCore& core)
    {
      std::vector<char> a;
      igl::serialize(core,"Core",a);
      return a;
    })

    .def("deserialize", [](igl::opengl::ViewerCore& core, const std::vector<char>& a)
    {
      igl::deserialize(core,"Core",a);
      return;
    })

    // TODO: wrap this!
    // Eigen::Quaternionf trackball_angle;
    ;

///////////////////////// VIEWER

// UI Enumerations
    py::class_<igl::opengl::glfw::Viewer> viewer_class(me, "Viewer");

//    #ifdef IGL_VIEWER_WITH_NANOGUI
//    py::object fh = (py::object) py::module::import("nanogui").attr("FormHelper");
//    py::class_<nanogui::FormHelper> form_helper_class(me, "FormHelper", fh);

//    py::object screen = (py::object) py::module::import("nanogui").attr("Screen");
//    py::class_<nanogui::Screen> screen_class(me, "Screen", screen);
//    #endif

    py::enum_<igl::opengl::glfw::Viewer::MouseButton>(viewer_class, "MouseButton")
        .value("Left", igl::opengl::glfw::Viewer::MouseButton::Left)
        .value("Middle", igl::opengl::glfw::Viewer::MouseButton::Middle)
        .value("Right", igl::opengl::glfw::Viewer::MouseButton::Right)
        .export_values();

    viewer_class
    .def(py::init<>())
    //.def_readwrite("data", &igl::opengl::glfw::Viewer::data)

    // .def_property("data",
    // [](igl::opengl::glfw::Viewer& viewer) {return viewer.data();},
    // [](igl::opengl::glfw::Viewer& viewer, const igl::opengl::ViewerData& data)
    // {
    //   viewer.data() = data;
    // })

    .def("data", &igl::opengl::glfw::Viewer::data,pybind11::return_value_policy::reference)

    .def_readwrite("core", &igl::opengl::glfw::Viewer::core)
    //.def_readwrite("opengl", &igl::opengl::glfw::Viewer::opengl)

    #ifdef IGL_VIEWER_WITH_NANOGUI
    .def_readwrite("ngui", &igl::opengl::glfw::Viewer::ngui)
    .def_readwrite("screen", &igl::opengl::glfw::Viewer::screen)
    #endif

    .def("launch", &igl::opengl::glfw::Viewer::launch, py::arg("resizable") = true, py::arg("fullscreen") = false)
    .def("launch_init", &igl::opengl::glfw::Viewer::launch_init, py::arg("resizable") = true, py::arg("fullscreen") = false)
    .def("launch_rendering", &igl::opengl::glfw::Viewer::launch_rendering, py::arg("loop") = true)
    .def("launch_shut", &igl::opengl::glfw::Viewer::launch_shut)
    .def("init", &igl::opengl::glfw::Viewer::init)
    .def("serialize", [](igl::opengl::glfw::Viewer& viewer)
    {
      std::vector<char> a;
      igl::serialize(viewer.core,"Core",a);
      //igl::serialize(viewer.data,"Data",a); TODO

      return a;
    })

    .def("deserialize", [](igl::opengl::glfw::Viewer& viewer, const std::vector<char>& a)
    {
      igl::deserialize(viewer.core,"Core",a);
      //igl::deserialize(viewer.data,"Data",a);
      return;
    })

    // Scene IO
    .def("load_scene", [](igl::opengl::glfw::Viewer& viewer)
    {
      viewer.load_scene();
    })

    .def("load_scene", [](igl::opengl::glfw::Viewer& viewer, std::string str)
    {
      viewer.load_scene(str);
    })

    .def("save_scene", [](igl::opengl::glfw::Viewer& viewer)
    {
      viewer.save_scene();
    })

    .def("save_scene", [](igl::opengl::glfw::Viewer& viewer, std::string str)
    {
      viewer.save_scene(str);
    })

    // Draw everything
    .def("draw", &igl::opengl::glfw::Viewer::draw)

    // OpenGL context resize
    .def("resize", &igl::opengl::glfw::Viewer::resize)

    // Helper functions
    .def("snap_to_canonical_quaternion", &igl::opengl::glfw::Viewer::snap_to_canonical_quaternion)
    .def("open_dialog_load_mesh", &igl::opengl::glfw::Viewer::open_dialog_load_mesh)
    .def("open_dialog_save_mesh", &igl::opengl::glfw::Viewer::open_dialog_save_mesh)

    // Input handling
    .def_readwrite("current_mouse_x", &igl::opengl::glfw::Viewer::current_mouse_x)
    .def_readwrite("current_mouse_y", &igl::opengl::glfw::Viewer::current_mouse_y)

    // Callbacks
    .def_readwrite("callback_init", &igl::opengl::glfw::Viewer::callback_init)
    .def_readwrite("callback_pre_draw", &igl::opengl::glfw::Viewer::callback_pre_draw)
    .def_readwrite("callback_post_draw", &igl::opengl::glfw::Viewer::callback_post_draw)
    .def_readwrite("callback_mouse_down", &igl::opengl::glfw::Viewer::callback_mouse_down)
    .def_readwrite("callback_mouse_up", &igl::opengl::glfw::Viewer::callback_mouse_up)
    .def_readwrite("callback_mouse_move", &igl::opengl::glfw::Viewer::callback_mouse_move)
    .def_readwrite("callback_mouse_scroll", &igl::opengl::glfw::Viewer::callback_mouse_scroll)
    .def_readwrite("callback_key_pressed", &igl::opengl::glfw::Viewer::callback_key_pressed)
    .def_readwrite("callback_key_down", &igl::opengl::glfw::Viewer::callback_key_down)
    .def_readwrite("callback_key_up", &igl::opengl::glfw::Viewer::callback_key_up)
    ;
}
