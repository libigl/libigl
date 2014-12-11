// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

#ifndef IGL_VIEWER_H
#define IGL_VIEWER_H
#ifndef IGL_OPENGL_4
#define IGL_OPENGL_4
#endif

#include <AntTweakBar.h>

#include <vector>
#include <string>
#include <cstdint>

#define IGL_MOD_SHIFT           0x0001
#define IGL_MOD_CONTROL         0x0002
#define IGL_MOD_ALT             0x0004
#define IGL_MOD_SUPER           0x0008

#ifdef ENABLE_XML_SERIALIZATION
  #include <igl/xml/XMLSerializer.h>
  #include <igl/xml/XMLSerialization.h>
#endif

#include <Eigen/Core>
#include <Eigen/Geometry>
#include <igl/viewer/OpenGL_shader.h>
#include <igl/viewer/ViewerData.h>
#include <igl/viewer/OpenGL_state.h>
#include <igl/viewer/ViewerPlugin.h>
#include <igl/viewer/ViewerCore.h>

namespace igl
{
  // GLFW-based mesh viewer
  class Viewer
  {
  public:

    int launch(std::string filename = "");
    void init();

    // Stores all the viewing options
    igl::ViewerCore core;

    // Stores all the data that should be visualized
    igl::ViewerData data;

    // Stores the vbos indices and opengl related settings
    igl::OpenGL_state opengl;

    // List of registered plugins
    std::vector<ViewerPlugin*> plugins;
    void init_plugins();
    void shutdown_plugins();

    // Temporary data stored when the mouse button is pressed
    Eigen::Vector4f down_rotation;
    int current_mouse_x;
    int current_mouse_y;
    int down_mouse_x;
    int down_mouse_y;
    float down_mouse_z;
    Eigen::Vector3f down_translation;
    bool down;
    bool hack_never_moved;

    // Anttweak bar
    TwBar* bar;

    // Keep track of the global position of the scrollwheel
    float scroll_position;

    // UI Enumerations
    enum MouseButton {IGL_LEFT, IGL_MIDDLE, IGL_RIGHT};
    enum MouseMode { NOTHING, ROTATION, ZOOM, PAN, TRANSLATE} mouse_mode;
    enum KeyModifier { NO_KEY = TW_KMOD_NONE, SHIFT = TW_KMOD_SHIFT, CTRL =TW_KMOD_CTRL, ALT = TW_KMOD_ALT } key_modifier;

    Viewer();
    ~Viewer();

    // Mesh IO
    bool load_mesh_from_file(const char* mesh_file_name);
    bool save_mesh_to_file(const char* mesh_file_name);

    // Callbacks
    bool key_down(unsigned char key, int modifier);
    bool key_up(unsigned char key, int modifier);

    bool mouse_down(MouseButton button, int modifier);
    bool mouse_up(MouseButton button, int modifier);

    bool mouse_move(int mouse_x, int mouse_y);
    bool mouse_scroll(float delta_y);

    // Scene IO
    bool load_scene();
    bool save_scene();

    // Draw everything
    void draw();

    // OpenGL context resize
    void resize(int w, int h);


    // C-style callbacks
    bool (*callback_init)(Viewer& viewer);
    bool (*callback_pre_draw)(Viewer& viewer);
    bool (*callback_post_draw)(Viewer& viewer);
    bool (*callback_mouse_down)(Viewer& viewer, int button, int modifier);
    bool (*callback_mouse_up)(Viewer& viewer, int button, int modifier);
    bool (*callback_mouse_move)(Viewer& viewer, int mouse_x, int mouse_y);
    bool (*callback_mouse_scroll)(Viewer& viewer, float delta_y);
    bool (*callback_key_down)(Viewer& viewer, unsigned char key, int modifiers);
    bool (*callback_key_up)(Viewer& viewer, unsigned char key, int modifiers);

    // Pointers to per-callback data
    void* callback_init_data;
    void* callback_pre_draw_data;
    void* callback_post_draw_data;
    void* callback_mouse_down_data;
    void* callback_mouse_up_data;
    void* callback_mouse_move_data;
    void* callback_mouse_scroll_data;
    void* callback_key_down_data;
    void* callback_key_up_data;


    /********* AntTweakBar callbacks *********/
    static void TW_CALL snap_to_canonical_quaternion_cb(void *clientData);
    static void TW_CALL save_scene_cb(void *clientData);
    static void TW_CALL load_scene_cb(void *clientData);
    static void TW_CALL open_dialog_mesh(void *clientData);
    static void TW_CALL align_camera_center_cb(void *clientData);
    static void TW_CALL set_face_based_cb(const void *param, void *clientData);
    static void TW_CALL get_face_based_cb(void *param, void *clientData);
    static void TW_CALL set_invert_normals_cb(const void *param, void *clientData);
    static void TW_CALL get_invert_normals_cb(void *param, void *clientData);
  public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  };

} // end namespace

#ifndef IGL_STATIC_LIBRARY
#  include "Viewer.cpp"
#endif

#endif
