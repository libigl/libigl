// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_VIEWER_VIEWER_H
#define IGL_VIEWER_VIEWER_H

#ifndef IGL_OPENGL_4
#define IGL_OPENGL_4
#endif

#include <vector>
#include <string>
#include <cstdint>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <igl/igl_inline.h>

#include "OpenGL_shader.h"
#include "OpenGL_state.h"
#include "ViewerCore.h"
#include "ViewerData.h"
#include "ViewerPlugin.h"

#define IGL_MOD_SHIFT           0x0001
#define IGL_MOD_CONTROL         0x0002
#define IGL_MOD_ALT             0x0004
#define IGL_MOD_SUPER           0x0008

#ifdef IGL_VIEWER_WITH_NANOGUI
namespace nanogui { class FormHelper; class Screen; class ComboBox; }
#endif

class GLFWwindow;

namespace igl
{
namespace viewer
{
  // GLFW-based mesh viewer
  class Viewer
  {
  public:
    GLFWwindow* window;

    IGL_INLINE int launch(bool resizable = true,bool fullscreen = false,int width = 1200,int height = 768);
    IGL_INLINE int launch_init(bool resizable = true,bool fullscreen = false, int width = 1200, int height = 768);
    IGL_INLINE bool launch_rendering(bool loop = true);
    IGL_INLINE void launch_shut();

    IGL_INLINE void init();

    // Stores windows properties
    bool window_maximized;
    Eigen::Vector2i window_position;
    Eigen::Vector2i window_size;

    // Stores all the viewing options
    ViewerCore core;

    // Stores all the data that should be visualized
    ViewerData data; // current acvtive mesh data
    int active_data_id;
    // Don't access this directly but rather use the functions get_mesh() or set_active_mesh()
    std::vector<ViewerData> data_buffer;
    std::vector<std::string> data_ids;

    // Stores the vbos indices and opengl related settings
    std::vector<OpenGL_state> opengl;

    // List of registered plugins
    std::vector<ViewerPlugin*> plugins;
    IGL_INLINE void init_plugins();
    IGL_INLINE void shutdown_plugins();

    // Temporary data stored when the mouse button is pressed
    Eigen::Quaternionf down_rotation;
    int current_mouse_x;
    int current_mouse_y;
    int down_mouse_x;
    int down_mouse_y;
    float down_mouse_z;
    Eigen::Vector3f down_translation;
    int down_modifier;
    bool down;
    bool hack_never_moved;
    bool key_pressed_handeld;

#ifdef IGL_VIEWER_WITH_NANOGUI
    nanogui::FormHelper* ngui;
    nanogui::Screen* screen;
    nanogui::ComboBox* currentDataCB;
#endif

    // Keep track of the global position of the scrollwheel
    float scroll_position;

    // UI Enumerations
    enum class MouseButton {Left, Middle, Right};
    enum class MouseMode { None, Rotation, Zoom, Pan, Translation} mouse_mode;

    Viewer();
    ~Viewer();

    // Add a new mesh instance
    IGL_INLINE unsigned int add_mesh(const std::string& id);
    IGL_INLINE unsigned int add_mesh(const char* mesh_file_name,const std::string& id);
    // Get the data of a specific mesh instance
    IGL_INLINE ViewerData& get_mesh(unsigned int data_id);
    // Remove the data of a specific mesh instance
    IGL_INLINE bool remove_mesh(unsigned int data_id);
    // Get/Set the active mesh instance exposed as a copy in data
    IGL_INLINE unsigned int get_active_mesh();
    IGL_INLINE bool set_active_mesh(unsigned int data_id);
    IGL_INLINE unsigned int get_mesh_count();

    // Mesh IO
    IGL_INLINE bool load_mesh_from_file(const char* mesh_file_name);
    IGL_INLINE bool load_mesh_from_file(const char* mesh_file_name,unsigned int data_id);
    IGL_INLINE bool save_mesh_to_file(const char* mesh_file_name);
    IGL_INLINE bool save_mesh_to_file(const char* mesh_file_name,unsigned int data_id);

    // Callbacks
    IGL_INLINE bool key_pressed(unsigned int unicode_key,int modifier);
    IGL_INLINE bool key_down(int key,int modifier);
    IGL_INLINE bool key_up(int key,int modifier);

    IGL_INLINE bool mouse_down(MouseButton button,int modifier);
    IGL_INLINE bool mouse_up(MouseButton button,int modifier);

    IGL_INLINE bool mouse_move(int mouse_x,int mouse_y);
    IGL_INLINE bool mouse_scroll(float delta_y);

    // Scene IO
    IGL_INLINE bool load_scene();
    IGL_INLINE bool load_scene(std::string fname);
    IGL_INLINE bool save_scene();

    // Draw everything
    IGL_INLINE void draw();

    // OpenGL context resize
    IGL_INLINE void resize(int w,int h);

    // Helper functions
    IGL_INLINE void snap_to_canonical_quaternion();
    IGL_INLINE void open_dialog_load_mesh();
    IGL_INLINE void open_dialog_save_mesh();

    // C++-style functions
    //
    // Returns **true** if action should be cancelled.
    std::function<bool(Viewer& viewer)> callback_init;
    std::function<bool(Viewer& viewer)> callback_pre_draw;
    std::function<bool(Viewer& viewer)> callback_post_draw;
    std::function<bool(Viewer& viewer, int button, int modifier)> callback_mouse_down;
    std::function<bool(Viewer& viewer, int button, int modifier)> callback_mouse_up;
    std::function<bool(Viewer& viewer, int mouse_x, int mouse_y)> callback_mouse_move;
    std::function<bool(Viewer& viewer, float delta_y)> callback_mouse_scroll;
    std::function<bool(Viewer& viewer, unsigned int key, int modifiers)> callback_key_pressed;
    // THESE SHOULD BE DEPRECATED:
    std::function<bool(Viewer& viewer, unsigned char key, int modifiers)> callback_key_down;
    std::function<bool(Viewer& viewer, unsigned char key, int modifiers)> callback_key_up;

    // Pointers to per-callback data
    void* callback_init_data;
    void* callback_pre_draw_data;
    void* callback_post_draw_data;
    void* callback_mouse_down_data;
    void* callback_mouse_up_data;
    void* callback_mouse_move_data;
    void* callback_mouse_scroll_data;
    void* callback_key_pressed_data;
    void* callback_key_down_data;
    void* callback_key_up_data;

  public:
      EIGEN_MAKE_ALIGNED_OPERATOR_NEW
  };

} // end namespace
} // end namespace

#ifndef IGL_STATIC_LIBRARY
#  include "Viewer.cpp"
#  include "ViewerSerialization.h"
#endif

#endif
