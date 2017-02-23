// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

// Must defined this before including Viewer.h
#define IGL_VIEWER_VIEWER_CPP
#include "Viewer.h"

#ifdef _WIN32
#  include <windows.h>
#  undef max
#  undef min
#endif

#include <chrono>
#include <thread>

#ifndef __APPLE__
#  define GLEW_STATIC
#  include <GL/glew.h>
#endif

#ifdef __APPLE__
#   include <OpenGL/gl3.h>
#   define __gl_h_ /* Prevent inclusion of the old gl.h */
#else
#   include <GL/gl.h>
#endif

#include <Eigen/LU>

//#define GLFW_INCLUDE_GLU
#if defined(__APPLE__)
#define GLFW_INCLUDE_GLCOREARB
#else
#define GL_GLEXT_PROTOTYPES
#endif

#include <GLFW/glfw3.h>

#include <cmath>
#include <cstdio>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>
#include <algorithm>
#include <limits>
#include <cassert>

#ifdef IGL_VIEWER_WITH_NANOGUI
#  include <nanogui/formhelper.h>
#  include <nanogui/screen.h>
#endif

#include <igl/project.h>
#include <igl/get_seconds.h>
#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/adjacency_list.h>
#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>
#include <igl/massmatrix.h>
#include <igl/file_dialog_open.h>
#include <igl/file_dialog_save.h>
#include <igl/quat_mult.h>
#include <igl/axis_angle_to_quat.h>
#include <igl/trackball.h>
#include <igl/two_axis_valuator_fixed_up.h>
#include <igl/snap_to_canonical_view_quat.h>
#include <igl/unproject.h>

#ifdef IGL_VIEWER_WITH_NANOGUI_SERIALIZATION
#include <igl/serialize.h>
#endif

// Internal global variables used for glfw event handling
static igl::viewer::Viewer * __viewer;
static double highdpi = 1;
static double scroll_x = 0;
static double scroll_y = 0;

static void glfw_mouse_press(GLFWwindow* window, int button, int action, int modifier)
{
  bool tw_used =
#ifdef IGL_VIEWER_WITH_NANOGUI
    __viewer->screen->mouseButtonCallbackEvent(button,action,modifier);
#else
    false;
#endif

  igl::viewer::Viewer::MouseButton mb;

  if (button == GLFW_MOUSE_BUTTON_1)
    mb = igl::viewer::Viewer::MouseButton::Left;
  else if (button == GLFW_MOUSE_BUTTON_2)
    mb = igl::viewer::Viewer::MouseButton::Right;
  else //if (button == GLFW_MOUSE_BUTTON_3)
    mb = igl::viewer::Viewer::MouseButton::Middle;

  if (action == GLFW_PRESS)
  {
    if(!tw_used)
    {
      __viewer->mouse_down(mb,modifier);
    }
  } else
  {
    // Always call mouse_up on up
    __viewer->mouse_up(mb,modifier);
  }

}

static void glfw_error_callback(int error, const char* description)
{
  fputs(description, stderr);
}

static void glfw_char_mods_callback(GLFWwindow* window, unsigned int codepoint, int modifier)
{
  // TODO: pass to nanogui (although it's also using physical key down/up
  // rather than character codes...
#ifdef IGL_VIEWER_WITH_NANOGUI
  if(! __viewer->screen->charCallbackEvent(codepoint) )
#endif
  {
    __viewer->key_pressed(codepoint, modifier);
  }
}

static void glfw_key_callback(GLFWwindow* window, int key, int scancode, int action, int modifier)
{
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    glfwSetWindowShouldClose(window, GL_TRUE);

#ifdef IGL_VIEWER_WITH_NANOGUI
  if (__viewer->screen->keyCallbackEvent(key,scancode,action,modifier) == false)
#endif
  {
    if (action == GLFW_PRESS || action == GLFW_REPEAT)
      __viewer->key_down(key, modifier);
    else if(action == GLFW_RELEASE)
      __viewer->key_up(key, modifier);
  }
}

static void glfw_window_position(GLFWwindow* window,int left,int top)
{
#ifdef GLFW_VERSION_HIGHER_THAN_3_1
  __viewer->window_maximized = glfwGetWindowAttrib(window,GLFW_MAXIMIZED);
#endif

  if(!__viewer->window_maximized)
  {
    __viewer->window_position(0) = left;
    __viewer->window_position(1) = top;
  }
}

static void glfw_window_size(GLFWwindow* window, int width, int height)
{
  int w = width*highdpi;
  int h = height*highdpi;

  __viewer->resize(w, h);

  // TODO: repositioning of the nanogui
}

static void glfw_mouse_move(GLFWwindow* window, double x, double y)
{
  if(
#ifdef IGL_VIEWER_WITH_NANOGUI
      __viewer->screen->cursorPosCallbackEvent(x,y) == false &&
#endif
      true
    )
  {
    __viewer->mouse_move(x*highdpi, y*highdpi);
  }
}

static void glfw_mouse_scroll(GLFWwindow* window, double x, double y)
{
  using namespace std;
  scroll_x += x;
  scroll_y += y;

#ifdef IGL_VIEWER_WITH_NANOGUI
  if (__viewer->screen->scrollCallbackEvent(x,y) == false)
#endif
  {
    __viewer->mouse_scroll(y);
  }
}

static void glfw_drop_callback(GLFWwindow *window,int count,const char **filenames)
{
#ifdef IGL_VIEWER_WITH_NANOGUI
  __viewer->screen->dropCallbackEvent(count,filenames);
#endif
}

#ifdef IGL_STATIC_LIBRARY
#  include "ViewerSerialization.h"
#endif

namespace igl
{
namespace viewer
{
  IGL_INLINE void Viewer::init()
  {
#ifdef IGL_VIEWER_WITH_NANOGUI
    using namespace nanogui;

    ngui->setFixedSize(Eigen::Vector2i(60,20));

    // Create nanogui widgets
    nanogui::Window *window = ngui->addWindow(Eigen::Vector2i(10,10),"libIGL-Viewer");

    // ---------------------- LOADING ----------------------

#ifdef IGL_VIEWER_WITH_NANOGUI_SERIALIZATION
    {
      ngui->addGroup("Workspace");

      Widget* container = new Widget(window);
      ngui->addWidget("",container);

      GridLayout* layout = new GridLayout(Orientation::Horizontal,2,Alignment::Fill);
      container->setLayout(layout);

      Button* loadButton = new Button(container,"Load");
      loadButton->setFixedHeight(25);
      loadButton->setCallback([&]() {this->load_scene();});

      Button* saveButton = new Button(container,"Save");
      saveButton->setFixedHeight(25);
      saveButton->setCallback([&]() {this->save_scene();});
    }
#endif

#ifdef IGL_VIEWER_WITH_NANOGUI_IO
    {
      ngui->addGroup("Mesh");

      Widget* container = new Widget(window);
      ngui->addWidget("",container);

      GridLayout* layout = new GridLayout(Orientation::Horizontal,2,Alignment::Fill);
      container->setLayout(layout);

      Button* loadButton = new Button(container,"Load");
      loadButton->setFixedHeight(25);
      loadButton->setCallback([&]() {this->open_dialog_load_mesh();});

      Button* saveButton = new Button(container,"Save");
      saveButton->setFixedHeight(25);
      saveButton->setCallback([&]() {this->open_dialog_save_mesh();});
    }
#endif

    ngui->addGroup("Viewing Options");

    ngui->addButton("Center object",[&]() {this->core.align_camera_center(this->data);});
    ngui->addButton("Canonical view",[&]()
    {
      this->snap_to_canonical_quaternion();
    });
    ngui->addVariable("Zoom",core.camera_zoom);
    ngui->addVariable("Orthographic view",core.orthographic);

    ngui->addVariable("Background",(nanogui::Color &) core.background_color);
    ngui->addVariable("Line color",(nanogui::Color &) core.line_color);
    ngui->addVariable("Shininess",core.shininess);

    ngui->addGroup("Draw options");

#ifdef IGL_VIEWER_WITH_NANOGUI_MULTIMESH
    currentDataCB = new ComboBox(window,data_ids);
    currentDataCB->setCallback([&](int id){
      if(id != active_data_id)
      {
        currentDataCB->setSelectedIndex(id);
        set_active_mesh(id);
      }
    });
    currentDataCB->setFixedHeight(20);
    currentDataCB->setFontSize(16);
    currentDataCB->setSelectedIndex(active_data_id);
    ngui->addWidget("Active Mesh",currentDataCB);
#endif

    ngui->addVariable<bool>("Visible",[&](bool checked)
    {
      data.visible = checked;
    },[&]()
    {
      return data.visible;
    });

    ngui->addVariable<bool>("Face-based",[&](bool checked)
    {
      data.set_face_based(checked);
    },[&]()
    {
      return data.face_based;
    });

    ngui->addVariable<bool>("Show texture",[&](bool checked) {
      data.show_texture = checked;
    },[&]() {
      return data.show_texture;
    });

    ngui->addVariable<bool>("Invert normals",[&](bool checked){
      data.dirty |= ViewerData::DIRTY_NORMAL;
      data.invert_normals = checked;
    },[&](){
      return data.invert_normals;
    });

    ngui->addVariable<bool>("Wireframe",[&](bool checked) {
      data.show_lines = checked;
    },[&]() {
      return data.show_lines;
    });

    ngui->addVariable<bool>("Fill",[&](bool checked) {
      data.show_faces = checked;
    },[&]() {
      return data.show_faces;
    });

    ngui->addVariable<bool>("Show overlay",[&](bool checked) {
      data.show_overlay = checked;
    },[&]() {
      return data.show_overlay;
    });

    ngui->addVariable<bool>("Show overlay depth",[&](bool checked) {
      data.show_overlay_depth = checked;
    },[&]() {
      return data.show_overlay_depth;
    });

    ngui->addVariable<bool>("Show vertex labels",[&](bool checked) {
      data.show_vertid = checked;
    },[&]() {
      return data.show_vertid;
    });

    ngui->addVariable<bool>("Show faces labels",[&](bool checked) {
      data.show_faceid = checked;
    },[&]() {
      return data.show_faceid;
    });

    screen->setVisible(true);
    screen->performLayout();

#endif

    core.init();

    if (callback_init)
      if (callback_init(*this))
        return;

    init_plugins();
  }

  IGL_INLINE Viewer::Viewer()
  {
#ifdef IGL_VIEWER_WITH_NANOGUI
    ngui = nullptr;
    screen = nullptr;
    currentDataCB = nullptr;
#endif

    __viewer = nullptr;
    window = nullptr;

    active_data_id = 0;
    data_buffer.push_back(ViewerData());
    data_ids.push_back("mesh");
    data = data_buffer[active_data_id];
    opengl.push_back(OpenGL_state());

    // Temporary variables initialization
    down = false;
    hack_never_moved = true;
    scroll_position = 0.0f;

    // Per face
    for(auto d : data_buffer)
      d.set_face_based(false);

    // C-style callbacks
    callback_init         = nullptr;
    callback_pre_draw     = nullptr;
    callback_post_draw    = nullptr;
    callback_mouse_down   = nullptr;
    callback_mouse_up     = nullptr;
    callback_mouse_move   = nullptr;
    callback_mouse_scroll = nullptr;
    callback_key_down     = nullptr;
    callback_key_up       = nullptr;

    callback_init_data          = nullptr;
    callback_pre_draw_data      = nullptr;
    callback_post_draw_data     = nullptr;
    callback_mouse_down_data    = nullptr;
    callback_mouse_up_data      = nullptr;
    callback_mouse_move_data    = nullptr;
    callback_mouse_scroll_data  = nullptr;
    callback_key_down_data      = nullptr;
    callback_key_up_data        = nullptr;

#ifndef IGL_VIEWER_VIEWER_QUIET
    const std::string usage(R"(igl::viewer::Viewer usage:
  [drag]  Rotate scene
  A,a     Toggle animation (tight draw loop)
  F,f     Toggle face based
  I,i     Toggle invert normals
  L,l     Toggle wireframe
  O,o     Toggle orthographic/perspective projection
  T,t     Toggle filled faces
  Z       Snap to canonical view
  [,]     Toggle between rotation control types (e.g. trackball, two-axis
          valuator with fixed up))"
#ifdef IGL_VIEWER_WITH_NANOGUI
		R"(
  ;       Toggle vertex labels
  :       Toggle face labels)"
#endif
);
    std::cout<<usage<<std::endl;
#endif
  }

  IGL_INLINE void Viewer::init_plugins()
  {
    // Init all plugins
    for (unsigned int i = 0; i<plugins.size(); ++i)
    {
      plugins[i]->init(this);
    }
  }

  IGL_INLINE Viewer::~Viewer()
  {
  }

  IGL_INLINE void Viewer::shutdown_plugins()
  {
    for (unsigned int i = 0; i<plugins.size(); ++i)
    {
      plugins[i]->shutdown();
    }
  }

  IGL_INLINE unsigned int Viewer::add_mesh(const std::string& id)
  {
    data_buffer.push_back(ViewerData());
    
    data_ids.push_back(id);
    
    opengl.push_back(OpenGL_state());

    if(window != nullptr)
    {
      opengl[opengl.size()-1].init();
#if defined(IGL_VIEWER_WITH_NANOGUI) && defined(IGL_VIEWER_WITH_NANOGUI_MULTIMESH)
      currentDataCB->setItems(data_ids);
      screen->performLayout();
#endif
    }
    else
    {
      set_active_mesh(opengl.size()-1);
    }

    return opengl.size()-1;
  }

  IGL_INLINE unsigned int Viewer::add_mesh(const char* mesh_file_name,const std::string& id)
  {
    unsigned int mid = add_mesh(id);
    load_mesh_from_file(mesh_file_name);

#if defined(IGL_VIEWER_WITH_NANOGUI) && defined(IGL_VIEWER_WITH_NANOGUI_MULTIMESH)
    if(window != nullptr)
    {
      screen->performLayout();
    }
#endif

    return mid;
  }

  IGL_INLINE ViewerData& Viewer::get_mesh(unsigned int data_id)
  {
    assert(data_buffer.size() > data_id && "data_id out of range");

    if(data_id == active_data_id)
      return this->data;
    else
      return data_buffer[data_id];
  }

  IGL_INLINE bool Viewer::remove_mesh(unsigned int data_id)
  {
    assert(data_buffer.size() > data_id && "data_id out of range");
    assert(active_data_id != data_id && "active mesh cannot be removed!");
    
    if(active_data_id > data_id)
      active_data_id--;

    data_buffer.erase(data_buffer.begin()+data_id);
    data_ids.erase(data_ids.begin()+data_id);
    opengl[data_id].free();
    opengl.erase(opengl.begin()+data_id);

#if defined(IGL_VIEWER_WITH_NANOGUI) && defined(IGL_VIEWER_WITH_NANOGUI_MULTIMESH)
    if(window != nullptr)
    {
      currentDataCB->setItems(data_ids);
      screen->performLayout();
    }
#endif

    return true;
  }

  IGL_INLINE unsigned int Viewer::get_active_mesh()
  {
    return active_data_id;
  }

  IGL_INLINE bool Viewer::set_active_mesh(unsigned int data_id)
  {
    assert(data_buffer.size() > data_id && "data_id out of range");

    data_buffer[active_data_id] = ViewerData(data);
    active_data_id = data_id;
    data = data_buffer[active_data_id];

#if defined(IGL_VIEWER_WITH_NANOGUI) && defined(IGL_VIEWER_WITH_NANOGUI_MULTIMESH)
    if(window != nullptr)
    {
      currentDataCB->setSelectedIndex(active_data_id);
    }
#endif

    return true;
  }

  IGL_INLINE unsigned int Viewer::get_mesh_count()
  {
    return data_buffer.size();
  }

  IGL_INLINE bool Viewer::load_mesh_from_file(const char* mesh_file_name)
  {
    return load_mesh_from_file(mesh_file_name,active_data_id);
  }

  IGL_INLINE bool Viewer::load_mesh_from_file(const char* mesh_file_name,unsigned int data_id)
  {
    assert(data_buffer.size() > data_id && "data_id out of range");

    std::string mesh_file_name_string = std::string(mesh_file_name);

    // first try to load it with a plugin
    for(unsigned int i = 0; i<plugins.size(); ++i)
    {
      if(plugins[i]->load(mesh_file_name_string))
      {
        return true;
      }
    }

    ViewerData& data = get_mesh(data_id);
    data.clear();

    size_t last_dot = mesh_file_name_string.rfind('.');
    if(last_dot == std::string::npos)
    {
      printf("Error: No file extension found in %s\n",mesh_file_name);
      return false;
    }

    std::string extension = mesh_file_name_string.substr(last_dot+1);

    if(extension == "off" || extension =="OFF")
    {
      Eigen::MatrixXd V;
      Eigen::MatrixXi F;
      if(!igl::readOFF(mesh_file_name_string,V,F))
        return false;
      data.set_mesh(V,F);
    }
	else if(extension == "obj" || extension =="OBJ")
    {
      Eigen::MatrixXd corner_normals;
      Eigen::MatrixXi fNormIndices;

      Eigen::MatrixXd UV_V;
      Eigen::MatrixXi UV_F;
      Eigen::MatrixXd V;
      Eigen::MatrixXi F;

      if(!(
        igl::readOBJ(
          mesh_file_name_string,
          V,UV_V,corner_normals,F,UV_F,fNormIndices)))
      {
	    return false;
      }
      
      data.set_mesh(V,F);

      if(UV_V.rows() > 0)
      {
        data.set_uv(UV_V,UV_F);
      }

    }
	else
    {
      // unrecognized file type
      printf("Error: %s is not a recognized file type.\n",extension.c_str());
      return false;
    }

    data.compute_normals();
    data.uniform_colors(Eigen::Vector3d(51.0/255.0,43.0/255.0,33.3/255.0),
                        Eigen::Vector3d(255.0/255.0,228.0/255.0,58.0/255.0),
                        Eigen::Vector3d(255.0/255.0,235.0/255.0,80.0/255.0));
    if(data.V_uv.rows() == 0)
    {
      data.grid_texture();
    }

    core.align_camera_center(data);

    for(unsigned int i = 0; i<plugins.size(); ++i)
      if(plugins[i]->post_load())
        return true;

    return true;
  }

  IGL_INLINE bool Viewer::save_mesh_to_file(const char* mesh_file_name)
  {
    return save_mesh_to_file(mesh_file_name,active_data_id);
  }

  IGL_INLINE bool Viewer::save_mesh_to_file(const char* mesh_file_name, unsigned int data_id)
  {
    assert(data_buffer.size() > data_id && "data_id out of range");

    std::string mesh_file_name_string(mesh_file_name);

    // first try to load it with a plugin
    for(unsigned int i = 0; i<plugins.size(); ++i)
      if(plugins[i]->save(mesh_file_name_string))
        return true;

    ViewerData& data = get_mesh(data_id);

    size_t last_dot = mesh_file_name_string.rfind('.');
    if(last_dot == std::string::npos)
    {
      // No file type determined
      printf("Error: No file extension found in %s\n",mesh_file_name);
      return false;
    }
    std::string extension = mesh_file_name_string.substr(last_dot+1);
    if(extension == "off" || extension =="OFF")
    {
      return igl::writeOFF(mesh_file_name_string,data.V,data.F);
    } else if(extension == "obj" || extension =="OBJ")
    {
      Eigen::MatrixXd corner_normals;
      Eigen::MatrixXi fNormIndices;

      Eigen::MatrixXd UV_V;
      Eigen::MatrixXi UV_F;

      return igl::writeOBJ(mesh_file_name_string,data.V,
                           data.F,corner_normals,fNormIndices,UV_V,UV_F);
    } else
    {
      // unrecognized file type
      printf("Error: %s is not a recognized file type.\n",extension.c_str());
      return false;
    }
    return true;
  }

  IGL_INLINE bool Viewer::key_pressed(unsigned int unicode_key,int modifiers)
  {
    if (callback_key_pressed)
      if (callback_key_pressed(*this,unicode_key,modifiers))
        return true;

    for (unsigned int i = 0; i<plugins.size(); ++i)
    {
      if (plugins[i]->key_pressed(unicode_key, modifiers))
      {
        return true;
      }
    }

    if(key_pressed_handeld == false)
    {
      switch(unicode_key)
      {
        case 'A':
        case 'a':
        {
          core.is_animating = !core.is_animating;
          return true;
        }
		case 'F':
        case 'f':
        {
          data.set_face_based(!data.face_based);
          return true;
        }
        case 'I':
        case 'i':
        {
          data.dirty |= ViewerData::DIRTY_NORMAL;
          data.invert_normals = !data.invert_normals;
          return true;
        }
        case 'L':
        case 'l':
        {
          data.show_lines = !data.show_lines;
          return true;
        }
        case 'O':
        case 'o':
        {
          core.orthographic = !core.orthographic;
          return true;
        }
        case 'T':
        case 't':
        {
          data.show_faces = !data.show_faces;
          return true;
        }
        case 'Z':
        {
          snap_to_canonical_quaternion();
          return true;
        }
        case '[':
        case ']':
        {
          if(core.rotation_type == ViewerCore::ROTATION_TYPE_TRACKBALL)
          {
            core.set_rotation_type(
              ViewerCore::ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP);
          } else
          {
            core.set_rotation_type(ViewerCore::ROTATION_TYPE_TRACKBALL);
          }
          return true;
        }
#ifdef IGL_VIEWER_WITH_NANOGUI
        case ';':
          data.show_vertid = !data.show_vertid;
          return true;
        case ':':
          data.show_faceid = !data.show_faceid;
          return true;
#endif        
        default: break;//do nothing
      }
    }
    
    return false;
  }

  IGL_INLINE bool Viewer::key_down(int key,int modifiers)
  {
    key_pressed_handeld = true;

    if (callback_key_down)
      if (callback_key_down(*this,key,modifiers))
        return true;
    for (unsigned int i = 0; i<plugins.size(); ++i)
      if (plugins[i]->key_down(key, modifiers))
        return true;

    key_pressed_handeld = false;
    return false;
  }

  IGL_INLINE bool Viewer::key_up(int key,int modifiers)
  {
    if (callback_key_up)
      if (callback_key_up(*this,key,modifiers))
        return true;

    for (unsigned int i = 0; i<plugins.size(); ++i)
      if (plugins[i]->key_up(key, modifiers))
        return true;

    return false;
  }

  IGL_INLINE bool Viewer::mouse_down(MouseButton button,int modifier)
  {
    // Remember mouse location at down even if used by callback/plugin
    down_mouse_x = current_mouse_x;
    down_mouse_y = current_mouse_y;

    if (callback_mouse_down)
      if (callback_mouse_down(*this,static_cast<int>(button),modifier))
        return true;

    for (unsigned int i = 0; i<plugins.size(); ++i)
      if(plugins[i]->mouse_down(static_cast<int>(button),modifier))
        return true;

    // Initialization code for the trackball
    Eigen::RowVector3d center;
    if (data.V.rows() == 0)
    {
      center << 0,0,0;
    }else
    {
      center = data.V.colwise().sum()/data.V.rows();
    }

    Eigen::Vector3f coord =
      igl::project(
        Eigen::Vector3f(center(0),center(1),center(2)),
        (core.view * data.model).eval(),
        core.proj,
        core.viewport);
    down_mouse_z = coord[2];

    down = true;

    down_modifier = modifier;
    if(down_modifier == GLFW_MOD_CONTROL)
    {
      down_translation = data.model_translation;
    } else
    {
      down_translation = core.global_translation;
    }

    down_rotation = core.trackball_angle;

    mouse_mode = MouseMode::Rotation;

    switch (button)
    {
      case MouseButton::Left:
        if(down_modifier == GLFW_MOD_SHIFT)
        {
          mouse_mode = MouseMode::Translation;
        }
        else
        {
          mouse_mode = MouseMode::Rotation;
        }
        break;

      case MouseButton::Right:
        mouse_mode = MouseMode::Translation;
        break;

      default:
        mouse_mode = MouseMode::None;
        break;
    }

    return true;
  }

  IGL_INLINE bool Viewer::mouse_up(MouseButton button,int modifier)
  {
    down = false;

    if (callback_mouse_up)
      if (callback_mouse_up(*this,static_cast<int>(button),modifier))
        return true;

    for (unsigned int i = 0; i<plugins.size(); ++i)
      if(plugins[i]->mouse_up(static_cast<int>(button),modifier))
          return true;

    mouse_mode = MouseMode::None;

    return true;
  }

  IGL_INLINE bool Viewer::mouse_move(int mouse_x,int mouse_y)
  {
    if(hack_never_moved)
    {
      down_mouse_x = mouse_x;
      down_mouse_y = mouse_y;
      hack_never_moved = false;
    }
    current_mouse_x = mouse_x;
    current_mouse_y = mouse_y;

    if (callback_mouse_move)
      if (callback_mouse_move(*this,mouse_x,mouse_y))
        return true;

    for (unsigned int i = 0; i<plugins.size(); ++i)
      if (plugins[i]->mouse_move(mouse_x, mouse_y))
        return true;

    if (down)
    {
      switch (mouse_mode)
      {
        case MouseMode::Rotation:
        {
          switch(core.rotation_type)
          {
            default:
              assert(false && "Unknown rotation type");
            case ViewerCore::ROTATION_TYPE_TRACKBALL:
            {
              igl::trackball(
                core.viewport(2),
                core.viewport(3),
                2.0f,
                down_rotation,
                down_mouse_x,
                down_mouse_y,
                mouse_x,
                mouse_y,
                core.trackball_angle);
            }
              break;
            case ViewerCore::ROTATION_TYPE_TWO_AXIS_VALUATOR_FIXED_UP:
              igl::two_axis_valuator_fixed_up(
                core.viewport(2),core.viewport(3),
                2.0,
                down_rotation,
                down_mouse_x, down_mouse_y, mouse_x, mouse_y,
                core.trackball_angle);
              break;
          }
          //Eigen::Vector4f snapq = core.trackball_angle;

          break;
        }

        case MouseMode::Translation:
        {
          //translation
          Eigen::Vector3f pos1 = igl::unproject(Eigen::Vector3f(mouse_x, core.viewport[3] - mouse_y, down_mouse_z), (core.view * data.model).eval(), core.proj, core.viewport);
          Eigen::Vector3f pos0 = igl::unproject(Eigen::Vector3f(down_mouse_x, core.viewport[3] - down_mouse_y, down_mouse_z), (core.view * data.model).eval(), core.proj, core.viewport);

          Eigen::Vector3f diff = pos1 - pos0;
          
          if(down_modifier == GLFW_MOD_CONTROL)
          {
            data.model_translation = down_translation + Eigen::Vector3f(diff[0],diff[1],diff[2]);
          }
          else
          {
            core.global_translation = down_translation + Eigen::Vector3f(diff[0],diff[1],diff[2]);
          }

          break;
        }
        case MouseMode::Zoom:
        {
          float delta = 0.001f * (mouse_x - down_mouse_x + mouse_y - down_mouse_y);
          core.camera_zoom *= 1 + delta;
          down_mouse_x = mouse_x;
          down_mouse_y = mouse_y;
          break;
        }

        default:
          break;
      }
    }
    return true;
  }

  IGL_INLINE bool Viewer::mouse_scroll(float delta_y)
  {
    scroll_position += delta_y;

    if (callback_mouse_scroll)
      if (callback_mouse_scroll(*this,delta_y))
        return true;

    for (unsigned int i = 0; i<plugins.size(); ++i)
      if (plugins[i]->mouse_scroll(delta_y))
        return true;

    // Only zoom if there's actually a change
    if(delta_y != 0)
    {
      float mult = (1.0+((delta_y>0)?1.:-1.)*0.05);
      const float min_zoom = 0.0001f;
      core.camera_zoom = (core.camera_zoom * mult > min_zoom ? core.camera_zoom * mult : min_zoom);
    }
    return true;
  }

  IGL_INLINE void Viewer::draw()
  {
    using namespace std;
    using namespace Eigen;

    core.clear_framebuffers();

    if (callback_pre_draw)
      if (callback_pre_draw(*this))
        return;

    for (unsigned int i = 0; i<plugins.size(); ++i)
      if (plugins[i]->pre_draw())
        return;

    for(int i=0;i<data_buffer.size();i++)
    {
      if(i != active_data_id)
        core.draw(data_buffer[i],opengl[i]);
      else
        core.draw(data,opengl[active_data_id]);
    }

    if (callback_post_draw)
      if (callback_post_draw(*this))
        return;

    for (unsigned int i = 0; i<plugins.size(); ++i)
      if (plugins[i]->post_draw())
        break;

#ifdef IGL_VIEWER_WITH_NANOGUI
  ngui->refresh();
	screen->drawContents();
	screen->drawWidgets();
#endif
  }

  IGL_INLINE bool Viewer::save_scene()
  {
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
      return false;

#ifdef IGL_VIEWER_WITH_NANOGUI_SERIALIZATION

    igl::serialize(window_maximized,"window_maximized",fname.c_str(),true);
    igl::serialize(window_position,"window_position",fname.c_str());
    igl::serialize(window_size,"window_size",fname.c_str());
    igl::serialize(core,"Core",fname.c_str());

#ifndef ENABLE_SERIALIZATION_CORE_ONLY
    data_buffer[active_data_id] = data;

    igl::serialize(active_data_id,"ActiveDataId",fname.c_str());
    igl::serialize(data_buffer,"Data",fname.c_str());
    for(unsigned int i = 0; i <plugins.size(); ++i)
      igl::serialize(*plugins[i],plugins[i]->plugin_name,fname.c_str());
#endif

#endif

    return true;
  }

  IGL_INLINE bool Viewer::load_scene()
  {
    std::string fname = igl::file_dialog_open();
    if(fname.length() == 0)
      return false;

    return load_scene(fname);
  }

  IGL_INLINE bool Viewer::load_scene(std::string fname)
  {
#ifdef IGL_VIEWER_WITH_NANOGUI_SERIALIZATION

    igl::deserialize(window_maximized,"window_maximized",fname.c_str());
    igl::deserialize(window_position,"window_position",fname.c_str());
    igl::deserialize(window_size,"window_size",fname.c_str());
    igl::deserialize(core,"Core",fname.c_str());

#ifndef ENABLE_SERIALIZATION_CORE_ONLY
    igl::deserialize(active_data_id,"ActiveDataId",fname.c_str());
    igl::deserialize(data_buffer,"Data",fname.c_str());
    for(unsigned int i = 0; i <plugins.size(); ++i)
      igl::deserialize(*plugins[i],plugins[i]->plugin_name,fname.c_str());
    
    data = data_buffer[active_data_id];
#endif

    if(window_position(0) < 0 || window_position(1) < 0 || window_size(0) < 0 || window_size(1) < 0)
    {
      window_position << 200,200;
      window_size << 640,480;
    }

    if(window_maximized)
    {
      glfwSetWindowPos(window,window_position(0),window_position(1));
      glfwSetWindowSize(window,window_size(0),window_size(1));
#ifdef GLFW_VERSION_HIGHER_THAN_3_1
      glfwMaximizeWindow(window);
#endif
    }
    else
    {
      glfwSetWindowPos(window,window_position(0),window_position(1));
      glfwSetWindowSize(window,window_size(0),window_size(1));
    }

#endif

    return true;
  }

  IGL_INLINE void Viewer::resize(int w,int h)
  {
#ifdef GLFW_VERSION_HIGHER_THAN_3_1
    window_maximized = glfwGetWindowAttrib(window,GLFW_MAXIMIZED);
#endif

    if(!window_maximized)
    {
      window_size << w,h;
    }

    core.viewport = Eigen::Vector4f(0,0,w,h);
  }

  IGL_INLINE void Viewer::snap_to_canonical_quaternion()
  {
    Eigen::Quaternionf snapq = this->core.trackball_angle;
    igl::snap_to_canonical_view_quat(snapq,1.0f,this->core.trackball_angle);
  }

  IGL_INLINE void Viewer::open_dialog_load_mesh()
  {
    std::string fname = igl::file_dialog_open();

    if (fname.length() == 0)
      return;

    this->load_mesh_from_file(fname.c_str());
  }

  IGL_INLINE void Viewer::open_dialog_save_mesh()
  {
    std::string fname = igl::file_dialog_save();

    if(fname.length() == 0)
      return;

    this->save_mesh_to_file(fname.c_str());
  }


  IGL_INLINE int  Viewer::launch_init(bool resizable,bool fullscreen,int width, int height)
  {
    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
      return EXIT_FAILURE;

    glfwWindowHint(GLFW_SAMPLES, 8);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 3);

    #ifdef __APPLE__
      glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
      glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    #endif

    if(fullscreen)
    {
      GLFWmonitor *monitor = glfwGetPrimaryMonitor();
      const GLFWvidmode *mode = glfwGetVideoMode(monitor);
      window = glfwCreateWindow(mode->width,mode->height,"libigl viewer",monitor,nullptr);
    }
    else
    {
      window = glfwCreateWindow(width,height,"libigl viewer",nullptr,nullptr);
    }

    if (!window)
    {
      glfwTerminate();
      return EXIT_FAILURE;
    }

    glfwMakeContextCurrent(window);

    #ifndef __APPLE__
      glewExperimental = true;
      GLenum err = glewInit();
      if(GLEW_OK != err)
      {
        /* Problem: glewInit failed, something is seriously wrong. */
       fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
      }
      glGetError(); // pull and savely ignonre unhandled errors like GL_INVALID_ENUM
      fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
    #endif

    #if defined(DEBUG) || defined(_DEBUG)
      int major, minor, rev;
      major = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MAJOR);
      minor = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MINOR);
      rev = glfwGetWindowAttrib(window, GLFW_CONTEXT_REVISION);
      printf("OpenGL version recieved: %d.%d.%d\n", major, minor, rev);
      printf("Supported OpenGL is %s\n", (const char*)glGetString(GL_VERSION));
      printf("Supported GLSL is %s\n", (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION));
    #endif

    glfwSetInputMode(window,GLFW_CURSOR,GLFW_CURSOR_NORMAL);

    // Initialize FormScreen
#ifdef IGL_VIEWER_WITH_NANOGUI
    screen = new nanogui::Screen();
    screen->initialize(window, false);
    ngui = new nanogui::FormHelper(screen);
#endif

    __viewer = this;

    // Register callbacks
    glfwSetKeyCallback(window, glfw_key_callback);
    glfwSetCursorPosCallback(window,glfw_mouse_move);
    glfwSetWindowPosCallback(window,glfw_window_position);
    glfwSetWindowSizeCallback(window,glfw_window_size);
    glfwSetMouseButtonCallback(window,glfw_mouse_press);
    glfwSetScrollCallback(window,glfw_mouse_scroll);
    glfwSetCharModsCallback(window,glfw_char_mods_callback);
    glfwSetDropCallback(window,glfw_drop_callback);

    // Handle retina displays (windows and mac)
    glfwGetFramebufferSize(window, &width, &height);
    
    glfwGetWindowPos(window,&window_position(0),&window_position(1));
    glfwGetWindowSize(window, &window_size(0), &window_size(1));

    highdpi = width/window_size(0);

    glfw_window_size(window,window_size(0),window_size(1));

    for(auto& v : opengl)
    {
      v.init();
    }

    core.align_camera_center(data);

    // Initialize IGL viewer
    init();
    return EXIT_SUCCESS;
  }

  IGL_INLINE bool Viewer::launch_rendering(bool loop)
  {
    // glfwMakeContextCurrent(window);

    // Rendering loop
    while (!glfwWindowShouldClose(window))
    {
      double tic = get_seconds();
      draw();

      glfwSwapBuffers(window);
      if(core.is_animating)
      {
        glfwPollEvents();
        // In microseconds
        double duration = 1000000.*(get_seconds()-tic);
        const double min_duration = 1000000./core.animation_max_fps;
        if(duration<min_duration)
        {
          std::this_thread::sleep_for(std::chrono::microseconds((int)(min_duration-duration)));
        }
      }
      else
      {
        glfwWaitEvents();
      }

      if (!loop)
        return !glfwWindowShouldClose(window);
    }
    return EXIT_SUCCESS;
  }

  IGL_INLINE void Viewer::launch_shut()
  {
    for(auto&& v : opengl)
      v.free();
    core.shut();

    shutdown_plugins();
#ifdef IGL_VIEWER_WITH_NANOGUI
    delete ngui;
    //delete screen;
    screen = nullptr;
    ngui = nullptr;
#endif

    glfwDestroyWindow(window);
    glfwTerminate();
    return;
  }

  IGL_INLINE int Viewer::launch(bool resizable,bool fullscreen,int width,int height)
  {
    // TODO return values are being ignored...
    launch_init(resizable,fullscreen,width,height);
    launch_rendering(true);
    launch_shut();
    return EXIT_SUCCESS;
  }
} // end namespace
}