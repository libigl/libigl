// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

// TODO: save_scene()/load_scene()

#include "Viewer.h"
#include <igl/get_seconds.h>

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

#define GLFW_INCLUDE_GLU
#include <GLFW/glfw3.h>

#include <cmath>
#include <cstdio>
#include <sstream>
#include <iomanip>
#include <iostream>
#include <fstream>

#include <algorithm>
#include <igl/project.h>

#include <limits>
#include <cassert>

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
#include <igl/snap_to_canonical_view_quat.h>
#include <igl/unproject.h>
#include <igl/viewer/TextRenderer.h>

#ifdef ENABLE_SERIALIZATION
#include <igl/serialize.h>
#endif

// Internal global variables used for glfw event handling
static igl::Viewer * __viewer;
static double highdpi = 1;
static double scroll_x = 0;
static double scroll_y = 0;
static int global_KMod = 0;

namespace {
void TW_CALL copy_str(std::string& dst, const std::string& src)
{
  dst = src;
}
}

static void glfw_mouse_press(GLFWwindow* window, int button, int action, int modifier)
{
  bool tw_used = TwEventMouseButtonGLFW(button, action);
  igl::Viewer::MouseButton mb;

  if (button == GLFW_MOUSE_BUTTON_1)
    mb = igl::Viewer::IGL_LEFT;
  else if (button == GLFW_MOUSE_BUTTON_2)
    mb = igl::Viewer::IGL_RIGHT;
  else //if (button == GLFW_MOUSE_BUTTON_3)
    mb = igl::Viewer::IGL_MIDDLE;

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

IGL_INLINE int TwEventKeyGLFW3(int glfwKey, int glfwAction)
{
  int handled = 0;

  // Register of modifiers state
  if (glfwAction==GLFW_PRESS)
  {
    switch (glfwKey)
    {
      case GLFW_KEY_LEFT_SHIFT:
      case GLFW_KEY_RIGHT_SHIFT:
        global_KMod |= TW_KMOD_SHIFT;
        break;
      case GLFW_KEY_LEFT_CONTROL:
      case GLFW_KEY_RIGHT_CONTROL:
        global_KMod |= TW_KMOD_CTRL;
        break;
      case GLFW_KEY_LEFT_ALT:
      case GLFW_KEY_RIGHT_ALT:
        global_KMod |= TW_KMOD_ALT;
        break;
    }
  }
  else
  {
    switch (glfwKey)
    {
      case GLFW_KEY_LEFT_SHIFT:
      case GLFW_KEY_RIGHT_SHIFT:
        global_KMod &= ~TW_KMOD_SHIFT;
        break;
      case GLFW_KEY_LEFT_CONTROL:
      case GLFW_KEY_RIGHT_CONTROL:
        global_KMod &= ~TW_KMOD_CTRL;
        break;
      case GLFW_KEY_LEFT_ALT:
      case GLFW_KEY_RIGHT_ALT:
        global_KMod &= ~TW_KMOD_ALT;
        break;
    }
  }

  // Process key pressed
  if (glfwAction==GLFW_PRESS)
  {
    int mod = global_KMod;
    int testkp = ((mod&TW_KMOD_CTRL) || (mod&TW_KMOD_ALT)) ? 1 : 0;

    if ((mod&TW_KMOD_CTRL) && glfwKey>0 && glfwKey<GLFW_KEY_ESCAPE   )   // CTRL cases
      handled = TwKeyPressed(glfwKey, mod);
    else if (glfwKey>=GLFW_KEY_ESCAPE  )
    {
      int k = 0;

      if (glfwKey>=GLFW_KEY_F1 && glfwKey<=GLFW_KEY_F15)
        k = TW_KEY_F1 + (glfwKey-GLFW_KEY_F1);
      else if (testkp && glfwKey>=GLFW_KEY_KP_0 && glfwKey<=GLFW_KEY_KP_9)
        k = '0' + (glfwKey-GLFW_KEY_KP_0);
      else
      {
        switch (glfwKey)
        {
          case GLFW_KEY_ESCAPE  :
            k = TW_KEY_ESCAPE;
            break;
          case GLFW_KEY_UP:
            k = TW_KEY_UP;
            break;
          case GLFW_KEY_DOWN:
            k = TW_KEY_DOWN;
            break;
          case GLFW_KEY_LEFT:
            k = TW_KEY_LEFT;
            break;
          case GLFW_KEY_RIGHT:
            k = TW_KEY_RIGHT;
            break;
          case GLFW_KEY_TAB:
            k = TW_KEY_TAB;
            break;
          case GLFW_KEY_ENTER:
            k = TW_KEY_RETURN;
            break;
          case GLFW_KEY_BACKSPACE:
            k = TW_KEY_BACKSPACE;
            break;
          case GLFW_KEY_INSERT:
            k = TW_KEY_INSERT;
            break;
          case GLFW_KEY_DELETE:
            k = TW_KEY_DELETE;
            break;
          case GLFW_KEY_PAGE_UP:
            k = TW_KEY_PAGE_UP;
            break;
          case GLFW_KEY_PAGE_DOWN:
            k = TW_KEY_PAGE_DOWN;
            break;
          case GLFW_KEY_HOME:
            k = TW_KEY_HOME;
            break;
          case GLFW_KEY_END:
            k = TW_KEY_END;
            break;
          case GLFW_KEY_KP_ENTER:
            k = TW_KEY_RETURN;
            break;
          case GLFW_KEY_KP_DIVIDE:
            if (testkp)
              k = '/';
            break;
          case GLFW_KEY_KP_MULTIPLY:
            if (testkp)
              k = '*';
            break;
          case GLFW_KEY_KP_SUBTRACT:
            if (testkp)
              k = '-';
            break;
          case GLFW_KEY_KP_ADD:
            if (testkp)
              k = '+';
            break;
          case GLFW_KEY_KP_DECIMAL:
            if (testkp)
              k = '.';
            break;
          case GLFW_KEY_KP_EQUAL:
            if (testkp)
              k = '=';
            break;
        }
      }

      if (k>0)
        handled = TwKeyPressed(k, mod);
    }
  }

  return handled;
}

static void glfw_key_callback(GLFWwindow* window, int key, int scancode, int action, int modifier)
{
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    glfwSetWindowShouldClose(window, GL_TRUE);

  if (!TwEventKeyGLFW3(key,action))
  {
    if (action == GLFW_PRESS)
      __viewer->key_down(key, modifier);
    else
      __viewer->key_up(key, modifier);
  }
}

static void glfw_window_size(GLFWwindow* window, int width, int height)
{
  int w = width*highdpi;
  int h = height*highdpi;

  __viewer->resize(w, h);

  TwWindowSize(w, h);
  const auto & bar = __viewer->bar;
  // Keep AntTweakBar on right side of screen and height == opengl height
  // get the current position of a bar
  int size[2];
  TwGetParam(bar, NULL, "size", TW_PARAM_INT32, 2, size);
  int pos[2];
  // Place bar on left side of opengl rect (padded by 10 pixels)
  pos[0] = 10;//max(10,(int)width - size[0] - 10);
  // place bar at top (padded by 10 pixels)
  pos[1] = 10;
  // Set height to new height of window (padded by 10 pixels on bottom)
  size[1] = highdpi*(height-pos[1]-10);
  TwSetParam(bar, NULL, "position", TW_PARAM_INT32, 2, pos);
  TwSetParam(bar, NULL, "size", TW_PARAM_INT32, 2,size);
}

static void glfw_mouse_move(GLFWwindow* window, double x, double y)
{
  if(!TwEventMousePosGLFW(x*highdpi,y*highdpi) || __viewer->down)
  {
    // Call if TwBar hasn't used or if down
    __viewer->mouse_move(x*highdpi, y*highdpi);
  }
}

static void glfw_mouse_scroll(GLFWwindow* window, double x, double y)
{
  using namespace std;
  scroll_x += x;
  scroll_y += y;

  if (!TwEventMouseWheelGLFW(scroll_y))
    __viewer->mouse_scroll(y);
}

static void glfw_char_callback(GLFWwindow* window, unsigned int c)
{
  if ((c & 0xff00)==0)
    TwKeyPressed(c, global_KMod);
}

namespace igl
{
  IGL_INLINE void Viewer::init()
  {
    // Create a tweak bar
    bar = TwNewBar("libIGL-Viewer");
    TwDefine(" libIGL-Viewer help='This is a simple 3D mesh viewer.' "); // Message added to the help bar->
    TwDefine(" libIGL-Viewer size='200 685'"); // change default tweak bar size
    TwDefine(" libIGL-Viewer color='76 76 127' "); // change default tweak bar color
    TwDefine(" libIGL-Viewer refresh=0.5"); // change refresh rate

    // ---------------------- LOADING ----------------------

    #ifdef ENABLE_SERIALIZATION
    TwAddButton(bar,"Load Scene", load_scene_cb,    this, "group='Workspace'");
    TwAddButton(bar,"Save Scene", save_scene_cb,    this, "group='Workspace'");
    #endif

    #ifdef ENABLE_IO
    TwAddButton(bar,"Load Mesh",  open_dialog_mesh, this, "group='Mesh' key=o");
    #endif

    // ---------------------- SCENE ----------------------

    TwAddButton(bar,"Center object", align_camera_center_cb, this,
                " group='Viewing Options'"
                " label='Center object' key=A help='Set the center of the camera to the mesh center.'");
    TwAddVarRW(bar, "Zoom", TW_TYPE_FLOAT, &(core.camera_zoom),
               " min=0.05 max=50 step=0.1 keyIncr=+ keyDecr=- help='Scale the object (1=original size).' group='Scene'");
    TwAddButton(bar,"SnapView", snap_to_canonical_quaternion_cb, this,
                " group='Scene'"
                " label='Snap to canonical view' key=Z "
                " help='Snaps view to nearest canonical view.'");
    TwAddVarRW(bar,"LightDir", TW_TYPE_DIR3F, core.light_position.data(),
               " group='Scene'"
               " label='Light direction' open help='Change the light direction.' ");

    // ---------------------- DRAW OPTIONS ----------------------
    TwAddVarRW(bar, "ToggleOrthographic", TW_TYPE_BOOLCPP, &(core.orthographic),
               " group='Viewing Options'"
               " label='Orthographic view' "
               " help='Toggles orthographic / perspective view. Default: perspective.'");
    TwAddVarRW(bar, "Rotation", TW_TYPE_QUAT4F, &(core.trackball_angle),
      " group='Viewing Options'"
      " label='Rotation'"
      " help='Rotates view.'");
    TwAddVarCB(bar,"Face-based Normals/Colors", TW_TYPE_BOOLCPP, set_face_based_cb, get_face_based_cb, this,
               " group='Draw options'"
               " label='Face-based' key=T help='Toggle per face shading/colors.' ");

    TwAddVarRW(bar,"Show texture", TW_TYPE_BOOLCPP, &(core.show_texture),
               " group='Draw options'");

    TwAddVarCB(bar,"Invert Normals", TW_TYPE_BOOLCPP, set_invert_normals_cb, get_invert_normals_cb, this,
               " group='Draw options'"
               " label='Invert normals' key=i help='Invert normal directions for inside out meshes.' ");

    TwAddVarRW(bar,"ShowOverlay", TW_TYPE_BOOLCPP, &(core.show_overlay),
               " group='Draw options'"
               " label='Show overlay' key=o help='Show the overlay layers.' ");
    TwAddVarRW(bar,"ShowOverlayDepth", TW_TYPE_BOOLCPP, &(core.show_overlay_depth),
               " group='Draw options'"
               " label='Show overlay depth test' help='Enable the depth test for overlay layer.' ");
    TwAddVarRW(bar,"Background color", TW_TYPE_COLOR3F,
               core.background_color.data(),
               " help='Select a background color' colormode=hls group='Draw options'");
    TwAddVarRW(bar, "LineColor", TW_TYPE_COLOR3F,
               core.line_color.data(),
               " label='Line color' help='Select a outline color' group='Draw options'");
    TwAddVarRW(bar,"Shininess",TW_TYPE_FLOAT,&core.shininess," group='Draw options'"
               " min=1 max=128");

    // ---------------------- Overlays ----------------------

    TwAddVarRW(bar,"Wireframe", TW_TYPE_BOOLCPP, &(core.show_lines),
               " group='Overlays'"
               " label='Wireframe' key=l help='Toggle wire frame of mesh'");
    TwAddVarRW(bar,"Fill", TW_TYPE_BOOLCPP, &(core.show_faces),
               " group='Overlays'"
               " label='Fill' key=t help='Display filled polygons of mesh'");
    TwAddVarRW(bar,"ShowVertexId", TW_TYPE_BOOLCPP, &(core.show_vertid),
               " group='Overlays'"
               " label='Show Vertex Labels' key=';' help='Toggle vertex indices'");
    TwAddVarRW(bar,"ShowFaceId", TW_TYPE_BOOLCPP, &(core.show_faceid),
               " group='Overlays'"
               " label='Show Faces Labels' key='CTRL+;' help='Toggle face"
               " indices'");

    core.init();

    if (callback_init)
      if (callback_init(*this))
        return;

    init_plugins();
  }

  IGL_INLINE Viewer::Viewer()
  {

    // Temporary variables initialization
    down = false;
    hack_never_moved = true;
    scroll_position = 0.0f;

    // Per face
    data.set_face_based(false);

    // C-style callbacks
    callback_init         = 0;
    callback_pre_draw     = 0;
    callback_post_draw    = 0;
    callback_mouse_down   = 0;
    callback_mouse_up     = 0;
    callback_mouse_move   = 0;
    callback_mouse_scroll = 0;
    callback_key_down     = 0;
    callback_key_up       = 0;

    callback_init_data          = 0;
    callback_pre_draw_data      = 0;
    callback_post_draw_data     = 0;
    callback_mouse_down_data    = 0;
    callback_mouse_up_data      = 0;
    callback_mouse_move_data    = 0;
    callback_mouse_scroll_data  = 0;
    callback_key_down_data      = 0;
    callback_key_up_data        = 0;

  }

  IGL_INLINE void Viewer::init_plugins()
  {
    // Init all plugins
    for (unsigned int i = 0; i<plugins.size(); ++i)
      plugins[i]->init(this);
  }

  IGL_INLINE Viewer::~Viewer()
  {
  }

  IGL_INLINE void Viewer::shutdown_plugins()
  {
    for (unsigned int i = 0; i<plugins.size(); ++i)
      plugins[i]->shutdown();
  }

  IGL_INLINE bool Viewer::load_mesh_from_file(const char* mesh_file_name)
  {
    std::string mesh_file_name_string = std::string(mesh_file_name);

    // first try to load it with a plugin
    for (unsigned int i = 0; i<plugins.size(); ++i)
    {
      if (plugins[i]->load(mesh_file_name_string))
      {
        return true;
      }
    }

    data.clear();

    size_t last_dot = mesh_file_name_string.rfind('.');
    if (last_dot == std::string::npos)
    {
      printf("Error: No file extension found in %s\n",mesh_file_name);
      return false;
    }

    std::string extension = mesh_file_name_string.substr(last_dot+1);

    if (extension == "off" || extension =="OFF")
    {
      Eigen::MatrixXd V;
      Eigen::MatrixXi F;
      if (!igl::readOFF(mesh_file_name_string, V, F))
        return false;
      data.set_mesh(V,F);
    }
    else if (extension == "obj" || extension =="OBJ")
    {
      Eigen::MatrixXd corner_normals;
      Eigen::MatrixXi fNormIndices;

      Eigen::MatrixXd UV_V;
      Eigen::MatrixXi UV_F;
      Eigen::MatrixXd V;
      Eigen::MatrixXi F;

      if (!(igl::readOBJ(mesh_file_name_string, V, F, corner_normals, fNormIndices, UV_V, UV_F)))
        return false;

      data.set_mesh(V,F);
      data.set_uv(UV_V,UV_F);

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
    if (data.V_uv.rows() == 0)
      data.grid_texture();

    core.align_camera_center(data.V,data.F);

    for (unsigned int i = 0; i<plugins.size(); ++i)
      if (plugins[i]->post_load())
        return true;

    return true;
  }

  IGL_INLINE bool Viewer::save_mesh_to_file(const char* mesh_file_name)
  {
    std::string mesh_file_name_string(mesh_file_name);

    // first try to load it with a plugin
    for (unsigned int i = 0; i<plugins.size(); ++i)
      if (plugins[i]->save(mesh_file_name_string))
        return true;

    size_t last_dot = mesh_file_name_string.rfind('.');
    if (last_dot == std::string::npos)
    {
      // No file type determined
      printf("Error: No file extension found in %s\n",mesh_file_name);
      return false;
    }
    std::string extension = mesh_file_name_string.substr(last_dot+1);
    if (extension == "off" || extension =="OFF")
    {
      return igl::writeOFF(mesh_file_name_string,data.V,data.F);
    }
    else if (extension == "obj" || extension =="OBJ")
    {
      Eigen::MatrixXd corner_normals;
      Eigen::MatrixXi fNormIndices;

      Eigen::MatrixXd UV_V;
      Eigen::MatrixXi UV_F;

      return igl::writeOBJ(mesh_file_name_string, data.V,
          data.F, corner_normals, fNormIndices, UV_V, UV_F);
    }
    else
    {
      // unrecognized file type
      printf("Error: %s is not a recognized file type.\n",extension.c_str());
      return false;
    }
    return true;
  }

  IGL_INLINE bool Viewer::key_down(int key,int modifiers)
  {
    if (callback_key_down)
      if (callback_key_down(*this,key,modifiers))
        return true;

    for (unsigned int i = 0; i<plugins.size(); ++i)
      if (plugins[i]->key_down(key, modifiers))
        return true;

    char k = key;

    if(key == GLFW_KEY_S && modifiers == GLFW_MOD_SHIFT)
      mouse_scroll(1);

    if(key == GLFW_KEY_A && modifiers == GLFW_MOD_SHIFT)
      mouse_scroll(-1);

    // Why aren't these handled via AntTweakBar?
    if(key == GLFW_KEY_Z) // Don't use 'Z' because that clobbers snap_to_canonical_view_quat
      core.trackball_angle << 0.0f, 0.0f, 0.0f, 1.0f;

    if(key == GLFW_KEY_Y)
      core.trackball_angle << -sqrt(2.0f)/2.0f, 0.0f, 0.0f, sqrt(2.0f)/2.0f;

    if(key == GLFW_KEY_X)
      core.trackball_angle << -0.5f, -0.5f, -0.5f, 0.5f;


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
    if (callback_mouse_down)
      if (callback_mouse_down(*this,button,modifier))
        return true;

    for (unsigned int i = 0; i<plugins.size(); ++i)
      if (plugins[i]->mouse_down(button,modifier))
        return true;

    down = true;

    down_mouse_x = current_mouse_x;
    down_mouse_y = current_mouse_y;
    down_translation = core.model_translation;


    // Initialization code for the trackball
    Eigen::RowVector3d center;
    if (data.V.rows() == 0)
      center << 0,0,0;
    else
      center = data.V.colwise().sum()/data.V.rows();

    Eigen::Vector3f coord = igl::project(Eigen::Vector3f(center(0),center(1),center(2)), core.view * core.model, core.proj, core.viewport);
    down_mouse_z = coord[2];
    down_rotation = core.trackball_angle;

    mouse_mode = ROTATION;

    switch (button)
    {
      case IGL_LEFT:
        mouse_mode = ROTATION;
        break;

      case IGL_RIGHT:
        mouse_mode = TRANSLATE;
        break;

      default:
        mouse_mode = NOTHING;
        break;
    }

    return true;
  }

  IGL_INLINE bool Viewer::mouse_up(MouseButton button,int modifier)
  {
    down = false;

    if (callback_mouse_up)
      if (callback_mouse_up(*this,button,modifier))
        return true;

    for (unsigned int i = 0; i<plugins.size(); ++i)
        if (plugins[i]->mouse_up(button,modifier))
          return true;

    mouse_mode = NOTHING;

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
        case ROTATION :
        {
          igl::trackball(core.viewport(2),
                         core.viewport(3),
                         2.0f,
                         down_rotation.data(),
                         down_mouse_x,
                         down_mouse_y,
                         mouse_x,
                         mouse_y,
                         core.trackball_angle.data());
          //Eigen::Vector4f snapq = core.trackball_angle;

          break;
        }

        case TRANSLATE:
        {
          //translation
          Eigen::Vector3f pos1 = igl::unproject(Eigen::Vector3f(mouse_x, core.viewport[3] - mouse_y, down_mouse_z), core.view * core.model, core.proj, core.viewport);
          Eigen::Vector3f pos0 = igl::unproject(Eigen::Vector3f(down_mouse_x, core.viewport[3] - down_mouse_y, down_mouse_z), core.view * core.model, core.proj, core.viewport);

          Eigen::Vector3f diff = pos1 - pos0;
          core.model_translation = down_translation + Eigen::Vector3f(diff[0],diff[1],diff[2]);

          break;
        }
        case ZOOM:
        {
          //float delta = 0.001f * (mouse_x - down_mouse_x + mouse_y - down_mouse_y);
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
      const float min_zoom = 0.1f;
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

    core.draw(data,opengl);

    if (callback_post_draw)
      if (callback_post_draw(*this))
        return;

    for (unsigned int i = 0; i<plugins.size(); ++i)
      if (plugins[i]->post_draw())
        break;

    TwDraw();
  }

  IGL_INLINE bool Viewer::save_scene()
  {
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
      return false;

#ifdef ENABLE_SERIALIZATION

    igl::serialize(core,"Core",fname.c_str(),true);
    igl::serialize(data,"Data",fname.c_str());

    for(unsigned int i = 0; i <plugins.size(); ++i)
      igl::serialize(*plugins[i],plugins[i]->plugin_name,fname.c_str());

#endif

    return true;
  }

  IGL_INLINE bool Viewer::load_scene()
  {
    std::string fname = igl::file_dialog_open();
    if (fname.length() == 0)
      return false;

#ifdef ENABLE_SERIALIZATION

    igl::deserialize(core,"Core",fname.c_str());
    igl::deserialize(data,"Data",fname.c_str());

    for(unsigned int i = 0; i <plugins.size(); ++i)
      igl::deserialize(*plugins[i],plugins[i]->plugin_name,fname.c_str());

#endif

    return true;
  }

  IGL_INLINE void Viewer::resize(int w,int h)
  {
    core.viewport = Eigen::Vector4f(0,0,w,h);
  }

  IGL_INLINE void TW_CALL Viewer::snap_to_canonical_quaternion_cb(void *clientData)
  {
    Eigen::Vector4f snapq = static_cast<Viewer *>(clientData)->core.trackball_angle;
    igl::snap_to_canonical_view_quat<float>(snapq.data(),1,static_cast<Viewer *>(clientData)->core.trackball_angle.data());
  }
  IGL_INLINE void TW_CALL Viewer::align_camera_center_cb(void *clientData)
  {
    static_cast<Viewer *>(clientData)->core.align_camera_center(
      static_cast<Viewer *>(clientData)->data.V,
      static_cast<Viewer *>(clientData)->data.F);
  }

  IGL_INLINE void TW_CALL Viewer::save_scene_cb(void *clientData)
  {
    static_cast<Viewer *>(clientData)->save_scene();
  }

  IGL_INLINE void TW_CALL Viewer::load_scene_cb(void *clientData)
  {
    static_cast<Viewer *>(clientData)->load_scene();
  }

  IGL_INLINE void TW_CALL Viewer::set_invert_normals_cb(const void *param,void *clientData)
  {
    Viewer *viewer = static_cast<Viewer *>(clientData);
    viewer->data.dirty |= ViewerData::DIRTY_NORMAL;
    viewer->core.invert_normals = *((bool *) param);
  }

  IGL_INLINE void TW_CALL Viewer::get_invert_normals_cb(void *param,void *clientData)
  {
    *((bool *) param) = static_cast<Viewer *>(clientData)->core.invert_normals;
  }

  IGL_INLINE void TW_CALL Viewer::set_face_based_cb(const void *param,void *clientData)
  {
    Viewer *viewer = static_cast<Viewer *>(clientData);
    viewer->data.set_face_based(*((bool *) param));
  }

  IGL_INLINE void TW_CALL Viewer::get_face_based_cb(void *param,void *clientData)
  {
    *((bool *) param) = static_cast<Viewer *>(clientData)->data.face_based;
  }

  IGL_INLINE void TW_CALL Viewer::open_dialog_mesh(void *clientData)
  {
    std::string fname = igl::file_dialog_open();

    if (fname.length() == 0)
      return;

    static_cast<Viewer *>(clientData)->load_mesh_from_file(fname.c_str());
  }

  IGL_INLINE int Viewer::launch(std::string filename)
  {
    GLFWwindow* window;

    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
      return EXIT_FAILURE;

    glfwWindowHint(GLFW_SAMPLES, 4);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    #ifdef __APPLE__
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);
    #endif
    window = glfwCreateWindow(1280, 800, "libigl viewer", NULL, NULL);
    if (!window)
    {
      glfwTerminate();
      return EXIT_FAILURE;
    }

  glfwMakeContextCurrent(window);

#ifndef __APPLE__
  glewExperimental = true;
  GLenum err = glewInit();
  if (GLEW_OK != err)
  {
    /* Problem: glewInit failed, something is seriously wrong. */
    fprintf(stderr, "Error: %s\n", glewGetErrorString(err));
  }
  fprintf(stdout, "Status: Using GLEW %s\n", glewGetString(GLEW_VERSION));
#endif

    #ifdef DEBUG
      int major, minor, rev;
      major = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MAJOR);
      minor = glfwGetWindowAttrib(window, GLFW_CONTEXT_VERSION_MINOR);
      rev = glfwGetWindowAttrib(window, GLFW_CONTEXT_REVISION);
      printf("OpenGL version recieved: %d.%d.%d\n", major, minor, rev);
      printf("Supported OpenGL is %s\n", (const char*)glGetString(GL_VERSION));
      printf("Supported GLSL is %s\n", (const char*)glGetString(GL_SHADING_LANGUAGE_VERSION));
    #endif

    glfwSetInputMode(window,GLFW_CURSOR,GLFW_CURSOR_NORMAL);

    // Initialize AntTweakBar
    TwInit(TW_OPENGL_CORE, NULL);
    TwCopyStdStringToClientFunc(static_cast<TwCopyStdStringToClient>(::copy_str));


    // Initialize IGL viewer
    init();
    __viewer = this;

    // Register callbacks
    glfwSetKeyCallback(window, glfw_key_callback);
    glfwSetCursorPosCallback(window,glfw_mouse_move);
    glfwSetWindowSizeCallback(window,glfw_window_size);
    glfwSetMouseButtonCallback(window,glfw_mouse_press);
    glfwSetScrollCallback(window,glfw_mouse_scroll);
    glfwSetCharCallback(window, glfw_char_callback);

    // Handle retina displays (windows and mac)
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);

    int width_window, height_window;
    glfwGetWindowSize(window, &width_window, &height_window);

    highdpi = width/width_window;

    glfw_window_size(window,width_window,height_window);

    opengl.init();

    // Alec: It seems silly to overload launch to take a filename as an
    // argument. load_mesh_from_file has many side effects so it makes
    // debugging launch difficult.

    // Load the mesh passed as input
    if (filename.size() > 0)
    {
      load_mesh_from_file(filename.c_str());
    }

    core.align_camera_center(data.V,data.F);

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
    }

    opengl.free();
    core.shut();

    shutdown_plugins();

    glfwDestroyWindow(window);
    glfwTerminate();
    return EXIT_SUCCESS;
  }
} // end namespace
