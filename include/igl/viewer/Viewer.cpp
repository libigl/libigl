// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2014 Daniele Panozzo <daniele.panozzo@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.

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

#define GLFW_INCLUDE_GLU
#ifndef _WIN32
  #define GLFW_INCLUDE_GLCOREARB
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

#include <nanogui/formscreen.h>

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

#ifdef ENABLE_SERIALIZATION
#include <igl/serialize.h>
#endif

// Internal global variables used for glfw event handling
static igl::viewer::Viewer * __viewer;
static double highdpi = 1;
static double scroll_x = 0;
static double scroll_y = 0;
static int global_KMod = 0;

static void glfw_mouse_press(GLFWwindow* window, int button, int action, int modifier)
{
  bool tw_used = __viewer->ngui->mouseButtonEvent(window,button,action,modifier);

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
  if(! __viewer->ngui->charEvent(window,codepoint) )
  {
    __viewer->key_pressed(codepoint, modifier);
  }
}

static void glfw_key_callback(GLFWwindow* window, int key, int scancode, int action, int modifier)
{
  if (key == GLFW_KEY_ESCAPE && action == GLFW_PRESS)
    glfwSetWindowShouldClose(window, GL_TRUE);

  if (__viewer->ngui->keyEvent(window,key,scancode,action,modifier) == false)
  {
    if (action == GLFW_PRESS)
      __viewer->key_down(key, modifier);
    else if(action == GLFW_RELEASE)
      __viewer->key_up(key, modifier);
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
  if(__viewer->ngui->cursorPosEvent(window,x,y) == false || __viewer->down)
  {
    __viewer->mouse_move(x*highdpi, y*highdpi);
  }
}

static void glfw_mouse_scroll(GLFWwindow* window, double x, double y)
{
  using namespace std;
  scroll_x += x;
  scroll_y += y;

  if (__viewer->ngui->scrollEvent(window,x,y) == false)
    __viewer->mouse_scroll(y);
}

static void glfw_drop_callback(GLFWwindow *window,int count,const char **filenames)
{
  __viewer->ngui->dropEvent(window,count,filenames);
}

namespace igl
{
namespace viewer
{
  IGL_INLINE void Viewer::init()
  {
    using namespace nanogui;

    ngui->setInputCellSize(Eigen::Vector2i(60,20));

    // Create nanogui widgets
    ngui->addNewWindow(Eigen::Vector2i(10,10),"libIGL-Viewer");

    // ---------------------- LOADING ----------------------

  #ifdef ENABLE_SERIALIZATION
    ngui->addNewGroup("Workspace",FormScreen::Layout::Horizontal);
    ngui->addButton("Load",[&](){this->load_scene();});
    ngui->addButton("Save",[&](){this->save_scene();});
  #endif

  #ifdef ENABLE_IO
    ngui->addNewGroup("Mesh",FormScreen::Layout::Horizontal);
    ngui->addButton("Load",[&](){this->open_dialog_load_mesh();});
    ngui->addButton("Save",[&](){this->open_dialog_save_mesh();});
  #endif

    ngui->addNewGroup("Viewing Options",FormScreen::Layout::Vertical);
    ngui->addButton("Center object",[&](){this->core.align_camera_center(this->data.V,this->data.F);});
    ngui->addButton("Snap canonical view",[&]()
    {
      this->snap_to_canonical_quaternion();
    });
    ngui->addVariable(core.camera_zoom,"Zoom");
    ngui->addVariable(core.orthographic,"Orthographic view");

    ngui->addNewGroup("Draw options",FormScreen::Layout::Vertical);

    ngui->addVariable([&](bool checked)
    {
      this->data.set_face_based(checked);
    },[&]()
    {
      return this->data.face_based;
    }, "Face-based",false);

    ngui->addVariable(core.show_texture,"Show texture");

    ngui->addVariable([&](bool checked)
    {
      this->data.dirty |= ViewerData::DIRTY_NORMAL;
      this->core.invert_normals = checked;
    },[&]()
    {
      return this->core.invert_normals;
    },
      "Invert normals",false);

    ngui->addVariable(core.show_overlay,"Show overlay");
    ngui->addVariable(core.show_overlay_depth,"Show overlay depth");
    ngui->addColorPicker(core.background_color,"Background");
    ngui->addColorPicker(core.line_color,"Line color");
    ngui->addVariable(core.shininess,"Shininess");

    ngui->addNewGroup("Overlays");
    ngui->addVariable(core.show_lines,"Wireframe");
    ngui->addVariable(core.show_faces,"Fill");
    ngui->addVariable(core.show_vertid,"Show vertex labels");
    ngui->addVariable(core.show_faceid,"Show faces labels");

    ngui->layout();

    core.init();

    if (callback_init)
      if (callback_init(*this))
        return;

    init_plugins();
  }

  IGL_INLINE Viewer::Viewer()
  {
    ngui = nullptr;

    // Temporary variables initialization
    down = false;
    hack_never_moved = true;
    scroll_position = 0.0f;

    // Per face
    data.set_face_based(false);

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

  }

  IGL_INLINE void Viewer::init_plugins()
  {
    // Init all plugins
    for (unsigned int i = 0; i<plugins.size(); ++i)
      plugins[i]->init(this);
  }

  IGL_INLINE Viewer::~Viewer()
  {
    if(!ngui) delete ngui;
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

      if (!(
            igl::readOBJ(
              mesh_file_name_string,
              V, UV_V, corner_normals, F, UV_F, fNormIndices)))
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

  IGL_INLINE bool Viewer::key_pressed(unsigned int unicode_key,int modifiers)
  {
    if (callback_key_pressed)
      if (callback_key_pressed(*this,unicode_key,modifiers))
        return true;

    for (unsigned int i = 0; i<plugins.size(); ++i)
      if (plugins[i]->key_pressed(unicode_key, modifiers))
        return true;
    return false;
  }

  IGL_INLINE bool Viewer::key_down(int key,int modifiers)
  {
    if (callback_key_down)
      if (callback_key_down(*this,key,modifiers))
        return true;
    for (unsigned int i = 0; i<plugins.size(); ++i)
      if (plugins[i]->key_down(key, modifiers))
        return true;
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
      if (callback_mouse_down(*this,static_cast<int>(button),modifier))
        return true;

    for (unsigned int i = 0; i<plugins.size(); ++i)
      if(plugins[i]->mouse_down(static_cast<int>(button),modifier))
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

    Eigen::Vector3f coord = igl::project(Eigen::Vector3f(center(0),center(1),center(2)), (core.view * core.model).eval(), core.proj, core.viewport);
    down_mouse_z = coord[2];
    down_rotation = core.trackball_angle;

    mouse_mode = MouseMode::Rotation;

    switch (button)
    {
      case MouseButton::Left:
        mouse_mode = MouseMode::Rotation;
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
              igl::trackball(core.viewport(2),
                             core.viewport(3),
                             2.0f,
                             down_rotation,
                             down_mouse_x,
                             down_mouse_y,
                             mouse_x,
                             mouse_y,
                             core.trackball_angle);
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
          Eigen::Vector3f pos1 = igl::unproject(Eigen::Vector3f(mouse_x, core.viewport[3] - mouse_y, down_mouse_z), (core.view * core.model).eval(), core.proj, core.viewport);
          Eigen::Vector3f pos0 = igl::unproject(Eigen::Vector3f(down_mouse_x, core.viewport[3] - down_mouse_y, down_mouse_z), (core.view * core.model).eval(), core.proj, core.viewport);

          Eigen::Vector3f diff = pos1 - pos0;
          core.model_translation = down_translation + Eigen::Vector3f(diff[0],diff[1],diff[2]);

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

    ngui->draw();
  }

  IGL_INLINE bool Viewer::save_scene()
  {
    std::string fname = igl::file_dialog_save();
    if (fname.length() == 0)
      return false;

#ifdef ENABLE_SERIALIZATION

    igl::serialize(core,"Core",fname.c_str(),true);

#ifndef ENABLE_SERIALIZATION_CORE_ONLY
    igl::serialize(data,"Data",fname.c_str());
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
#ifdef ENABLE_SERIALIZATION

    igl::deserialize(core,"Core",fname.c_str());

#ifndef ENABLE_SERIALIZATION_CORE_ONLY
    igl::deserialize(data,"Data",fname.c_str());
    for(unsigned int i = 0; i <plugins.size(); ++i)
      igl::deserialize(*plugins[i],plugins[i]->plugin_name,fname.c_str());
#endif

#endif

    return true;
  }

  IGL_INLINE void Viewer::resize(int w,int h)
  {
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

  IGL_INLINE int Viewer::launch(bool resizable,bool fullscreen)
  {
    GLFWwindow* window;

    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
      return EXIT_FAILURE;

    glfwWindowHint(GLFW_SAMPLES, 8);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);

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
      window = glfwCreateWindow(1280,800,"libigl viewer",nullptr,nullptr);
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
    ngui = new nanogui::FormScreen();
    ngui->init(window);

    __viewer = this;

    // Register callbacks
    glfwSetKeyCallback(window, glfw_key_callback);
    glfwSetCursorPosCallback(window,glfw_mouse_move);
    glfwSetWindowSizeCallback(window,glfw_window_size);
    glfwSetMouseButtonCallback(window,glfw_mouse_press);
    glfwSetScrollCallback(window,glfw_mouse_scroll);
    glfwSetCharModsCallback(window,glfw_char_mods_callback);
    glfwSetDropCallback(window,glfw_drop_callback);

    // Handle retina displays (windows and mac)
    int width, height;
    glfwGetFramebufferSize(window, &width, &height);

    int width_window, height_window;
    glfwGetWindowSize(window, &width_window, &height_window);

    highdpi = width/width_window;

    glfw_window_size(window,width_window,height_window);

    opengl.init();

    core.align_camera_center(data.V,data.F);

    // Initialize IGL viewer
    init();

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
}
