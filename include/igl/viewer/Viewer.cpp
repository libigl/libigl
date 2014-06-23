#include "Viewer.h"

#ifdef _WIN32
#  include <windows.h>
#  undef max
#  undef min
#  include <GL/glew.h>
#endif

#ifdef __APPLE__
#   include <OpenGL/gl3.h>
#   define __gl_h_ /* Prevent inclusion of the old gl.h */
#else
#   ifdef _WIN32
#       include <windows.h>
#   endif
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

//OK NV
Eigen::Vector3f project(const Eigen::Vector3f&  obj,
                        const Eigen::Matrix4f& model,
                        const Eigen::Matrix4f& proj,
                        const Eigen::Vector4f&  viewport)
{
  Eigen::Vector4f tmp;
  tmp << obj,1;

  tmp = model * tmp;

  tmp = proj * tmp;

  tmp = tmp.array() / tmp(3);
  tmp = tmp.array() * 0.5f + 0.5f;
  tmp(0) = tmp(0) * viewport(2) + viewport(0);
  tmp(1) = tmp(1) * viewport(3) + viewport(1);

  return tmp.head(3);
}

Eigen::Vector3f unproject(const Eigen::Vector3f& win,
                          const Eigen::Matrix4f& model,
                          const Eigen::Matrix4f& proj,
                          const Eigen::Vector4f& viewport)
{
  Eigen::Matrix4f Inverse = (proj * model).inverse();

  Eigen::Vector4f tmp;
  tmp << win, 1;
  tmp(0) = (tmp(0) - viewport(0)) / viewport(2);
  tmp(1) = (tmp(1) - viewport(1)) / viewport(3);
  tmp = tmp.array() * 2.0f - 1.0f;

  Eigen::Vector4f obj = Inverse * tmp;
  obj /= obj(3);

  return obj.head(3);
}

Eigen::Matrix4f lookAt (
                        const Eigen::Vector3f& eye,
                        const Eigen::Vector3f& center,
                        const Eigen::Vector3f& up)
{
  Eigen::Vector3f f = (center - eye).normalized();
  Eigen::Vector3f s = f.cross(up).normalized();
  Eigen::Vector3f u = s.cross(f);

  Eigen::Matrix4f Result = Eigen::Matrix4f::Identity();
  Result(0,0) = s(0);
  Result(0,1) = s(1);
  Result(0,2) = s(2);
  Result(1,0) = u(0);
  Result(1,1) = u(1);
  Result(1,2) = u(2);
  Result(2,0) =-f(0);
  Result(2,1) =-f(1);
  Result(2,2) =-f(2);
  Result(0,3) =-s.transpose() * eye;
  Result(1,3) =-u.transpose() * eye;
  Result(2,3) = f.transpose() * eye;
  return Result;
}

Eigen::Matrix4f ortho (
                       const float left,
                       const float right,
                       const float bottom,
                       const float top,
                       const float zNear,
                       const float zFar
                       )
{
  Eigen::Matrix4f Result = Eigen::Matrix4f::Identity();
  Result(0,0) = 2.0f / (right - left);
  Result(1,1) = 2.0f / (top - bottom);
  Result(2,2) = - 2.0f / (zFar - zNear);
  Result(0,3) = - (right + left) / (right - left);
  Result(1,3) = - (top + bottom) / (top - bottom);
  Result(2,3) = - (zFar + zNear) / (zFar - zNear);
  return Result;
}

Eigen::Matrix4f frustum (
                         const float left,
                         const float right,
                         const float bottom,
                         const float top,
                         const float nearVal,
                         const float farVal)
{
  Eigen::Matrix4f Result = Eigen::Matrix4f::Zero();
  Result(0,0) = (2.0f * nearVal) / (right - left);
  Result(1,1) = (2.0f * nearVal) / (top - bottom);
  Result(0,2) = (right + left) / (right - left);
  Result(1,2) = (top + bottom) / (top - bottom);
  Result(2,2) = -(farVal + nearVal) / (farVal - nearVal);
  Result(3,2) = -1.0f;
  Result(2,3) = -(2.0f * farVal * nearVal) / (farVal - nearVal);
  return Result;
}

Eigen::Matrix4f scale (const Eigen::Matrix4f& m,
                       const Eigen::Vector3f& v)
{
  Eigen::Matrix4f Result;
  Result.col(0) = m.col(0).array() * v(0);
  Result.col(1) = m.col(1).array() * v(1);
  Result.col(2) = m.col(2).array() * v(2);
  Result.col(3) = m.col(3);
  return Result;
}

Eigen::Matrix4f translate(
                          const Eigen::Matrix4f& m,
                          const Eigen::Vector3f& v)
{
  Eigen::Matrix4f Result = m;
  Result.col(3) = m.col(0).array() * v(0) + m.col(1).array() * v(1) + m.col(2).array() * v(2) + m.col(3).array();
  return Result;
}


#include <limits>
#include <cassert>

#ifdef ENABLE_XML_SERIALIZATION
  #include "igl/xml/XMLSerializer.h"
#endif

#include <igl/readOBJ.h>
#include <igl/readOFF.h>
#include <igl/per_face_normals.h>
#include <igl/per_vertex_normals.h>
#include <igl/per_corner_normals.h>
#include <igl/vf.h>
#include <igl/adjacency_list.h>
#include <igl/writeOBJ.h>
#include <igl/writeOFF.h>
#include <igl/massmatrix.h>
#include <igl/file_dialog_open.h>
#include <igl/file_dialog_save.h>
#include <igl/quat_to_mat.h>
#include <igl/quat_mult.h>
#include <igl/axis_angle_to_quat.h>
#include <igl/trackball.h>
#include <igl/snap_to_canonical_view_quat.h>
#include <TwOpenGLCore.h>

// Plugin manager (exported to other compilation units)
igl::Plugin_manager igl_viewer_plugin_manager;

// Internal global variables used for glfw event handling
static igl::Viewer * __viewer;
static double highdpi = 1;
static double scroll_x = 0;
static double scroll_y = 0;

/* This class extends the font rendering code in AntTweakBar
   so that it can be used to render text at arbitrary 3D positions */
class TextRenderer : public CTwGraphOpenGLCore {
public:
  TextRenderer() : m_shaderHandleBackup(0) { }

  virtual int Init()
  {
    int retval = CTwGraphOpenGLCore::Init();
    if (retval == 1)
    {
      std::string vertexShader =
          "#version 150\n"
          "uniform vec2 offset;"
          "uniform vec2 wndSize;"
          "uniform vec4 color;"
          "uniform float depth;"
          "in vec2 vertex;"
          "in vec2 uv;"
          "out vec4 fcolor;"
          "out vec2 fuv;"
          "void main() {"
          "  gl_Position = vec4(2.0*(vertex.x+offset.x-0.5)/wndSize.x - 1.0,"
          "                     1.0 - 2.0*(vertex.y+offset.y-0.5)/wndSize.y,"
          "                     depth, 1);"
          " fuv = uv;"
          " fcolor = color;"
          "}";

      std::string fragmentShader =
        "#version 150\n"
        "uniform sampler2D tex;"
        "in vec2 fuv;"
        "in vec4 fcolor;"
        "out vec4 outColor;"
        "void main() { outColor.rgb = fcolor.bgr; outColor.a = fcolor.a * texture(tex, fuv).r; }";

      if (!m_shader.init(vertexShader, fragmentShader, "outColor"))
        return 0;

      /* Adjust location bindings */
      glBindAttribLocation(m_shader.program_shader, 0, "vertex");
      glBindAttribLocation(m_shader.program_shader, 1, "uv");
      glBindAttribLocation(m_shader.program_shader, 2, "color");
      glLinkProgram(m_shader.program_shader);

      m_shaderHandleBackup = m_TriTexUniProgram;
      m_TriTexUniProgram = m_shader.program_shader;
      m_TriTexUniLocationOffset = m_shader.uniform("offset");
      m_TriTexUniLocationWndSize = m_shader.uniform("wndSize");
      m_TriTexUniLocationColor = m_shader.uniform("color");
      m_TriTexUniLocationTexture = m_shader.uniform("tex");
      m_TriTexUniLocationDepth = m_shader.uniform("depth");
    }
    return retval;
  }

  virtual int Shut()
  {
    for (auto kv : m_textObjects)
      DeleteTextObj(kv.second);
    m_shader.free();
    m_TriTexUniProgram = m_shaderHandleBackup;
    return CTwGraphOpenGLCore::Shut();
  }

  void BeginDraw(const Eigen::Matrix4f &view, const Eigen::Matrix4f &proj,
    const Eigen::Vector4f &_viewport, float _object_scale)
  {
    viewport = _viewport;
    proj_matrix = proj;
    view_matrix = view;
    CTwGraphOpenGLCore::BeginDraw(viewport[2], viewport[3]);
    glEnable(GL_DEPTH_TEST);
    glDepthMask(GL_FALSE);
    object_scale = _object_scale;
  }

  void EndDraw()
  {
    /* Limit the number of cached text objects */
    for (auto it = m_textObjects.cbegin(); it != m_textObjects.cend(); )
    {
      if (m_textObjects.size() < 1000000)
        break;
      DeleteTextObj(it->second);
      m_textObjects.erase(it++);
    }

    glDepthMask(GL_TRUE);
    CTwGraphOpenGLCore::EndDraw();
  }

  void DrawText(Eigen::Vector3d pos, Eigen::Vector3d normal, const std::string &text)
  {
    pos += normal * 0.005f * object_scale;
    Eigen::Vector3f coord = project(Eigen::Vector3f(pos(0), pos(1), pos(2)),
        view_matrix, proj_matrix, viewport);
    auto it = m_textObjects.find(text);
    void *text_obj = nullptr;
    if (it == m_textObjects.end())
    {
      text_obj = NewTextObj();
      BuildText(text_obj, &text, NULL, NULL, 1, g_DefaultNormalFont, 0, 0);
      m_textObjects[text] = text_obj;
    } else {
      text_obj = it->second;
    }
    m_shader.bind();
    glUniform1f(m_TriTexUniLocationDepth, 2*(coord(2)-0.5f));
    CTwGraphOpenGLCore::DrawText(text_obj, coord[0], viewport[3] - coord[1], COLOR32_BLUE, 0);
  }
protected:
  igl::Viewer::OpenGL_shader m_shader;
  std::map<std::string, void *> m_textObjects;
  GLuint m_shaderHandleBackup;
  GLuint m_TriTexUniLocationDepth;
  Eigen::Matrix4f view_matrix, proj_matrix;
  Eigen::Vector4f viewport;
  float object_scale;
};

static TextRenderer __font_renderer;

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

int global_KMod = 0;

int TwEventKeyGLFW3(int glfwKey, int glfwAction)
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
  size[1] = height-pos[1]-10;
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
  void Viewer::init(Plugin_manager* pm)
  {
    plugin_manager = pm;

    // Create a tweak bar
    bar = TwNewBar("libIGL-Viewer");
    TwDefine(" libIGL-Viewer help='This is a simple 3D mesh viewer.' "); // Message added to the help bar->
    TwDefine(" libIGL-Viewer size='200 685'"); // change default tweak bar size
    TwDefine(" libIGL-Viewer color='76 76 127' "); // change default tweak bar color
    TwDefine(" libIGL-Viewer refresh=0.5"); // change refresh rate

    // ---------------------- LOADING ----------------------

    #ifdef ENABLE_XML_SERIALIZATION
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
    TwAddVarRW(bar, "Zoom", TW_TYPE_FLOAT, &(options.camera_zoom),
               " min=0.05 max=50 step=0.1 keyIncr=+ keyDecr=- help='Scale the object (1=original size).' group='Scene'");
    TwAddButton(bar,"SnapView", snap_to_canonical_quaternion_cb, this,
                " group='Scene'"
                " label='Snap to canonical view' key=Z "
                " help='Snaps view to nearest canonical view.'");
    TwAddVarRW(bar,"LightDir", TW_TYPE_DIR3F, options.light_position.data(),
               " group='Scene'"
               " label='Light direction' open help='Change the light direction.' ");

    // ---------------------- DRAW OPTIONS ----------------------
    TwAddVarRW(bar, "Toggle Orthographic/Perspective", TW_TYPE_BOOLCPP, &(options.orthographic),
               " group='Viewing Options'"
               " label='Orthographic view' "
               " help='Toggles orthographic / perspective view. Default: perspective.'");

    TwAddVarCB(bar,"Face-based Normals/Colors", TW_TYPE_BOOLCPP, set_face_based_cb, get_face_based_cb, this,
               " group='Draw options'"
               " label='Face-based' key=T help='Toggle per face shading/colors.' ");

    TwAddVarRW(bar,"Show texture", TW_TYPE_BOOLCPP, &(options.show_texture),
               " group='Draw options'");

    TwAddVarCB(bar,"Invert Normals", TW_TYPE_BOOLCPP, set_invert_normals_cb, get_invert_normals_cb, this,
               " group='Draw options'"
               " label='Invert normals' key=i help='Invert normal directions for inside out meshes.' ");

    TwAddVarRW(bar,"ShowOverlay", TW_TYPE_BOOLCPP, &(options.show_overlay),
               " group='Draw options'"
               " label='Show overlay' key=o help='Show the overlay layers.' ");
    TwAddVarRW(bar,"ShowOverlayDepth", TW_TYPE_BOOLCPP, &(options.show_overlay_depth),
               " group='Draw options'"
               " label='Show overlay depth test' help='Enable the depth test for overlay layer.' ");
    TwAddVarRW(bar,"Background color", TW_TYPE_COLOR3F,
               options.background_color.data(),
               " help='Select a background color' colormode=hls group='Draw options'");
    TwAddVarRW(bar, "LineColor", TW_TYPE_COLOR3F,
               options.line_color.data(),
               " label='Line color' help='Select a outline color' group='Draw options'");
    TwAddVarRW(bar,"Shininess",TW_TYPE_FLOAT,&options.shininess," group='Draw options'"
               " min=1 max=128");

    // ---------------------- Overlays ----------------------

    TwAddVarRW(bar,"Wireframe", TW_TYPE_BOOLCPP, &(options.show_lines),
               " group='Overlays'"
               " label='Wireframe' key=l help='Toggle wire frame of mesh'");
    TwAddVarRW(bar,"Fill", TW_TYPE_BOOLCPP, &(options.show_faces),
               " group='Overlays'"
               " label='Fill' key=t help='Display filled polygons of mesh'");
    TwAddVarRW(bar,"ShowVertexId", TW_TYPE_BOOLCPP, &(options.show_vertid),
               " group='Overlays'"
               " label='Show Vertex Labels' key=';' help='Toggle vertex indices'");
    TwAddVarRW(bar,"ShowFaceId", TW_TYPE_BOOLCPP, &(options.show_faceid),
               " group='Overlays'"
               " label='Show Faces Labels' key='CTRL+;' help='Toggle face"
               " indices'");

    __font_renderer.Init();

    init_plugins();
  }

  Viewer::Viewer()
  {
    plugin_manager = 0;

    // Default shininess
    options.shininess = 35.0f;

    // Default colors
    options.background_color << 0.3f, 0.3f, 0.5f;
    options.line_color << 0.0f, 0.0f, 0.0f;

    // Default lights settings
    options.light_position << 0.0f, -0.30f, -5.0f;

    // Default trackball
    options.trackball_angle << 0.0f, 0.0f, 0.0f, 1.0f;

    // Defalut model viewing parameters
    options.model_zoom = 1.0f;
    options.model_translation << 0,0,0;

    // Camera parameters
    options.camera_zoom = 1.0f;
    options.orthographic = false;
    options.camera_view_angle = 45.0;
    options.camera_dnear = 1.0;
    options.camera_dfar = 100.0;
    options.camera_eye << 0, 0, 5;
    options.camera_center << 0, 0, 0;
    options.camera_up << 0, 1, 0;

    // Default visualization options
    options.show_faces = true;
    options.show_lines = true;
    options.invert_normals = false;
    options.show_overlay = true;
    options.show_overlay_depth = true;
    options.show_vertid = false;
    options.show_faceid = false;
    options.show_texture = false;

    // Default point size / line width
    options.point_size = 15;
    options.line_width = 0.5f;

    // Temporary variables initialization
    down = false;
    scroll_position = 0.0f;

    // Per face
    set_face_based(false);

    // C-style callbacks
    callback_pre_draw     = 0;
    callback_post_draw    = 0;
    callback_mouse_down   = 0;
    callback_mouse_up     = 0;
    callback_mouse_move   = 0;
    callback_mouse_scroll = 0;
    callback_key_down     = 0;
    callback_key_up       = 0;

    callback_pre_draw_data = 0;
    callback_post_draw     = 0;
    callback_mouse_down    = 0;
    callback_mouse_up      = 0;
    callback_mouse_move    = 0;
    callback_mouse_scroll  = 0;
    callback_key_down      = 0;
    callback_key_up        = 0;
  }

  void Viewer::init_plugins()
  {
    // Init all plugins
    if (plugin_manager)
      for (unsigned int i = 0; i<plugin_manager->plugin_list.size(); ++i)
        plugin_manager->plugin_list[i]->init(this);
  }

  Viewer::~Viewer()
  {
  }

  void Viewer::shutdown_plugins()
  {
    if (plugin_manager)
      for (unsigned int i = 0; i<plugin_manager->plugin_list.size(); ++i)
        plugin_manager->plugin_list[i]->shutdown();
  }

  bool Viewer::load_mesh_from_file(const char* mesh_file_name)
  {
    std::string mesh_file_name_string = std::string(mesh_file_name);

    // first try to load it with a plugin
    if (plugin_manager)
      for (unsigned int i = 0; i<plugin_manager->plugin_list.size(); ++i)
        if (plugin_manager->plugin_list[i]->load(mesh_file_name_string))
          return true;

    clear_mesh();

    size_t last_dot = mesh_file_name_string.rfind('.');
    if (last_dot == std::string::npos)
    {
      printf("Error: No file extension found in %s\n",mesh_file_name);
      return false;
    }

    std::string extension = mesh_file_name_string.substr(last_dot+1);

    if (extension == "off" || extension =="OFF")
    {
      if (!igl::readOFF(mesh_file_name_string, data.V, data.F))
        return false;
    }
    else if (extension == "obj" || extension =="OBJ")
    {
      Eigen::MatrixXd corner_normals;
      Eigen::MatrixXi fNormIndices;

      Eigen::MatrixXd UV_V;
      Eigen::MatrixXi UV_F;

      if (!(igl::readOBJ(mesh_file_name_string, data.V, data.F, corner_normals, fNormIndices, UV_V, UV_F)))
        return false;
    }
    else
    {
      // unrecognized file type
      printf("Error: %s is not a recognized file type.\n",extension.c_str());
      return false;
    }

    compute_normals();
    uniform_colors(Eigen::Vector3d(51.0/255.0,43.0/255.0,33.3/255.0),
                   Eigen::Vector3d(255.0/255.0,228.0/255.0,58.0/255.0),
                   Eigen::Vector3d(255.0/255.0,235.0/255.0,80.0/255.0));
    if (data.V_uv.rows() == 0)
      grid_texture();

    align_camera_center();

    if (plugin_manager)
      for (unsigned int i = 0; i<plugin_manager->plugin_list.size(); ++i)
        if (plugin_manager->plugin_list[i]->post_load())
          return true;

    return true;
  }

  void Viewer::compute_normals()
  {
    igl::per_face_normals(data.V, data.F, data.F_normals);
    igl::per_vertex_normals(data.V, data.F, data.F_normals, data.V_normals);
    data.dirty |= DIRTY_NORMAL;
  }

  void Viewer::uniform_colors(Eigen::Vector3d ambient, Eigen::Vector3d diffuse, Eigen::Vector3d specular)
  {
    data.V_material_ambient.resize(data.V.rows(),3);
    data.V_material_diffuse.resize(data.V.rows(),3);
    data.V_material_specular.resize(data.V.rows(),3);

    for (unsigned i=0; i<data.V.rows();++i)
    {
      data.V_material_ambient.row(i) = ambient;
      data.V_material_diffuse.row(i) = diffuse;
      data.V_material_specular.row(i) = specular;
    }

    data.F_material_ambient.resize(data.F.rows(),3);
    data.F_material_diffuse.resize(data.F.rows(),3);
    data.F_material_specular.resize(data.F.rows(),3);

    for (unsigned i=0; i<data.F.rows();++i)
    {
      data.F_material_ambient.row(i) = ambient;
      data.F_material_diffuse.row(i) = diffuse;
      data.F_material_specular.row(i) = specular;
    }
    data.dirty |= DIRTY_SPECULAR | DIRTY_DIFFUSE | DIRTY_AMBIENT;
  }

  void Viewer::grid_texture()
  {
    if (data.V_uv.rows() == 0)
    {
      data.V_uv = data.V.block(0, 0, data.V.rows(), 2);
      data.V_uv.col(0) = data.V_uv.col(0).array() - data.V_uv.col(0).minCoeff();
      data.V_uv.col(0) = data.V_uv.col(0).array() / data.V_uv.col(0).maxCoeff();
      data.V_uv.col(1) = data.V_uv.col(1).array() - data.V_uv.col(1).minCoeff();
      data.V_uv.col(1) = data.V_uv.col(1).array() / data.V_uv.col(1).maxCoeff();
      data.V_uv = data.V_uv.array() * 10;
      data.dirty |= DIRTY_TEXTURE;
    }

    unsigned size = 128;
    unsigned size2 = size/2;
    data.texture_R.resize(size, size);
    for (unsigned i=0; i<size; ++i)
    {
      for (unsigned j=0; j<size; ++j)
      {
        data.texture_R(i,j) = 0;
        if ((i<size2 && j<size2) || (i>=size2 && j>=size2))
          data.texture_R(i,j) = 255;
      }
    }

    data.texture_G = data.texture_R;
    data.texture_B = data.texture_R;
    data.dirty |= DIRTY_TEXTURE;
  }

  bool Viewer::save_mesh_to_file(const char* mesh_file_name)
  {
    std::string mesh_file_name_string(mesh_file_name);

    // first try to load it with a plugin
    if (plugin_manager)
      for (unsigned int i = 0; i<plugin_manager->plugin_list.size(); ++i)
        if (plugin_manager->plugin_list[i]->save(mesh_file_name_string))
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

  void Viewer::clear_mesh()
  {
    data.V                       = Eigen::MatrixXd (0,3);
    data.F                       = Eigen::MatrixXi (0,3);

    data.F_material_ambient      = Eigen::MatrixXd (0,3);
    data.F_material_diffuse      = Eigen::MatrixXd (0,3);
    data.F_material_specular     = Eigen::MatrixXd (0,3);

    data.V_material_ambient      = Eigen::MatrixXd (0,3);
    data.V_material_diffuse      = Eigen::MatrixXd (0,3);
    data.V_material_specular     = Eigen::MatrixXd (0,3);

    data.F_normals               = Eigen::MatrixXd (0,3);
    data.V_normals               = Eigen::MatrixXd (0,3);

    data.V_uv                    = Eigen::MatrixXd (0,2);
    data.F_uv                    = Eigen::MatrixXi (0,3);

    data.lines                   = Eigen::MatrixXd (0,9);
    data.points                  = Eigen::MatrixXd (0,6);
    data.labels_positions        = Eigen::MatrixXd (0,3);
    data.labels_strings.clear();
  }

  bool Viewer::key_down(unsigned char key, int modifiers)
  {
    if (callback_key_down)
      if (callback_key_down(*this,key,modifiers))
        return true;

    if (plugin_manager)
      for (unsigned int i = 0; i <plugin_manager->plugin_list.size(); ++i)
        if (plugin_manager->plugin_list[i]->key_down(key, modifiers))
          return true;

    if (key == 'S')
      mouse_scroll(1);

    if (key == 'A')
      mouse_scroll(-1);

    // Why aren't these handled view AntTweakBar?
    if (key == 'z') // Don't use 'Z' because that clobbers snap_to_canonical_view_quat
      options.trackball_angle << 0.0f, 0.0f, 0.0f, 1.0f;

    if (key == 'y')
      options.trackball_angle << -sqrt(2.0f)/2.0f, 0.0f, 0.0f, sqrt(2.0f)/2.0f;

    if (key == 'x')
      options.trackball_angle << -0.5f, -0.5f, -0.5f, 0.5f;


    return false;
  }

  bool Viewer::key_up(unsigned char key, int modifiers)
  {
    if (callback_key_up)
      if (callback_key_up(*this,key,modifiers))
        return true;

    if (plugin_manager)
      for (unsigned int i = 0; i <plugin_manager->plugin_list.size(); ++i)
        if (plugin_manager->plugin_list[i]->key_up(key, modifiers))
          return true;

    return false;
  }

  bool Viewer::mouse_down(MouseButton button, int modifier)
  {
    if (callback_mouse_down)
      if (callback_mouse_down(*this,button,modifier))
        return true;

    if (plugin_manager)
      for (unsigned int i = 0; i < plugin_manager->plugin_list.size(); ++i)
        if (plugin_manager->plugin_list[i]->mouse_down(button,modifier))
          return true;

    down = true;

    down_mouse_x = current_mouse_x;
    down_mouse_y = current_mouse_y;
    down_translation = options.model_translation;


    // Initialization code for the trackball
    Eigen::RowVector3d center;
    if (data.V.rows() == 0)
      center << 0,0,0;
    else
      center = data.V.colwise().sum()/data.V.rows();

    Eigen::Vector3f coord = project(Eigen::Vector3f(center(0),center(1),center(2)), view * model, proj, viewport);
    down_mouse_z = coord[2];
    down_rotation = options.trackball_angle;

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

  bool Viewer::mouse_up(MouseButton button, int modifier)
  {
    down = false;

    if (callback_mouse_up)
      if (callback_mouse_up(*this,button,modifier))
        return true;

    if (plugin_manager)
      for (unsigned int i = 0; i <plugin_manager->plugin_list.size(); ++i)
        if (plugin_manager->plugin_list[i]->mouse_up(button,modifier))
          return true;

    mouse_mode = NOTHING;

    return true;
  }

  bool Viewer::mouse_move(int mouse_x, int mouse_y)
  {
    current_mouse_x = mouse_x;
    current_mouse_y = mouse_y;

    if (callback_mouse_move)
      if (callback_mouse_move(*this,mouse_x,mouse_y))
        return true;

    if (plugin_manager)
      for (unsigned int i = 0; i < plugin_manager->plugin_list.size(); ++i)
        if (plugin_manager->plugin_list[i]->mouse_move(mouse_x, mouse_y))
          return true;

    if (down)
    {
      switch (mouse_mode)
      {
        case ROTATION :
        {
          igl::trackball(width,
                         height,
                         2.0f,
                         down_rotation.data(),
                         down_mouse_x,
                         down_mouse_y,
                         mouse_x,
                         mouse_y,
                         options.trackball_angle.data());
          //Eigen::Vector4f snapq = options.trackball_angle;

          break;
        }

        case TRANSLATE:
        {
          //translation
          Eigen::Vector3f pos1 = unproject(Eigen::Vector3f(mouse_x, viewport[3] - mouse_y, down_mouse_z), view * model, proj, viewport);
          Eigen::Vector3f pos0 = unproject(Eigen::Vector3f(down_mouse_x, viewport[3] - down_mouse_y, down_mouse_z), view * model, proj, viewport);

          Eigen::Vector3f diff = pos1 - pos0;
          options.model_translation = down_translation + Eigen::Vector3f(diff[0],diff[1],diff[2]);

          break;
        }
        case ZOOM:
        {
          //float delta = 0.001f * (mouse_x - down_mouse_x + mouse_y - down_mouse_y);
          float delta = 0.001f * (mouse_x - down_mouse_x + mouse_y - down_mouse_y);
          options.camera_zoom *= 1 + delta;
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

  bool Viewer::mouse_scroll(float delta_y)
  {
    scroll_position += delta_y;

    if (callback_mouse_scroll)
      if (callback_mouse_scroll(*this,delta_y))
        return true;

    if (plugin_manager)
      for (unsigned int i = 0; i <plugin_manager->plugin_list.size(); ++i)
        if (plugin_manager->plugin_list[i]->mouse_scroll(delta_y))
          return true;

    // Only zoom if there's actually a change
    if(delta_y != 0)
    {
      float mult = (1.0+((delta_y>0)?1.:-1.)*0.05);
      const float min_zoom = 0.1f;
      options.camera_zoom = (options.camera_zoom * mult > min_zoom ? options.camera_zoom * mult : min_zoom);
    }
    return true;
  }

  static GLuint create_shader_helper(GLint type, const std::string &shader_string)
  {
    using namespace std;
    if (shader_string.empty())
      return (GLuint) 0;

    GLuint id = glCreateShader(type);
    const char *shader_string_const = shader_string.c_str();
    glShaderSource(id, 1, &shader_string_const, NULL);
    glCompileShader(id);

    GLint status;
    glGetShaderiv(id, GL_COMPILE_STATUS, &status);

    if (status != GL_TRUE)
    {
      char buffer[512];
      if (type == GL_VERTEX_SHADER)
        cerr << "Vertex shader:" << endl;
      else if (type == GL_FRAGMENT_SHADER)
        cerr << "Fragment shader:" << endl;
      else if (type == GL_GEOMETRY_SHADER)
        cerr << "Geometry shader:" << endl;
      cerr << shader_string << endl << endl;
      glGetShaderInfoLog(id, 512, NULL, buffer);
      cerr << "Error: " << endl << buffer << endl;
      return (GLuint) 0;
    }

    return id;
  }

  bool Viewer::OpenGL_shader::init_from_files(
    const std::string &vertex_shader_filename,
    const std::string &fragment_shader_filename,
    const std::string &fragment_data_name,
    const std::string &geometry_shader_filename,
    int geometry_shader_max_vertices)
  {
    auto file_to_string = [](const std::string &filename)
    {
      std::ifstream t(filename);
      return std::string((std::istreambuf_iterator<char>(t)),
                          std::istreambuf_iterator<char>());
    };

    return init(
      file_to_string(vertex_shader_filename),
      file_to_string(fragment_shader_filename),
      fragment_data_name,
      file_to_string(geometry_shader_filename),
      geometry_shader_max_vertices
   );
  }

  bool Viewer::OpenGL_shader::init(
    const std::string &vertex_shader_string,
    const std::string &fragment_shader_string,
    const std::string &fragment_data_name,
    const std::string &geometry_shader_string,
    int geometry_shader_max_vertices)
  {
    using namespace std;
    vertex_shader = create_shader_helper(GL_VERTEX_SHADER, vertex_shader_string);
    geometry_shader = create_shader_helper(GL_GEOMETRY_SHADER, geometry_shader_string);
    fragment_shader = create_shader_helper(GL_FRAGMENT_SHADER, fragment_shader_string);

    if (!vertex_shader || !fragment_shader)
      return false;

    program_shader = glCreateProgram();

    glAttachShader(program_shader, vertex_shader);
    glAttachShader(program_shader, fragment_shader);

    if (geometry_shader)
    {
      glAttachShader(program_shader, geometry_shader);

      /* This covers only basic cases and may need to be modified */
      glProgramParameteri(program_shader, GL_GEOMETRY_INPUT_TYPE, GL_TRIANGLES);
      glProgramParameteri(program_shader, GL_GEOMETRY_OUTPUT_TYPE, GL_TRIANGLES);
      glProgramParameteri(program_shader, GL_GEOMETRY_VERTICES_OUT, geometry_shader_max_vertices);
    }

    glBindFragDataLocation(program_shader, 0, fragment_data_name.c_str());
    glLinkProgram(program_shader);

    GLint status;
    glGetProgramiv(program_shader, GL_LINK_STATUS, &status);

    if (status != GL_TRUE)
    {
      char buffer[512];
      glGetProgramInfoLog(program_shader, 512, NULL, buffer);
      cerr << "Linker error: " << endl << buffer << endl;
      program_shader = 0;
      return false;
    }

    return true;
  }

  void Viewer::OpenGL_shader::bind()
  {
    glUseProgram(program_shader);
  }

  GLint Viewer::OpenGL_shader::attrib(const std::string &name) const
  {
    return glGetAttribLocation(program_shader, name.c_str());
  }

  GLint Viewer::OpenGL_shader::uniform(const std::string &name) const
  {
    return glGetUniformLocation(program_shader, name.c_str());
  }

  GLint Viewer::OpenGL_shader::bindVertexAttribArray(
    const std::string &name, GLuint bufferID, const Eigen::MatrixXf &M, bool refresh) const
  {
    GLint id = attrib(name);
    if (id < 0)
      return id;
    if (M.size() == 0)
    {
      glDisableVertexAttribArray(id);
      return id;
    }
    glBindBuffer(GL_ARRAY_BUFFER, bufferID);
    if (refresh)
      glBufferData(GL_ARRAY_BUFFER, sizeof(float)*M.size(), M.data(), GL_DYNAMIC_DRAW);
    glVertexAttribPointer(id, M.rows(), GL_FLOAT, GL_FALSE, 0, 0);
    glEnableVertexAttribArray(id);
    return id;
  }

  void Viewer::OpenGL_shader::free()
  {
    if (program_shader)
    {
      glDeleteProgram(program_shader);
      program_shader = 0;
    }
    if (vertex_shader)
    {
      glDeleteShader(vertex_shader);
      vertex_shader = 0;
    }
    if (fragment_shader)
    {
      glDeleteShader(fragment_shader);
      fragment_shader = 0;
    }
    if (geometry_shader)
    {
      glDeleteShader(geometry_shader);
      geometry_shader = 0;
    }
  }

  void Viewer::OpenGL_state::init()
  {
    // Mesh: Vertex Array Object & Buffer objects
    glGenVertexArrays(1, &vao_mesh);
    glBindVertexArray(vao_mesh);
    glGenBuffers(1, &vbo_V);
    glGenBuffers(1, &vbo_V_normals);
    glGenBuffers(1, &vbo_V_ambient);
    glGenBuffers(1, &vbo_V_diffuse);
    glGenBuffers(1, &vbo_V_specular);
    glGenBuffers(1, &vbo_V_uv);
    glGenBuffers(1, &vbo_F);
    glGenTextures(1, &vbo_tex);

    // Line overlay
    glGenVertexArrays(1, &vao_overlay_lines);
    glBindVertexArray(vao_overlay_lines);
    glGenBuffers(1, &vbo_lines_F);
    glGenBuffers(1, &vbo_lines_V);
    glGenBuffers(1, &vbo_lines_V_colors);

    // Point overlay
    glGenVertexArrays(1, &vao_overlay_points);
    glBindVertexArray(vao_overlay_points);
    glGenBuffers(1, &vbo_points_F);
    glGenBuffers(1, &vbo_points_V);
    glGenBuffers(1, &vbo_points_V_colors);

    dirty = DIRTY_ALL;
  }

  void Viewer::OpenGL_state::free()
  {
    glDeleteVertexArrays(1, &vao_mesh);
    glDeleteVertexArrays(1, &vao_overlay_lines);
    glDeleteVertexArrays(1, &vao_overlay_points);

    glDeleteBuffers(1, &vbo_V);
    glDeleteBuffers(1, &vbo_V_normals);
    glDeleteBuffers(1, &vbo_V_ambient);
    glDeleteBuffers(1, &vbo_V_diffuse);
    glDeleteBuffers(1, &vbo_V_specular);
    glDeleteBuffers(1, &vbo_V_uv);
    glDeleteBuffers(1, &vbo_F);
    glDeleteBuffers(1, &vbo_lines_F);
    glDeleteBuffers(1, &vbo_lines_V);
    glDeleteBuffers(1, &vbo_lines_V_colors);
    glDeleteBuffers(1, &vbo_points_F);
    glDeleteBuffers(1, &vbo_points_V);
    glDeleteBuffers(1, &vbo_points_V_colors);

    glDeleteTextures(1, &vbo_tex);
  }

  void Viewer::OpenGL_state::set_data(const Data &data, bool face_based, bool invert_normals)
  {
    bool per_corner_uv = (data.F_uv.rows() == data.F.rows());
    bool per_corner_normals = (data.F_normals.rows() == 3 * data.F.rows());

    dirty |= data.dirty;

    if (!face_based)
    {
      if (!per_corner_uv)
      {
        // Vertex positions
        if (dirty & DIRTY_POSITION)
          V_vbo = (data.V.transpose()).cast<float>();

        // Vertex normals
        if (dirty & DIRTY_NORMAL)
        {
          V_normals_vbo = (data.V_normals.transpose()).cast<float>();
          if (invert_normals)
            V_normals_vbo = -V_normals_vbo;
        }

        // Per-vertex material settings
        if (dirty & DIRTY_AMBIENT)
          V_ambient_vbo = (data.V_material_ambient.transpose()).cast<float>();
        if (dirty & DIRTY_DIFFUSE)
          V_diffuse_vbo = (data.V_material_diffuse.transpose()).cast<float>();
        if (dirty & DIRTY_SPECULAR)
          V_specular_vbo = (data.V_material_specular.transpose()).cast<float>();

        // Face indices
        if (dirty & DIRTY_FACE)
          F_vbo = (data.F.transpose()).cast<unsigned>();

        // Texture coordinates
        if (dirty & DIRTY_UV)
          V_uv_vbo = (data.V_uv.transpose()).cast<float>();
      }
      else
      {
        // Per vertex properties with per corner UVs
        if (dirty & DIRTY_POSITION)
        {
          V_vbo.resize(3,data.F.rows()*3);
          for (unsigned i=0; i<data.F.rows();++i)
            for (unsigned j=0;j<3;++j)
              V_vbo.col(i*3+j) = data.V.row(data.F(i,j)).transpose().cast<float>();
        }

        if (dirty & DIRTY_AMBIENT)
        {
          V_ambient_vbo.resize(3,data.F.rows()*3);
          for (unsigned i=0; i<data.F.rows();++i)
            for (unsigned j=0;j<3;++j)
              V_ambient_vbo.col (i*3+j) = data.V_material_ambient.row(data.F(i,j)).transpose().cast<float>();
        }

        if (dirty & DIRTY_DIFFUSE)
        {
          V_diffuse_vbo.resize(3,data.F.rows()*3);
          for (unsigned i=0; i<data.F.rows();++i)
            for (unsigned j=0;j<3;++j)
              V_diffuse_vbo.col (i*3+j) = data.V_material_diffuse.row(data.F(i,j)).transpose().cast<float>();
        }

        if (dirty & DIRTY_SPECULAR)
        {
          V_specular_vbo.resize(3,data.F.rows()*3);
          for (unsigned i=0; i<data.F.rows();++i)
            for (unsigned j=0;j<3;++j)
              V_specular_vbo.col(i*3+j) = data.V_material_specular.row(data.F(i,j)).transpose().cast<float>();
        }

        if (dirty & DIRTY_NORMAL)
        {
          V_normals_vbo.resize(3,data.F.rows()*3);
          for (unsigned i=0; i<data.F.rows();++i)
            for (unsigned j=0;j<3;++j)
              V_normals_vbo.col (i*3+j) = data.V_normals.row(data.F(i,j)).transpose().cast<float>();

          if (invert_normals)
            V_normals_vbo = -V_normals_vbo;
        }

        if (dirty & DIRTY_FACE)
        {
          F_vbo.resize(3,data.F.rows());
          for (unsigned i=0; i<data.F.rows();++i)
            F_vbo.col(i) << i*3+0, i*3+1, i*3+2;
        }

        if (dirty & DIRTY_UV)
        {
          V_uv_vbo.resize(2,data.F.rows()*3);
          for (unsigned i=0; i<data.F.rows();++i)
            for (unsigned j=0;j<3;++j)
              V_uv_vbo.col(i*3+j) = data.V_uv.row(data.F(i,j)).transpose().cast<float>();
        }
      }
    }
    else
    {
      if (dirty & DIRTY_POSITION)
      {
        V_vbo.resize(3,data.F.rows()*3);
        for (unsigned i=0; i<data.F.rows();++i)
          for (unsigned j=0;j<3;++j)
            V_vbo.col(i*3+j) = data.V.row(data.F(i,j)).transpose().cast<float>();
      }

      if (dirty & DIRTY_AMBIENT)
      {
        V_ambient_vbo.resize(3,data.F.rows()*3);
        for (unsigned i=0; i<data.F.rows();++i)
          for (unsigned j=0;j<3;++j)
            V_ambient_vbo.col (i*3+j) = data.F_material_ambient.row(i).transpose().cast<float>();
      }

      if (dirty & DIRTY_DIFFUSE)
      {
        V_diffuse_vbo.resize(3,data.F.rows()*3);
        for (unsigned i=0; i<data.F.rows();++i)
          for (unsigned j=0;j<3;++j)
            V_diffuse_vbo.col (i*3+j) = data.F_material_diffuse.row(i).transpose().cast<float>();
      }

      if (dirty & DIRTY_SPECULAR)
      {
        V_specular_vbo.resize(3,data.F.rows()*3);
        for (unsigned i=0; i<data.F.rows();++i)
          for (unsigned j=0;j<3;++j)
            V_specular_vbo.col(i*3+j) = data.F_material_specular.row(i).transpose().cast<float>();
      }

      if (dirty & DIRTY_NORMAL)
      {
        V_normals_vbo.resize(3,data.F.rows()*3);
        for (unsigned i=0; i<data.F.rows();++i)
          for (unsigned j=0;j<3;++j)
            V_normals_vbo.col (i*3+j) =
               per_corner_normals ?
                 data.F_normals.row(i*3+j).transpose().cast<float>() :
                 data.F_normals.row(i).transpose().cast<float>();

        if (invert_normals)
          V_normals_vbo = -V_normals_vbo;
      }

      if (dirty & DIRTY_FACE)
      {
        F_vbo.resize(3,data.F.rows());
        for (unsigned i=0; i<data.F.rows();++i)
          F_vbo.col(i) << i*3+0, i*3+1, i*3+2;
      }

      if (dirty & DIRTY_UV)
      {
          V_uv_vbo.resize(2,data.F.rows()*3);
          for (unsigned i=0; i<data.F.rows();++i)
            for (unsigned j=0;j<3;++j)
              V_uv_vbo.col(i*3+j) = data.V_uv.row(per_corner_uv ? data.F_uv(i,j) : data.F(i,j)).transpose().cast<float>();
      }
    }

    if (dirty & DIRTY_TEXTURE)
    {
      tex_u = data.texture_R.rows();
      tex_v = data.texture_R.cols();
      tex.resize(data.texture_R.size()*3);
      for (unsigned i=0;i<data.texture_R.size();++i)
      {
        tex(i*3+0) = data.texture_R(i);
        tex(i*3+1) = data.texture_G(i);
        tex(i*3+2) = data.texture_B(i);
      }
    }

    if (dirty & DIRTY_OVERLAY_LINES)
    {
      lines_V_vbo.resize(3, data.lines.rows()*2);
      lines_V_colors_vbo.resize(3, data.lines.rows()*2);
      lines_F_vbo.resize(1, data.lines.rows()*2);
      for (unsigned i=0; i<data.lines.rows();++i)
      {
        lines_V_vbo.col(2*i+0) = data.lines.block<1, 3>(i, 0).transpose().cast<float>();
        lines_V_vbo.col(2*i+1) = data.lines.block<1, 3>(i, 3).transpose().cast<float>();
        lines_V_colors_vbo.col(2*i+0) = data.lines.block<1, 3>(i, 6).transpose().cast<float>();
        lines_V_colors_vbo.col(2*i+1) = data.lines.block<1, 3>(i, 6).transpose().cast<float>();
        lines_F_vbo(2*i+0) = 2*i+0;
        lines_F_vbo(2*i+1) = 2*i+1;
      }
    }

    if (dirty & DIRTY_OVERLAY_POINTS)
    {
      points_V_vbo.resize(3, data.points.rows());
      points_V_colors_vbo.resize(3, data.points.rows());
      points_F_vbo.resize(1, data.points.rows());
      for (unsigned i=0; i<data.points.rows();++i)
      {
        points_V_vbo.col(i) = data.points.block<1, 3>(i, 0).transpose().cast<float>();
        points_V_colors_vbo.col(i) = data.points.block<1, 3>(i, 3).transpose().cast<float>();
        points_F_vbo(i) = i;
      }
    }
  }

  void Viewer::OpenGL_state::bind_mesh()
  {
    glBindVertexArray(vao_mesh);
    shader_mesh.bind();
    shader_mesh.bindVertexAttribArray("position", vbo_V, V_vbo, dirty & DIRTY_POSITION);
    shader_mesh.bindVertexAttribArray("normal", vbo_V_normals, V_normals_vbo, dirty & DIRTY_NORMAL);
    shader_mesh.bindVertexAttribArray("Ka", vbo_V_ambient, V_ambient_vbo, dirty & DIRTY_AMBIENT);
    shader_mesh.bindVertexAttribArray("Kd", vbo_V_diffuse, V_diffuse_vbo, dirty & DIRTY_DIFFUSE);
    shader_mesh.bindVertexAttribArray("Ks", vbo_V_specular, V_specular_vbo, dirty & DIRTY_SPECULAR);
    shader_mesh.bindVertexAttribArray("texcoord", vbo_V_uv, V_uv_vbo, dirty & DIRTY_UV);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_F);
    if (dirty & DIRTY_FACE)
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned)*F_vbo.size(), F_vbo.data(), GL_DYNAMIC_DRAW);

    glActiveTexture(GL_TEXTURE0);
    glBindTexture(GL_TEXTURE_2D, vbo_tex);
    if (dirty & DIRTY_TEXTURE)
    {
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
      glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
      glTexImage2D(GL_TEXTURE_2D, 0, GL_RGB, tex_u, tex_v, 0, GL_RGB, GL_UNSIGNED_BYTE, tex.data());
    }
    glUniform1i(shader_mesh.uniform("tex"), 0);
    dirty &= ~DIRTY_MESH;
  }

  void Viewer::OpenGL_state::bind_overlay_lines()
  {
    bool is_dirty = dirty & DIRTY_OVERLAY_LINES;

    glBindVertexArray(vao_overlay_lines);
    shader_overlay_lines.bind();
    shader_overlay_lines.bindVertexAttribArray("position", vbo_lines_V, lines_V_vbo, is_dirty);
    shader_overlay_lines.bindVertexAttribArray("color", vbo_lines_V_colors, lines_V_colors_vbo, is_dirty);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_lines_F);
    if (is_dirty)
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned)*lines_F_vbo.size(), lines_F_vbo.data(), GL_DYNAMIC_DRAW);

    dirty &= ~DIRTY_OVERLAY_LINES;
  }

  void Viewer::OpenGL_state::bind_overlay_points()
  {
    bool is_dirty = dirty & DIRTY_OVERLAY_POINTS;

    glBindVertexArray(vao_overlay_points);
    shader_overlay_points.bind();
    shader_overlay_points.bindVertexAttribArray("position", vbo_points_V, points_V_vbo, is_dirty);
    shader_overlay_points.bindVertexAttribArray("color", vbo_points_V_colors, points_V_colors_vbo, is_dirty);

    glBindBuffer(GL_ELEMENT_ARRAY_BUFFER, vbo_points_F);
    if (is_dirty)
      glBufferData(GL_ELEMENT_ARRAY_BUFFER, sizeof(unsigned)*points_F_vbo.size(), points_F_vbo.data(), GL_DYNAMIC_DRAW);

    dirty &= ~DIRTY_OVERLAY_POINTS;
  }

  void Viewer::OpenGL_state::draw_mesh(bool solid)
  {
    glPolygonMode(GL_FRONT_AND_BACK, solid ? GL_FILL : GL_LINE);

    /* Avoid Z-buffer fighting between filled triangles & wireframe lines */
    if (solid)
    {
      glEnable(GL_POLYGON_OFFSET_FILL);
      glPolygonOffset(1.0, 1.0);
    }
    glDrawElements(GL_TRIANGLES, 3*F_vbo.cols(), GL_UNSIGNED_INT, 0);

    glDisable(GL_POLYGON_OFFSET_FILL);
    glPolygonMode(GL_FRONT_AND_BACK, GL_FILL);
  }

  void Viewer::OpenGL_state::draw_overlay_lines()
  {
    glDrawElements(GL_LINES, lines_F_vbo.cols(), GL_UNSIGNED_INT, 0);
  }

  void Viewer::OpenGL_state::draw_overlay_points()
  {
    glDrawElements(GL_POINTS, points_F_vbo.cols(), GL_UNSIGNED_INT, 0);
  }

  void Viewer::init_opengl()
  {
    std::string mesh_vertex_shader_string =
    "#version 150\n"
    "uniform mat4 model;"
    "uniform mat4 view;"
    "uniform mat4 proj;"
    "in vec3 position;"
    "in vec3 normal;"
    "out vec3 position_eye;"
    "out vec3 normal_eye;"
    "in vec3 Ka;"
    "in vec3 Kd;"
    "in vec3 Ks;"
    "in vec2 texcoord;"
    "out vec2 texcoordi;"
    "out vec3 Kai;"
    "out vec3 Kdi;"
    "out vec3 Ksi;"

    "void main()"
    "{"
    "  position_eye = vec3 (view * model * vec4 (position, 1.0));"
    "  normal_eye = vec3 (view * model * vec4 (normal, 0.0));"
    "  normal_eye = normalize(normal_eye);"
    "  gl_Position = proj * vec4 (position_eye, 1.0);" //proj * view * model * vec4(position, 1.0);"
    "  Kai = Ka;"
    "  Kdi = Kd;"
    "  Ksi = Ks;"
    "  texcoordi = texcoord;"
    "}";

    std::string mesh_fragment_shader_string =
    "#version 150\n"
    "uniform mat4 model;"
    "uniform mat4 view;"
    "uniform mat4 proj;"
    "uniform vec4 fixed_color;"
    "in vec3 position_eye;"
    "in vec3 normal_eye;"
    "uniform vec3 light_position_world;"
    "vec3 Ls = vec3 (1, 1, 1);"
    "vec3 Ld = vec3 (1, 1, 1);"
    "vec3 La = vec3 (1, 1, 1);"
    "in vec3 Ksi;"
    "in vec3 Kdi;"
    "in vec3 Kai;"
    "in vec2 texcoordi;"
    "uniform sampler2D tex;"
    "uniform float specular_exponent;"
    "uniform float lighting_factor;"
    "uniform float texture_factor;"
    "out vec4 outColor;"
    "void main()"
    "{"
    "vec3 Ia = La * Kai;"    // ambient intensity

    "vec3 light_position_eye = vec3 (view * vec4 (light_position_world, 1.0));"
    "vec3 distance_to_light_eye = light_position_eye - position_eye;"
    "vec3 direction_to_light_eye = normalize (distance_to_light_eye);"
    "float dot_prod = dot (direction_to_light_eye, normal_eye);"
    "dot_prod = max (dot_prod, 0.0);"
    "vec3 Id = Ld * Kdi * dot_prod;"    // Diffuse intensity

    "vec3 reflection_eye = reflect (-direction_to_light_eye, normal_eye);"
    "vec3 surface_to_viewer_eye = normalize (-position_eye);"
    "float dot_prod_specular = dot (reflection_eye, surface_to_viewer_eye);"
    "dot_prod_specular = max (dot_prod_specular, 0.0);"
    "float specular_factor = pow (dot_prod_specular, specular_exponent);"
    "vec3 Is = Ls * Ksi * specular_factor;"    // specular intensity
    "vec4 color = vec4(lighting_factor * (Is + Id) + Ia, 1.0);"
    "outColor = mix(vec4(1,1,1,1), texture(tex, texcoordi), texture_factor) * color;"
    "if (fixed_color != vec4(0.0)) outColor = fixed_color;"
    "}";

    std::string overlay_vertex_shader_string =
    "#version 150\n"
    "uniform mat4 model;"
    "uniform mat4 view;"
    "uniform mat4 proj;"
    "in vec3 position;"
    "in vec3 color;"
    "out vec3 color_frag;"

    "void main()"
    "{"
    "  gl_Position = proj * view * model * vec4 (position, 1.0);"
    "  color_frag = color;"
    "}";

    std::string overlay_fragment_shader_string =
    "#version 150\n"
    "in vec3 color_frag;"
    "out vec4 outColor;"
    "void main()"
    "{"
    "  outColor = vec4(color_frag, 1.0);"
    "}";

    std::string overlay_point_fragment_shader_string =
    "#version 150\n"
    "in vec3 color_frag;"
    "out vec4 outColor;"
    "void main()"
    "{"
    "  if (length(gl_PointCoord - vec2(0.5)) > 0.5)"
    "    discard;"
    "  outColor = vec4(color_frag, 1.0);"
    "}";

    opengl.init();
    opengl.shader_mesh.init(mesh_vertex_shader_string,
        mesh_fragment_shader_string, "outColor");
    opengl.shader_overlay_lines.init(overlay_vertex_shader_string,
        overlay_fragment_shader_string, "outColor");
    opengl.shader_overlay_points.init(overlay_vertex_shader_string,
        overlay_point_fragment_shader_string, "outColor");
  }

  void Viewer::free_opengl()
  {
    opengl.shader_mesh.free();
    opengl.shader_overlay_lines.free();
    opengl.shader_overlay_points.free();
    opengl.free();
    __font_renderer.Shut();
  }

  void Viewer::draw()
  {
    using namespace std;
    using namespace Eigen;
    glClearColor(options.background_color[0],
                 options.background_color[1],
                 options.background_color[2],
                 1.0f);
    glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

    glEnable(GL_DEPTH_TEST);

    if (callback_pre_draw)
      if (callback_pre_draw(*this))
        return;

    if (plugin_manager)
      for (unsigned int i = 0; i <plugin_manager->plugin_list.size(); ++i)
        if (plugin_manager->plugin_list[i]->pre_draw())
          return;

    /* Bind and potentially refresh mesh/line/point data */
    if (data.dirty)
    {
      opengl.set_data(data, options.face_based, options.invert_normals);
      data.dirty = DIRTY_NONE;
    }
    opengl.bind_mesh();

    // Initialize uniform
    glViewport(0, 0, width, height);

    model = Eigen::Matrix4f::Identity();
    view  = Eigen::Matrix4f::Identity();
    proj  = Eigen::Matrix4f::Identity();

    // Set view
    view = lookAt(Eigen::Vector3f(options.camera_eye[0], options.camera_eye[1], options.camera_eye[2]),
                  Eigen::Vector3f(options.camera_center[0], options.camera_center[1], options.camera_center[2]),
                  Eigen::Vector3f(options.camera_up[0], options.camera_up[1], options.camera_up[2]));

    // Set projection
    if (options.orthographic)
    {
      float length = (options.camera_eye - options.camera_center).norm();
      float h = tan(options.camera_view_angle/360.0 * M_PI) * (length);
      proj = ortho(-h*width/height, h*width/height, -h, h, options.camera_dnear, options.camera_dfar);
    }
    else
    {
      float fH = tan(options.camera_view_angle / 360.0 * M_PI) * options.camera_dnear;
      float fW = fH * (double)width/(double)height;
      proj = frustum(-fW, fW, -fH, fH, options.camera_dnear, options.camera_dfar);
    }
    // end projection

    // Set model transformation
    float mat[16];
    igl::quat_to_mat(options.trackball_angle.data(), mat);

    for (unsigned i=0;i<4;++i)
      for (unsigned j=0;j<4;++j)
        model(i,j) = mat[i+4*j];

    model = scale(model, Eigen::Vector3f(options.camera_zoom,options.camera_zoom,options.camera_zoom));
    model = scale(model, Eigen::Vector3f(options.model_zoom,options.model_zoom,options.model_zoom));
    model = translate(model, Eigen::Vector3f(options.model_translation[0],options.model_translation[1],options.model_translation[2]));

    // Send transformations to the GPU
    GLint modeli = opengl.shader_mesh.uniform("model");
    GLint viewi  = opengl.shader_mesh.uniform("view");
    GLint proji  = opengl.shader_mesh.uniform("proj");
    glUniformMatrix4fv(modeli, 1, GL_FALSE, model.data());
    glUniformMatrix4fv(viewi, 1, GL_FALSE, view.data());
    glUniformMatrix4fv(proji, 1, GL_FALSE, proj.data());

    // Light parameters
    GLint specular_exponenti    = opengl.shader_mesh.uniform("specular_exponent");
    GLint light_position_worldi = opengl.shader_mesh.uniform("light_position_world");
    GLint lighting_factori      = opengl.shader_mesh.uniform("lighting_factor");
    GLint fixed_colori          = opengl.shader_mesh.uniform("fixed_color");
    GLint texture_factori       = opengl.shader_mesh.uniform("texture_factor");

    glUniform1f(specular_exponenti, options.shininess);
    Vector3f rev_light = -1.*options.light_position;
    glUniform3fv(light_position_worldi, 1, rev_light.data());
    glUniform1f(lighting_factori, 1.0f); // enables lighting
    glUniform4f(fixed_colori, 0.0, 0.0, 0.0, 0.0);

    if (data.V.rows()>0)
    {
      // Render fill
      if (options.show_faces)
      {
        // Texture
        glUniform1f(texture_factori, options.show_texture ? 1.0f : 0.0f);
        opengl.draw_mesh(true);
        glUniform1f(texture_factori, 0.0f);
      }

      // Render wireframe
      if (options.show_lines)
      {
        glLineWidth(options.line_width);
        glUniform4f(fixed_colori, options.line_color[0], options.line_color[1],
          options.line_color[2], 1.0f);
        opengl.draw_mesh(false);
        glUniform4f(fixed_colori, 0.0f, 0.0f, 0.0f, 0.0f);
      }

      if (options.show_vertid)
      {
        __font_renderer.BeginDraw(view*model, proj, viewport, data.object_scale);
        for (int i=0; i<data.V.rows(); ++i)
          __font_renderer.DrawText(data.V.row(i), data.V_normals.row(i), to_string(i));
        __font_renderer.EndDraw();
      }

      if (options.show_faceid)
      {
        __font_renderer.BeginDraw(view*model, proj, viewport, data.object_scale);

        for (int i=0; i<data.F.rows(); ++i)
        {
          Eigen::RowVector3d p = Eigen::RowVector3d::Zero();
          for (int j=0;j<data.F.cols();++j)
            p += data.V.row(data.F(i,j));
          p /= data.F.cols();

          __font_renderer.DrawText(p, data.F_normals.row(i), to_string(i));
        }
        __font_renderer.EndDraw();
      }
    }

    if (options.show_overlay)
    {
      if (options.show_overlay_depth)
        glEnable(GL_DEPTH_TEST);
      else
        glDisable(GL_DEPTH_TEST);

      if (data.lines.rows() > 0)
      {
        opengl.bind_overlay_lines();
        modeli = opengl.shader_overlay_lines.uniform("model");
        viewi  = opengl.shader_overlay_lines.uniform("view");
        proji  = opengl.shader_overlay_lines.uniform("proj");

        glUniformMatrix4fv(modeli, 1, GL_FALSE, model.data());
        glUniformMatrix4fv(viewi, 1, GL_FALSE, view.data());
        glUniformMatrix4fv(proji, 1, GL_FALSE, proj.data());
        glLineWidth(options.line_width);

        opengl.draw_overlay_lines();
      }

      if (data.points.rows() > 0)
      {
        opengl.bind_overlay_points();
        modeli = opengl.shader_overlay_points.uniform("model");
        viewi  = opengl.shader_overlay_points.uniform("view");
        proji  = opengl.shader_overlay_points.uniform("proj");

        glUniformMatrix4fv(modeli, 1, GL_FALSE, model.data());
        glUniformMatrix4fv(viewi, 1, GL_FALSE, view.data());
        glUniformMatrix4fv(proji, 1, GL_FALSE, proj.data());
        glPointSize(options.point_size);

        opengl.draw_overlay_points();
      }

      if (data.labels_positions.rows() > 0)
      {
        __font_renderer.BeginDraw(view*model, proj, viewport, data.object_scale);
        for (int i=0; i<data.labels_positions.rows(); ++i)
          __font_renderer.DrawText(data.labels_positions.row(i), Eigen::Vector3d(0.0,0.0,0.0),
              data.labels_strings[i]);
        __font_renderer.EndDraw();
      }

      glEnable(GL_DEPTH_TEST);
    }

    if (callback_post_draw)
      if (callback_post_draw(*this))
        return;

    if (plugin_manager)
      for (unsigned int i = 0; i <plugin_manager->plugin_list.size(); ++i)
        if (plugin_manager->plugin_list[i]->post_draw())
          break;

    TwDraw();
  }

  bool Viewer::save_scene()
  {
    #ifdef ENABLE_XML_SERIALIZATION
    string fname = igl::file_dialog_save();
    if (fname.length() == 0)
      return false;

    ::igl::XMLSerializer serializer("Viewer");
    serializer.Add(data,"Data");
    serializer.Add(options,"Options");

    if (plugin_manager)
      for (unsigned int i = 0; i <plugin_manager->plugin_list.size(); ++i)
        serializer.Add(*(plugin_manager->plugin_list[i]),plugin_manager->plugin_list[i]->plugin_name);

    serializer.Save(fname.c_str(),true);

    #endif
    return true;
  }

  bool Viewer::load_scene()
  {
    #ifdef ENABLE_XML_SERIALIZATION
    string fname = igl::file_dialog_open();
    if (fname.length() == 0)
      return false;

    ::igl::XMLSerializer serializer("Viewer");
    serializer.Add(data,"Data");
    serializer.Add(options,"Options");

    if (plugin_manager)
      for (unsigned int i = 0; i <plugin_manager->plugin_list.size(); ++i)
        serializer.Add(*(plugin_manager->plugin_list[i]),plugin_manager->plugin_list[i]->plugin_name);

    serializer.Load(fname.c_str());

    #endif
    return true;
  }

  void Viewer::align_camera_center()
  {
    get_scale_and_shift_to_fit_mesh(data.V,data.F,options.model_zoom,options.model_translation);
    data.object_scale = (data.V.colwise().maxCoeff() - data.V.colwise().minCoeff()).norm();
  }

  void Viewer::resize(int w, int h)
  {
    width = w;
    height = h;
    viewport = Eigen::Vector4f(0,0,width,height);
  }

  void Viewer::get_scale_and_shift_to_fit_mesh(
    const Eigen::MatrixXd& V,
    const Eigen::MatrixXi& F,
    float& zoom,
    Eigen::Vector3f& shift)
  {
    if (V.rows() == 0)
      return;
    //Compute mesh centroid
    Eigen::SparseMatrix<double> M;
    igl::massmatrix(V,F,igl::MASSMATRIX_VORONOI,M);
    const auto & MV = M*V;
    Eigen::RowVector3d centroid  = MV.colwise().sum()/M.diagonal().sum();
    Eigen::RowVector3d min_point = V.colwise().minCoeff();
    Eigen::RowVector3d max_point = V.colwise().maxCoeff();

    shift = -centroid.cast<float>();
    double x_scale = fabs(max_point[0] - min_point[0]);
    double y_scale = fabs(max_point[1] - min_point[1]);
    double z_scale = fabs(max_point[2] - min_point[2]);
    zoom = 2.0/ std::max(z_scale,std::max(x_scale,y_scale));
  }

  void TW_CALL Viewer::snap_to_canonical_quaternion_cb(void *clientData)
  {
    Eigen::Vector4f snapq = static_cast<Viewer *>(clientData)->options.trackball_angle;
    igl::snap_to_canonical_view_quat<float>(snapq.data(),1,static_cast<Viewer *>(clientData)->options.trackball_angle.data());
  }
  void TW_CALL Viewer::align_camera_center_cb(void *clientData)
  {
    static_cast<Viewer *>(clientData)->align_camera_center();
  }

  void TW_CALL Viewer::save_scene_cb(void *clientData)
  {
    static_cast<Viewer *>(clientData)->save_scene();
  }

  void TW_CALL Viewer::load_scene_cb(void *clientData)
  {
    static_cast<Viewer *>(clientData)->load_scene();
  }

  void TW_CALL Viewer::set_invert_normals_cb(const void *param, void *clientData)
  {
    Viewer *viewer = static_cast<Viewer *>(clientData);
    viewer->data.dirty |= Viewer::DIRTY_NORMAL;
    viewer->options.invert_normals = *((bool *) param);
  }

  void TW_CALL Viewer::get_invert_normals_cb(void *param, void *clientData)
  {
    *((bool *) param) = static_cast<Viewer *>(clientData)->options.invert_normals;
  }

  void TW_CALL Viewer::set_face_based_cb(const void *param, void *clientData)
  {
    Viewer *viewer = static_cast<Viewer *>(clientData);
    viewer->set_face_based(*((bool *) param));
  }

  void TW_CALL Viewer::get_face_based_cb(void *param, void *clientData)
  {
    *((bool *) param) = static_cast<Viewer *>(clientData)->options.face_based;
  }

  void TW_CALL Viewer::open_dialog_mesh(void *clientData)
  {
    std::string fname = igl::file_dialog_open();

    if (fname.length() == 0)
      return;

    static_cast<Viewer *>(clientData)->load_mesh_from_file(fname.c_str());
  }

  // Serialization
  void Viewer::Options::InitSerialization()
  {
    #ifdef ENABLE_XML_SERIALIZATION
    xmlSerializer->Add(shininess, "shininess");
    xmlSerializer->Add(background_color, "background_color");
    xmlSerializer->Add(line_color, "line_color");
    xmlSerializer->Add(light_position, "light_position");
    xmlSerializer->Add(trackball_angle, "trackball_angle");
    xmlSerializer->Add(model_zoom, "model_zoom");
    xmlSerializer->Add(model_translation, "model_translation");
    xmlSerializer->Add(model_zoom_uv, "model_zoom_uv");
    xmlSerializer->Add(model_translation_uv, "model_translation_uv");
    xmlSerializer->Add(camera_zoom, "camera_zoom");
    xmlSerializer->Add(orthographic, "orthographic");
    xmlSerializer->Add(camera_eye, "camera_eye");
    xmlSerializer->Add(camera_up, "camera_up");
    xmlSerializer->Add(camera_center, "camera_center");
    xmlSerializer->Add(camera_view_angle, "camera_view_angle");
    xmlSerializer->Add(camera_dnear, "camera_dnear");
    xmlSerializer->Add(camera_dfar, "camera_dfar");
    xmlSerializer->Add(show_overlay, "show_overlay");
    xmlSerializer->Add(show_overlay_depth, "show_overlay_depth");
    xmlSerializer->Add(show_texture, "show_texture");
    xmlSerializer->Add(show_faces, "show_faces");
    xmlSerializer->Add(show_lines, "show_lines");
    xmlSerializer->Add(show_vertid, "show_vertid");
    xmlSerializer->Add(show_faceid, "show_faceid");
    xmlSerializer->Add(point_size, "point_size");
    xmlSerializer->Add(line_width, "line_width");
    xmlSerializer->Add(invert_normals, "invert_normals");
    xmlSerializer->Add(face_based, "face_based");
    #endif
  }

  void Viewer::Data::InitSerialization()
  {
    #ifdef ENABLE_XML_SERIALIZATION
    xmlSerializer->Add(V,"V");
    xmlSerializer->Add(F,"F");
    xmlSerializer->Add(F_normals,"F_normals");

    xmlSerializer->Add(F_material_ambient,"F_material_ambient");
    xmlSerializer->Add(F_material_diffuse,"F_material_diffuse");
    xmlSerializer->Add(F_material_specular,"F_material_specular");

    xmlSerializer->Add(V_normals,"V_normals");
    xmlSerializer->Add(V_material_ambient,"V_material_ambient");
    xmlSerializer->Add(V_material_diffuse,"V_material_diffuse");
    xmlSerializer->Add(V_material_specular,"V_material_specular");

    xmlSerializer->Add(V_uv,"V_uv");
    xmlSerializer->Add(F_uv,"F_uv");
    xmlSerializer->Add(texture_R,"texture_R");
    xmlSerializer->Add(texture_G,"texture_G");
    xmlSerializer->Add(texture_B,"texture_B");
    xmlSerializer->Add(lines,"lines");
    xmlSerializer->Add(points,"points");

    xmlSerializer->Add(labels_positions,"labels_positions");
    xmlSerializer->Add(labels_strings,"labels_strings");
    #endif
  }

  void Viewer::set_face_based(bool newvalue)
  {
    if (options.face_based != newvalue)
    {
      options.face_based = newvalue;
      data.dirty = DIRTY_ALL;
    }
  }

  // Helpers that draws the most common meshes
  void Viewer::set_mesh(const Eigen::MatrixXd& V, const Eigen::MatrixXi& F)
  {
    using namespace std;

    Eigen::MatrixXd V_temp;

    // If V only has two columns, pad with a column of zeros
    if (V.cols() == 2)
    {
      V_temp = Eigen::MatrixXd::Zero(V.rows(),3);
      V_temp.block(0,0,V.rows(),2) = V;
    }
    else
      V_temp = V;

    if (data.V.rows() == 0 && data.F.rows() == 0)
    {
      clear_mesh();
      data.V = V_temp;
      data.F = F;

      compute_normals();
      uniform_colors(Eigen::Vector3d(51.0/255.0,43.0/255.0,33.3/255.0),
                     Eigen::Vector3d(255.0/255.0,228.0/255.0,58.0/255.0),
                     Eigen::Vector3d(255.0/255.0,235.0/255.0,80.0/255.0));

      grid_texture();
      align_camera_center();
    }
    else
    {
      if (data.V.rows() == V.rows() && data.F.rows() == F.rows())
      {
        data.V = V_temp;
        data.F = F;
        align_camera_center();
      }
      else
        cerr << "ERROR (set_mesh): The new mesh has a different number of vertices/faces. Please clear the mesh before plotting.";
    }
    data.dirty |= DIRTY_FACE | DIRTY_POSITION;
  }

  void Viewer::set_vertices(const Eigen::MatrixXd& V)
  {
    data.V = V;
    assert(data.F.size() == 0 || data.F.maxCoeff() < data.V.rows());
    data.dirty |= DIRTY_POSITION;
  }

  void Viewer::set_normals(const Eigen::MatrixXd& N)
  {
    using namespace std;
    if (N.rows() == data.V.rows())
    {
      set_face_based(false);
      data.V_normals = N;
    }
    else if (N.rows() == data.F.rows() || N.rows() == data.F.rows()*3)
    {
      set_face_based(true);
      data.F_normals = N;
    }
    else
      cerr << "ERROR (set_normals): Please provide a normal per face, per corner or per vertex.";
    data.dirty |= DIRTY_NORMAL;
  }

  void Viewer::set_colors(const Eigen::MatrixXd &C)
  {
    using namespace std;
    using namespace Eigen;
    // Ambient color should be darker color
    const auto ambient = [](const MatrixXd & C)->MatrixXd
    {
      return 0.1*C;
    };
    // Specular color should be a less saturated and darker color: dampened
    // highlights
    const auto specular = [](const MatrixXd & C)->MatrixXd
    {
      const double grey = 0.3;
      return grey+0.1*(C.array()-grey);
    };
    if (C.rows() == 1)
    {
      for (unsigned i=0;i<data.V_material_diffuse.rows();++i)
      {
        data.V_material_diffuse.row(i) = C.row(0);
      }
      data.V_material_ambient = ambient(data.V_material_diffuse);
      data.V_material_specular = specular(data.V_material_diffuse);

      for (unsigned i=0;i<data.F_material_diffuse.rows();++i)
      {
        data.F_material_diffuse.row(i) = C.row(0);
      }
      data.F_material_ambient = ambient(data.F_material_diffuse);
      data.F_material_specular = specular(data.F_material_diffuse);
    }
    else if (C.rows() == data.V.rows())
    {
      set_face_based(false);
      data.V_material_diffuse = C;
      data.V_material_ambient = ambient(data.V_material_diffuse);
      data.V_material_specular = specular(data.V_material_diffuse);
    }
    else if (C.rows() == data.F.rows())
    {
      set_face_based(true);
      data.F_material_diffuse = C;
      data.F_material_ambient = ambient(data.F_material_diffuse);
      data.F_material_specular = specular(data.F_material_diffuse);
    }
    else
      cerr << "ERROR (set_colors): Please provide a single color, or a color per face or per vertex.";
    data.dirty |= DIRTY_DIFFUSE;

  }

  void Viewer::set_uv(const Eigen::MatrixXd& UV)
  {
    using namespace std;
    if (UV.rows() == data.V.rows())
    {
      set_face_based(false);
      data.V_uv = UV;
    }
    else
      cerr << "ERROR (set_UV): Please provide uv per vertex.";
    data.dirty |= DIRTY_UV;
  }

  void Viewer::set_uv(const Eigen::MatrixXd& UV_V, const Eigen::MatrixXi& UV_F)
  {
    set_face_based(true);
    data.V_uv = UV_V;
    data.F_uv = UV_F;
    data.dirty |= DIRTY_UV;
  }


  void Viewer::set_texture(
    const Eigen::Matrix<char,Eigen::Dynamic,Eigen::Dynamic>& R,
    const Eigen::Matrix<char,Eigen::Dynamic,Eigen::Dynamic>& G,
    const Eigen::Matrix<char,Eigen::Dynamic,Eigen::Dynamic>& B)
  {
    data.texture_R = R;
    data.texture_G = G;
    data.texture_B = B;
    data.dirty |= DIRTY_TEXTURE;
  }

  void Viewer::add_points(const Eigen::MatrixXd& P,  const Eigen::MatrixXd& C)
  {
    Eigen::MatrixXd P_temp;
    
    // If P only has two columns, pad with a column of zeros
    if (P.cols() == 2)
    {
      P_temp = Eigen::MatrixXd::Zero(P.rows(),3);
      P_temp.block(0,0,P.rows(),2) = P;
    }
    else
      P_temp = P;

    int lastid = data.points.rows();
    data.points.conservativeResize(data.points.rows() + P_temp.rows(),6);
    for (unsigned i=0; i<P_temp.rows(); ++i)
      data.points.row(lastid+i) << P_temp.row(i), i<C.rows() ? C.row(i) : C.row(C.rows()-1);

    data.dirty |= DIRTY_OVERLAY_POINTS;
  }

  void Viewer::add_edges(const Eigen::MatrixXd& P1, const Eigen::MatrixXd& P2, const Eigen::MatrixXd& C)
  {
    Eigen::MatrixXd P1_temp,P2_temp;
    
    // If P1 only has two columns, pad with a column of zeros
    if (P1.cols() == 2)
    {
      P1_temp = Eigen::MatrixXd::Zero(P1.rows(),3);
      P1_temp.block(0,0,P1.rows(),2) = P1;
      P2_temp = Eigen::MatrixXd::Zero(P2.rows(),3);
      P2_temp.block(0,0,P2.rows(),2) = P2;
    }
    else
    {
      P1_temp = P1;
      P2_temp = P2;
    }

    int lastid = data.lines.rows();
    data.lines.conservativeResize(data.lines.rows() + P1_temp.rows(),9);
    for (unsigned i=0; i<P1_temp.rows(); ++i)
      data.lines.row(lastid+i) << P1_temp.row(i), P2_temp.row(i), i<C.rows() ? C.row(i) : C.row(C.rows()-1);

    data.dirty |= DIRTY_OVERLAY_LINES;
  }

  void Viewer::add_label(const Eigen::VectorXd& P,  const std::string& str)
  {
    Eigen::MatrixXd P_temp;
    
    // If P only has two columns, pad with a column of zeros
    if (P.cols() == 2)
    {
      P_temp = Eigen::MatrixXd::Zero(P.rows(),3);
      P_temp.block(0,0,P.rows(),2) = P;
    }
    else
      P_temp = P;

    int lastid = data.labels_positions.rows();
    data.labels_positions.conservativeResize(lastid+1, 3);
    data.labels_positions.row(lastid) = P_temp;
    data.labels_strings.push_back(str);
  }

  int Viewer::launch(std::string filename)
  {
    GLFWwindow* window;

    glfwSetErrorCallback(glfw_error_callback);
    if (!glfwInit())
      return EXIT_FAILURE;

    glfwWindowHint(GLFW_SAMPLES, 16);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MAJOR, 3);
    glfwWindowHint(GLFW_CONTEXT_VERSION_MINOR, 2);
    glfwWindowHint(GLFW_OPENGL_PROFILE, GLFW_OPENGL_CORE_PROFILE);
    glfwWindowHint(GLFW_OPENGL_FORWARD_COMPAT, GL_TRUE);

    window = glfwCreateWindow(1280, 800, "IGL Viewer", NULL, NULL);
    if (!window)
    {
      glfwTerminate();
      return EXIT_FAILURE;
    }

	glfwMakeContextCurrent(window);

#ifdef _WIN32
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

    // Initialize IGL viewer
    init(&igl_viewer_plugin_manager);
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

    init_opengl();

    // Load the mesh passed as input
    if (filename.size() > 0)
      load_mesh_from_file(filename.c_str());

    // Rendering loop
    while (!glfwWindowShouldClose(window))
    {
      draw();

      glfwSwapBuffers(window);
      //glfwPollEvents();
      glfwWaitEvents();
    }

    free_opengl();
    shutdown_plugins();

    glfwDestroyWindow(window);
    glfwTerminate();
    return EXIT_SUCCESS;
  }
} // end namespace
