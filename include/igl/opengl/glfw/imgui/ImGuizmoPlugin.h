#ifndef IGL_OPENGL_GFLW_IMGUI_IMGUIZMOPLUGIN_H
#define IGL_OPENGL_GFLW_IMGUI_IMGUIZMOPLUGIN_H
#include "../../../igl_inline.h"
#include "ImGuiMenu.h"
#include <imgui/imgui.h>
#include <imgui/imgui_internal.h>
#include <imguizmo/ImGuizmo.h>
#include <Eigen/Dense>

namespace igl{ namespace opengl{ namespace glfw{ namespace imgui{

class ImGuizmoPlugin : public igl::opengl::glfw::imgui::ImGuiMenu
{
public:
  // callback(T) called when the stored transform T changes
  std::function<void(const Eigen::Matrix4f &)> callback;
  std::function<bool(ImGuizmoPlugin & plugin, int button, int modifier)> callback_mouse_down = nullptr;
  std::function<bool(ImGuizmoPlugin & plugin, int button, int modifier)> callback_mouse_up = nullptr;
  std::function<bool(ImGuizmoPlugin & plugin, int mouse_x,int mouse_y )> callback_mouse_move = nullptr;
  // Whether to display
  bool visible = true;
  // whether rotating, translating or scaling
  ImGuizmo::OPERATION operation;
  // stored transformation
  Eigen::Matrix4f T;
  // Initilize with rotate operation on an identity transform (at origin)
  ImGuizmoPlugin():operation(ImGuizmo::ROTATE),T(Eigen::Matrix4f::Identity()){};
  /////////////////////////////////////////////////////////////////////////////
  // Boilerplate
  IGL_INLINE virtual void init(igl::opengl::glfw::Viewer *_viewer) override;
  IGL_INLINE virtual bool pre_draw() override;
  IGL_INLINE virtual bool post_draw() override;
  IGL_INLINE virtual bool mouse_down(int button, int modifier) override;
  IGL_INLINE virtual bool mouse_up(int button, int modifier) override;
  IGL_INLINE virtual bool mouse_move(int mouse_x, int mouse_y) override;
  IGL_INLINE virtual bool mouse_scroll(float delta_y) override;
};

}}}}

#ifndef IGL_STATIC_LIBRARY
#  include "ImGuizmoPlugin.cpp"
#endif

#endif
