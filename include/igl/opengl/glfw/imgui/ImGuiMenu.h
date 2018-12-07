// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Jérémie Dumas <jeremie.dumas@ens-lyon.org>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef IGL_OPENGL_GLFW_IMGUI_IMGUIMENU_H
#define IGL_OPENGL_GLFW_IMGUI_IMGUIMENU_H

////////////////////////////////////////////////////////////////////////////////
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/ViewerPlugin.h>
#include <igl/opengl/gl.h>
#include <igl/igl_inline.h>
#include <memory>
#include <exception>
////////////////////////////////////////////////////////////////////////////////

// Forward declaration
struct ImGuiContext;
struct ImGui_ImplGlfw_Data;

namespace igl
{
namespace opengl
{
namespace glfw
{
namespace imgui
{

class ImGuiMenu : public igl::opengl::glfw::ViewerPlugin
{
  // RAII for ImGuiContext
  struct ScopedContext
  {
    ImGuiContext *old_context_;
    ImGui_ImplGlfw_Data *old_binding_data_;

    ScopedContext(ImGuiContext *ctx, ImGui_ImplGlfw_Data *binding_data);
    ~ScopedContext();

    ScopedContext(ScopedContext &&) = delete;
    ScopedContext &operator=(ScopedContext &&) = delete;
    ScopedContext(const ScopedContext &) = delete;
    ScopedContext &operator=(ScopedContext &) = delete;
  };

public:
  ImGuiMenu() = default;
  virtual ~ImGuiMenu() = default;

private:
  // Delete copy and move constructors, since we need to manage the internal ImGuiContext *
  ImGuiMenu(ImGuiMenu&&) = delete;
  ImGuiMenu& operator=(ImGuiMenu&&) = delete;
  ImGuiMenu(const ImGuiMenu&) = delete;
  ImGuiMenu& operator=(const ImGuiMenu&) = delete;

protected:
  // Hidpi scaling to be used for text rendering.
  float hidpi_scaling_;

  // Ratio between the framebuffer size and the window size.
  // May be different from the hipdi scaling!
  float pixel_ratio_;

  // ImGui context
  ImGuiContext * context_ = nullptr;
  ImGuiContext * old_context_ = nullptr;

  // ImGui binding data
  ImGui_ImplGlfw_Data * binding_data_ = nullptr;
  ImGui_ImplGlfw_Data * old_binding_data_ = nullptr;

public:
  IGL_INLINE virtual void init(igl::opengl::glfw::Viewer *_viewer) override;

  IGL_INLINE void init_imgui();

  IGL_INLINE virtual void reload_font(int font_size = 13);

  IGL_INLINE virtual void shutdown() override;

  IGL_INLINE virtual bool pre_draw() override;

  IGL_INLINE  virtual bool post_draw() override;

  IGL_INLINE virtual void post_resize(int width, int height) override;

  // Mouse IO
  IGL_INLINE virtual bool mouse_down(int button, int modifier) override;

  IGL_INLINE virtual bool mouse_up(int button, int modifier) override;

  IGL_INLINE virtual bool mouse_move(int mouse_x, int mouse_y) override;

  IGL_INLINE virtual bool mouse_scroll(float delta_y) override;

  // Keyboard IO
  IGL_INLINE virtual bool key_pressed(unsigned int key, int modifiers) override;

  IGL_INLINE virtual bool key_down(int key, int modifiers) override;

  IGL_INLINE virtual bool key_up(int key, int modifiers) override;

  // Draw menu
  IGL_INLINE virtual void draw_menu();

  // Can be overwritten by `callback_draw_viewer_window`
  IGL_INLINE virtual void draw_viewer_window();

  // Can be overwritten by `callback_draw_viewer_menu`
  IGL_INLINE virtual void draw_viewer_menu();

  // Can be overwritten by `callback_draw_custom_window`
  IGL_INLINE virtual void draw_custom_window() { }

  // Easy-to-customize callbacks
  std::function<void(void)> callback_draw_viewer_window;
  std::function<void(void)> callback_draw_viewer_menu;
  std::function<void(void)> callback_draw_custom_window;

  IGL_INLINE void draw_labels_window();

  IGL_INLINE void draw_labels(const igl::opengl::ViewerData &data);

  IGL_INLINE void draw_text(Eigen::Vector3d pos, Eigen::Vector3d normal, const std::string &text);

  IGL_INLINE float pixel_ratio();

  IGL_INLINE float hidpi_scaling();

  float menu_scaling() { return hidpi_scaling_ / pixel_ratio_; }
};

} // end namespace
} // end namespace
} // end namespace
} // end namespace

#ifndef IGL_STATIC_LIBRARY
#  include "ImGuiMenu.cpp"
#endif

#endif // IGL_OPENGL_GLFW_IMGUI_IMGUIMENU_H
