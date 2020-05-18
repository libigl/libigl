// This file is part of libigl, a simple c++ geometry processing library.
//
// Copyright (C) 2018 Jérémie Dumas <jeremie.dumas@ens-lyon.org>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
////////////////////////////////////////////////////////////////////////////////
#include "ImGuiMenu.h"
#include "ImGuiHelpers.h"
#include <igl/project.h>
#include <igl/unproject_ray.h>
#include <imgui/imgui.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <imgui_fonts_droid_sans.h>
#include <GLFW/glfw3.h>
#include <iostream>
////////////////////////////////////////////////////////////////////////////////

namespace igl
{
namespace opengl
{
namespace glfw
{
namespace imgui
{

IGL_INLINE void ImGuiMenu::init(igl::opengl::glfw::Viewer *_viewer)
{
  ViewerPlugin::init(_viewer);
  // Setup ImGui binding
  if (_viewer)
  {
    IMGUI_CHECKVERSION();
    if (!context_)
    {
      // Single global context by default, but can be overridden by the user
      static ImGuiContext * __global_context = ImGui::CreateContext();
      context_ = __global_context;
    }
    const char* glsl_version = "#version 150";
    ImGui_ImplGlfw_InitForOpenGL(viewer->window, false);
    ImGui_ImplOpenGL3_Init(glsl_version);
    ImGui::GetIO().IniFilename = nullptr;
    ImGui::StyleColorsDark();
    ImGuiStyle& style = ImGui::GetStyle();
    style.FrameRounding = 5.0f;
    reload_font();
  }
}

IGL_INLINE void ImGuiMenu::reload_font(int font_size)
{
  hidpi_scaling_ = hidpi_scaling();
  pixel_ratio_ = pixel_ratio();
  ImGuiIO& io = ImGui::GetIO();
  io.Fonts->Clear();
  io.Fonts->AddFontFromMemoryCompressedTTF(droid_sans_compressed_data,
    droid_sans_compressed_size, font_size * hidpi_scaling_);
  io.FontGlobalScale = 1.0 / pixel_ratio_;
}

IGL_INLINE void ImGuiMenu::shutdown()
{
  // Cleanup
  ImGui_ImplOpenGL3_Shutdown();
  ImGui_ImplGlfw_Shutdown();
  // User is responsible for destroying context if a custom context is given
  // ImGui::DestroyContext(*context_);
}

IGL_INLINE bool ImGuiMenu::pre_draw()
{
  glfwPollEvents();

  // Check whether window dpi has changed
  float scaling = hidpi_scaling();
  if (std::abs(scaling - hidpi_scaling_) > 1e-5)
  {
    reload_font();
    ImGui_ImplOpenGL3_DestroyDeviceObjects();
  }

  ImGui_ImplOpenGL3_NewFrame();
  ImGui_ImplGlfw_NewFrame();
  ImGui::NewFrame();
  return false;
}

IGL_INLINE bool ImGuiMenu::post_draw()
{
  draw_menu();
  ImGui::Render();
  ImGui_ImplOpenGL3_RenderDrawData(ImGui::GetDrawData());
  return false;
}

IGL_INLINE void ImGuiMenu::post_resize(int width, int height)
{
  if (context_)
  {
    ImGui::GetIO().DisplaySize.x = float(width);
    ImGui::GetIO().DisplaySize.y = float(height);
  }
}

// Mouse IO
IGL_INLINE bool ImGuiMenu::mouse_down(int button, int modifier)
{
  ImGui_ImplGlfw_MouseButtonCallback(viewer->window, button, GLFW_PRESS, modifier);
  return ImGui::GetIO().WantCaptureMouse;
}

IGL_INLINE void ImGuiMenu::filterLabelsByDepth()
{

  float closestDepth = -100.;
  float colinearityThreshold = -0.1;
  float depthThreshold = 1;

  int width = viewer->core().viewport(2);
  int height = viewer->core().viewport(3);
  bool faceVisible;
  std::vector<int> candidateFaceIds;
  viewer->data().vertex_label_mask = Eigen::MatrixXi::Zero(viewer->data().V.rows(), 1);
  viewer->data().face_label_mask = Eigen::MatrixXi::Zero(viewer->data().F.rows(), 1);

  Eigen::MatrixXi viableVertIds = Eigen::MatrixXi::Zero(viewer->data().V.rows(), 1);

  // Get unprojected ray shooting into the middle of the screen
  Eigen::Vector3f s,dir;
  igl::unproject_ray(
    Eigen::Vector2f(width/2, height/2),
    viewer->core().view, 
    viewer->core().proj, 
    viewer->core().viewport,
    s,
    dir);

  // We use unprojected face normals to detemrmine which are pointing towards the ray
  Eigen::MatrixXf faceNormals = (viewer->data().F_normals.cast <float> ()).transpose();
  // Original transformation matrix for projecting the face normals
  // Eigen::MatrixXf normalMatrix = (((viewer->core().view).block<3,3>(0,0)).inverse()).transpose();

  // See which face ids are in the viewport AND are facing the camera
  for(int f=0; f<viewer->data().F.rows(); f++)
  {
    float dotProduct = ((dir).normalized()).dot(faceNormals.col(f).normalized());
    // Ray and face normal are coliniear and opposite direction
    if(dotProduct <= colinearityThreshold)
    {
      faceVisible = true; // Reset boolean
      // Get screen coordinates of each vertex of current face
      Eigen::Vector3f currScreenCoords;
      for(int i=0; i<=2; i++)
      {
        int vertIdx = viewer->data().F.row(f)(i);
        if(viableVertIds(vertIdx, 0)) continue;

        // Eigen::Vector3d p = viewer->data().V.row(vertIdx) + viewer->data().V_normals.row(vertIdx) * 0.005f * viewer->core().object_scale;
        Eigen::Vector3d p = viewer->data().V.row(vertIdx);
        currScreenCoords = igl::project(
          Eigen::Vector3f(p.cast<float>()), 
          viewer->core().view, 
          viewer->core().proj, 
          viewer->core().viewport
        );
        float x = currScreenCoords(0);
        float y = currScreenCoords(1);
        if(x < 0 || x >= width || y < 0 || y >= height) 
        {
          faceVisible = false;
          break;
        }

        // Obtain depth of vertex relative to the camera
        Eigen::Vector4f currentVertex(
          viewer->data().V.row(vertIdx)(0), 
          viewer->data().V.row(vertIdx)(1), 
          viewer->data().V.row(vertIdx)(2), 
          1.0);
        Eigen::MatrixXf viewSpaceVertex = viewer->core().view*(currentVertex);
        if(viewSpaceVertex(2) > closestDepth)
          closestDepth = viewSpaceVertex(2);
      }
      if(faceVisible)
      {
        candidateFaceIds.push_back(f);
        // Cache visited verts to optimize this nested for loop
        for(int i=0;i<=2;i++)
          viableVertIds(viewer->data().F.row(f)(i), 0) = 1;
      }
    }
  }

  // Of the candidate faces only select ones that are closest to screen.
  // This calculation is fast enough so that this function can
  // be invoked on mouse_move, not just on mouse_up.
  // Alternatively we use the unproject_onto_mesh function which 
  // is more accurate but much slower.
  for(int f=0; f < candidateFaceIds.size(); f++)
  {
    int v1 = viewer->data().F.row(candidateFaceIds[f])(0) ; 
    int v2 = viewer->data().F.row(candidateFaceIds[f])(1) ; 
    int v3 = viewer->data().F.row(candidateFaceIds[f])(2) ;

    Eigen::Vector4f currentVertex(
      viewer->data().V.row(v1)(0), 
      viewer->data().V.row(v1)(1), 
      viewer->data().V.row(v1)(2), 
      1.0);
    Eigen::MatrixXf viewSpaceVertex = viewer->core().view*( currentVertex );
    if(viewSpaceVertex(2) >= closestDepth-depthThreshold)
    {
      viewer->data().face_label_mask(candidateFaceIds[f], 0) = 1;
      viewer->data().vertex_label_mask(v1,0)=1;
      viewer->data().vertex_label_mask(v2,0)=1;
      viewer->data().vertex_label_mask(v3,0)=1;
    }
  }
}

IGL_INLINE bool ImGuiMenu::mouse_up(int button, int modifier)
{
  //return ImGui::GetIO().WantCaptureMouse;
  // !! Should not steal mouse up
  return false;
}

IGL_INLINE bool ImGuiMenu::mouse_move(int mouse_x, int mouse_y)
{
  return ImGui::GetIO().WantCaptureMouse;
}

IGL_INLINE bool ImGuiMenu::mouse_scroll(float delta_y)
{
  ImGui_ImplGlfw_ScrollCallback(viewer->window, 0.f, delta_y);
  return ImGui::GetIO().WantCaptureMouse;
}

// Keyboard IO
IGL_INLINE bool ImGuiMenu::key_pressed(unsigned int key, int modifiers)
{
  ImGui_ImplGlfw_CharCallback(nullptr, key);
  return ImGui::GetIO().WantCaptureKeyboard;
}

IGL_INLINE bool ImGuiMenu::key_down(int key, int modifiers)
{
  ImGui_ImplGlfw_KeyCallback(viewer->window, key, 0, GLFW_PRESS, modifiers);
  return ImGui::GetIO().WantCaptureKeyboard;
}

IGL_INLINE bool ImGuiMenu::key_up(int key, int modifiers)
{
  ImGui_ImplGlfw_KeyCallback(viewer->window, key, 0, GLFW_RELEASE, modifiers);
  return ImGui::GetIO().WantCaptureKeyboard;
}

// Draw menu
IGL_INLINE void ImGuiMenu::draw_menu()
{
  // Text labels
  draw_labels_window();

  // Viewer settings
  if (callback_draw_viewer_window) { callback_draw_viewer_window(); }
  else { draw_viewer_window(); }

  // Other windows
  if (callback_draw_custom_window) { callback_draw_custom_window(); }
  else { draw_custom_window(); }
}

IGL_INLINE void ImGuiMenu::draw_viewer_window()
{
  float menu_width = 180.f * menu_scaling();
  ImGui::SetNextWindowPos(ImVec2(0.0f, 0.0f), ImGuiCond_FirstUseEver);
  ImGui::SetNextWindowSize(ImVec2(0.0f, 0.0f), ImGuiCond_FirstUseEver);
  ImGui::SetNextWindowSizeConstraints(ImVec2(menu_width, -1.0f), ImVec2(menu_width, -1.0f));
  bool _viewer_menu_visible = true;
  ImGui::Begin(
      "Viewer", &_viewer_menu_visible,
      ImGuiWindowFlags_NoSavedSettings
      | ImGuiWindowFlags_AlwaysAutoResize
  );
  ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.4f);
  if (callback_draw_viewer_menu) { callback_draw_viewer_menu(); }
  else { draw_viewer_menu(); }
  ImGui::PopItemWidth();
  ImGui::End();
}

IGL_INLINE void ImGuiMenu::draw_viewer_menu()
{
  // Workspace
  if (ImGui::CollapsingHeader("Workspace", ImGuiTreeNodeFlags_DefaultOpen))
  {
    float w = ImGui::GetContentRegionAvailWidth();
    float p = ImGui::GetStyle().FramePadding.x;
    if (ImGui::Button("Load##Workspace", ImVec2((w-p)/2.f, 0)))
    {
      viewer->load_scene();
    }
    ImGui::SameLine(0, p);
    if (ImGui::Button("Save##Workspace", ImVec2((w-p)/2.f, 0)))
    {
      viewer->save_scene();
    }
  }

  // Mesh
  if (ImGui::CollapsingHeader("Mesh", ImGuiTreeNodeFlags_DefaultOpen))
  {
    float w = ImGui::GetContentRegionAvailWidth();
    float p = ImGui::GetStyle().FramePadding.x;
    if (ImGui::Button("Load##Mesh", ImVec2((w-p)/2.f, 0)))
    {
      viewer->open_dialog_load_mesh();
    }
    ImGui::SameLine(0, p);
    if (ImGui::Button("Save##Mesh", ImVec2((w-p)/2.f, 0)))
    {
      viewer->open_dialog_save_mesh();
    }
  }

  // Viewing options
  if (ImGui::CollapsingHeader("Viewing Options", ImGuiTreeNodeFlags_DefaultOpen))
  {
    if (ImGui::Button("Center object", ImVec2(-1, 0)))
    {
      viewer->core().align_camera_center(viewer->data().V, viewer->data().F);
    }
    if (ImGui::Button("Snap canonical view", ImVec2(-1, 0)))
    {
      viewer->snap_to_canonical_quaternion();
    }

    // Zoom
    ImGui::PushItemWidth(80 * menu_scaling());
    ImGui::DragFloat("Zoom", &(viewer->core().camera_zoom), 0.05f, 0.1f, 20.0f);

    // Select rotation type
    int rotation_type = static_cast<int>(viewer->core().rotation_type);
    static Eigen::Quaternionf trackball_angle = Eigen::Quaternionf::Identity();
    static bool orthographic = true;
    if (ImGui::Combo("Camera Type", &rotation_type, "Trackball\0Two Axes\0002D Mode\0\0"))
    {
      using RT = igl::opengl::ViewerCore::RotationType;
      auto new_type = static_cast<RT>(rotation_type);
      if (new_type != viewer->core().rotation_type)
      {
        if (new_type == RT::ROTATION_TYPE_NO_ROTATION)
        {
          trackball_angle = viewer->core().trackball_angle;
          orthographic = viewer->core().orthographic;
          viewer->core().trackball_angle = Eigen::Quaternionf::Identity();
          viewer->core().orthographic = true;
        }
        else if (viewer->core().rotation_type == RT::ROTATION_TYPE_NO_ROTATION)
        {
          viewer->core().trackball_angle = trackball_angle;
          viewer->core().orthographic = orthographic;
        }
        viewer->core().set_rotation_type(new_type);
      }
    }

    // Orthographic view
    ImGui::Checkbox("Orthographic view", &(viewer->core().orthographic));
    ImGui::PopItemWidth();
  }

  // Helper for setting viewport specific mesh options
  auto make_checkbox = [&](const char *label, unsigned int &option)
  {
    return ImGui::Checkbox(label,
      [&]() { return viewer->core().is_set(option); },
      [&](bool value) { return viewer->core().set(option, value); }
    );
  };

  // Draw options
  if (ImGui::CollapsingHeader("Draw Options", ImGuiTreeNodeFlags_DefaultOpen))
  {
    if (ImGui::Checkbox("Face-based", &(viewer->data().face_based)))
    {
      viewer->data().dirty = MeshGL::DIRTY_ALL;
    }
    make_checkbox("Show texture", viewer->data().show_texture);
    if (ImGui::Checkbox("Invert normals", &(viewer->data().invert_normals)))
    {
      viewer->data().dirty |= igl::opengl::MeshGL::DIRTY_NORMAL;
    }
    make_checkbox("Show overlay", viewer->data().show_overlay);
    make_checkbox("Show overlay depth", viewer->data().show_overlay_depth);
    ImGui::ColorEdit4("Background", viewer->core().background_color.data(),
        ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel);
    ImGui::ColorEdit4("Line color", viewer->data().line_color.data(),
        ImGuiColorEditFlags_NoInputs | ImGuiColorEditFlags_PickerHueWheel);
    ImGui::PushItemWidth(ImGui::GetWindowWidth() * 0.3f);
    ImGui::DragFloat("Shininess", &(viewer->data().shininess), 0.05f, 0.0f, 100.0f);
    ImGui::PopItemWidth();
  }

  // Overlays
  if (ImGui::CollapsingHeader("Overlays", ImGuiTreeNodeFlags_DefaultOpen))
  {
    make_checkbox("Wireframe", viewer->data().show_lines);
    make_checkbox("Fill", viewer->data().show_faces);
    ImGui::Checkbox("Show vertex labels", &(viewer->data().show_vertid));
    ImGui::Checkbox("Show faces labels", &(viewer->data().show_faceid));
    ImGui::Checkbox("Show extra labels", &(viewer->data().show_labels));
    if(viewer->data().show_vertid 
    || viewer->data().show_faceid
    || viewer->data().show_labels)
    {
      ImGui::Checkbox("Use depth testing", &(viewer->data().filter_labels));
    }
  }
}

IGL_INLINE void ImGuiMenu::draw_labels_window()
{
  // Text labels
  ImGui::SetNextWindowPos(ImVec2(0,0), ImGuiCond_Always);
  ImGui::SetNextWindowSize(ImGui::GetIO().DisplaySize, ImGuiCond_Always);
  bool visible = true;
  ImGui::PushStyleColor(ImGuiCol_WindowBg, ImVec4(0,0,0,0));
  ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0);
  ImGui::Begin("ViewerLabels", &visible,
      ImGuiWindowFlags_NoTitleBar
      | ImGuiWindowFlags_NoResize
      | ImGuiWindowFlags_NoMove
      | ImGuiWindowFlags_NoScrollbar
      | ImGuiWindowFlags_NoScrollWithMouse
      | ImGuiWindowFlags_NoCollapse
      | ImGuiWindowFlags_NoSavedSettings
      | ImGuiWindowFlags_NoInputs);
  for (const auto & data : viewer->data_list)
  {
    draw_labels(data);
  }
  ImGui::End();
  ImGui::PopStyleColor();
  ImGui::PopStyleVar();
}

IGL_INLINE void ImGuiMenu::draw_labels(const igl::opengl::ViewerData &data)
{

  // Perform label occlusion detection
  if(viewer->data().filter_labels)
  {
    filterLabelsByDepth();
  }

  // Alec: How can we get these to respect (optionally) the depth of the scene?
  if (data.show_vertid)
  {
    for (int i = 0; i < data.V.rows(); ++i)
    {
      int show_label  = data.filter_labels ? viewer->data().vertex_label_mask(i,0) : 1;
      if(show_label) {
        draw_text(
          data.V.row(i), 
          data.V_normals.row(i), 
          std::to_string(i),
          data.label_color);
      }
    }
  }

  if (data.show_faceid)
  {
    for (int i = 0; i < data.F.rows(); ++i)
    {
      int show_label  = data.filter_labels ? viewer->data().face_label_mask(i,0) : 1;
      if(show_label) {
        Eigen::RowVector3d p = Eigen::RowVector3d::Zero();
        for (int j = 0; j < data.F.cols(); ++j)
        {
          p += data.V.row(data.F(i,j));
        }
        p /= (double) data.F.cols();

        draw_text(
          p, 
          data.F_normals.row(i), 
          std::to_string(i),
          data.label_color);
      }
    }
  }

  if (data.show_labels)
  {
    for (int i = 0; i < data.labels_positions.rows(); ++i)
    {
      draw_text(
        data.labels_positions.row(i), 
        Eigen::Vector3d(0.0,0.0,0.0),
        data.labels_strings[i],
        data.label_color);
    }
  }
}

IGL_INLINE void ImGuiMenu::draw_text(
  Eigen::Vector3d pos, 
  Eigen::Vector3d normal, 
  const std::string &text,
  const Eigen::Vector4f color)
{
  pos += normal * 0.005f * viewer->core().object_scale;
  Eigen::Vector3f coord = igl::project(Eigen::Vector3f(pos.cast<float>()),
    viewer->core().view, viewer->core().proj, viewer->core().viewport);

  // Draw text labels slightly bigger than normal text
  ImDrawList* drawList = ImGui::GetWindowDrawList();
  drawList->AddText(ImGui::GetFont(), ImGui::GetFontSize() * 1.2,
      ImVec2(coord[0]/pixel_ratio_, (viewer->core().viewport[3] - coord[1])/pixel_ratio_),
      ImGui::GetColorU32(ImVec4(
        color(0),
        color(1),
        color(2),
        color(3))),
      &text[0], &text[0] + text.size());
}

IGL_INLINE float ImGuiMenu::pixel_ratio()
{
  // Computes pixel ratio for hidpi devices
  int buf_size[2];
  int win_size[2];
  GLFWwindow* window = glfwGetCurrentContext();
  glfwGetFramebufferSize(window, &buf_size[0], &buf_size[1]);
  glfwGetWindowSize(window, &win_size[0], &win_size[1]);
  return (float) buf_size[0] / (float) win_size[0];
}

IGL_INLINE float ImGuiMenu::hidpi_scaling()
{
  // Computes scaling factor for hidpi devices
  float xscale, yscale;
  GLFWwindow* window = glfwGetCurrentContext();
  glfwGetWindowContentScale(window, &xscale, &yscale);
  return 0.5 * (xscale + yscale);
}

} // end namespace
} // end namespace
} // end namespace
} // end namespace
