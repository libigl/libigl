#include <igl/kelvinlets.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/readOFF.h>
#include <igl/unproject.h>
#include <igl/unproject_onto_mesh.h>
#include <imgui/imgui.h>
#include <iostream>

namespace {
void ShowHelpMarker(const char* desc)
{
  ImGui::SameLine();
  ImGui::TextDisabled("(?)");
  if (ImGui::IsItemHovered()) {
    ImGui::BeginTooltip();
    ImGui::PushTextWrapPos(450.0f);
    ImGui::TextUnformatted(desc);
    ImGui::PopTextWrapPos();
    ImGui::EndTooltip();
  }
}
}

int main()
{
  Eigen::MatrixXd V1, OrigV;
  Eigen::MatrixXi F1, OrigF;

  igl::readOFF(TUTORIAL_SHARED_PATH "/bumpy.off", OrigV, OrigF);
  std::cout << "1 View original mesh\n";
  std::cout << "2 Switch to deformed mesh\n";
  V1 = OrigV;
  F1 = OrigF;

  igl::opengl::glfw::Viewer viewer;
  igl::opengl::glfw::imgui::ImGuiMenu menu;

  viewer.plugins.push_back(&menu);

  auto brushRadius = 1.;
  auto brushType = igl::BrushType::GRAB;
  auto scale = 1;
  menu.callback_draw_custom_window = [&]() {
    ImGui::SetNextWindowPos(ImVec2(180.f * menu.menu_scaling(), 10),
                            ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(200, 160), ImGuiCond_FirstUseEver);
    ImGui::Begin(
      "Kelvinlet Brushes", nullptr, ImGuiWindowFlags_NoSavedSettings);

    ImGui::InputDouble("Brush Radius", &brushRadius, 0, 0, "%.4f");
    ImGui::Combo("Brush type",
                 reinterpret_cast<int*>(&brushType),
                 "Grab\0Scale\0Twist\0Pinch\0\0");
    ImGui::InputInt("Falloff", &scale);
    ShowHelpMarker("Defines how localized the stroke is {1,2,3}");
    ImGui::End();
  };

  Eigen::Vector3d posStart(0, 0, 0);
  Eigen::Vector3d posEnd;
  decltype(OrigV) result;
  auto min_point = V1.colwise().minCoeff();
  auto max_point = V1.colwise().maxCoeff();
  // to multiply brush force proportional to size of mesh
  auto brush_strength = (max_point - min_point).norm();
  Eigen::Matrix3d twist, pinch;
  twist << 0, 1, -1, -1, 0, 1, 1, -1, 0; // skew-symmetric
  pinch << 0, 1, 1, 1, 0, 1, 1, 1, 0;    // symmetric

  viewer.callback_key_down =
    [&](igl::opengl::glfw::Viewer& viewer, unsigned char key, int) {
      std::cout << "Key: " << key << " " << (unsigned int)key << std::endl;
      if (key == '1') {
        viewer.data().clear();
        viewer.data().set_mesh(OrigV, OrigF);
        viewer.core().align_camera_center(OrigV, OrigF);
      } else if (key == '2') {
        viewer.data().clear();
        viewer.data().set_mesh(V1, F1);
        viewer.core().align_camera_center(V1, F1);
      }
      return false;
    };

  viewer.callback_mouse_down =
    [&](igl::opengl::glfw::Viewer& viewer, int, int) -> bool {
    Eigen::Vector3f bc;
    int fid;
    auto x = viewer.current_mouse_x;
    auto y =
      viewer.core().viewport(3) - static_cast<float>(viewer.current_mouse_y);
    if (igl::unproject_onto_mesh(Eigen::Vector2f(x, y),
                                 viewer.core().view,
                                 viewer.core().proj,
                                 viewer.core().viewport,
                                 V1,
                                 F1,
                                 fid,
                                 bc)) {
      posStart = igl::unproject(Eigen::Vector3f(x, y, viewer.down_mouse_z),
                                viewer.core().view,
                                viewer.core().proj,
                                viewer.core().viewport)
                   .template cast<double>();
      return true;
    }
    return false;
  };

  viewer.callback_mouse_move =
    [&](igl::opengl::glfw::Viewer& viewer, int, int) -> bool {
    if (!posStart.isZero() && !posStart.hasNaN()) {
      posEnd = igl::unproject(
                 Eigen::Vector3f(viewer.current_mouse_x,
                                 viewer.core().viewport[3] -
                                   static_cast<float>(viewer.current_mouse_y),
                                 viewer.down_mouse_z),
                 viewer.core().view,
                 viewer.core().proj,
                 viewer.core().viewport)
                 .template cast<double>();

      // exaggerate the force by a little bit
      Eigen::Vector3d forceVec = (posEnd - posStart) * brush_strength;

      int scaleFactor = forceVec.norm();
      if (posEnd.x() < posStart.x()) {
        // probably not the best way to determine direction.
        scaleFactor = -scaleFactor;
      }
      Eigen::Matrix3d mat;
      switch (brushType) {
        case igl::BrushType::GRAB:
          mat.setZero();
          break;
        case igl::BrushType::SCALE:
          mat = Eigen::Matrix3d::Identity() * scaleFactor;
          break;
        case igl::BrushType::TWIST:
          mat = twist * scaleFactor;
          break;
        case igl::BrushType::PINCH:
          mat = pinch * scaleFactor;
          break;
      }

      igl::kelvinlets(
        V1,
        posStart,
        forceVec,
        mat,
        igl::KelvinletParams<double>(brushRadius, scale, brushType),
        result);
      viewer.data().set_mesh(result, F1);
      viewer.core().align_camera_center(result, F1);
      return true;
    }
    return false;
  };

  viewer.callback_mouse_up =
    [&](igl::opengl::glfw::Viewer& viewer, int, int) -> bool {
    if (!posStart.isZero()) {
      V1 = result;
      posStart.setZero();
      return true;
    }
    return false;
  };

  viewer.data().set_mesh(V1, F1);
  viewer.core().align_camera_center(V1, F1);
  viewer.launch();
}
