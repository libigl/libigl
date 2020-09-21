#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuizmo.h>
#include <imgui_impl_glfw.h>
#include <imgui_impl_opengl3.h>
#include <GLFW/glfw3.h>
#include <imgui/imgui.h>
#include "tutorial_shared_path.h"

class SlicingPlugin : public igl::opengl::glfw::imgui::ImGuiMenu
{
  igl::opengl::ViewerData data;

  const Eigen::MatrixXd OV = (Eigen::MatrixXd(4, 3) <<
    -0.5, -0.5, 0.0,
    -0.5,  0.5, 0.0,
    0.5,  0.5, 0.0,
    0.5, -0.5, 0.0).finished();

  const Eigen::MatrixXi OF = (Eigen::MatrixXi(2, 3) <<
    0, 2, 1,
    0, 3, 2).finished();

  virtual void init(igl::opengl::glfw::Viewer *_viewer) override {
    ImGuiMenu::init(_viewer);

    // Load slicing plane into viewer (a square mesh)
    data.set_mesh(SlicingPlugin::OV, SlicingPlugin::OF);
    data.set_face_based(true);
    data.set_colors(Eigen::RowVector4d(224, 86, 253, 128)/255.0);
    data.show_lines = false;
  }

  virtual bool pre_draw() override {
    ImGuiMenu::pre_draw();
    ImGui::PushStyleVar(ImGuiStyleVar_WindowBorderSize, 0);
    ImGuizmo::BeginFrame();
    ImGui::PopStyleVar();
    return false;
  }

  virtual bool post_draw() override {
    viewer->core().draw(data);
    ImGuiMenu::post_draw();
    return false;
  }

  virtual void draw_custom_window() override {
    float menu_width = 180.f * menu_scaling();
    ImGui::SetNextWindowPos(ImVec2(menu_width, 0.0f), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSize(ImVec2(0.0f, 0.0f), ImGuiCond_FirstUseEver);
    ImGui::SetNextWindowSizeConstraints(ImVec2(1.5f * menu_width, -1.0f), ImVec2(1.5f * menu_width, -1.0f));
    bool _guizmo_menu_visible = true;
    ImGui::Begin(
      "ImGuizmo Tools", &_guizmo_menu_visible,
      ImGuiWindowFlags_NoSavedSettings
      | ImGuiWindowFlags_AlwaysAutoResize
    );

    draw_imguizmo_menu();

    ImGui::End();
  }

  virtual void draw_imguizmo_menu() {
    static Eigen::Matrix4f matrix = Eigen::Matrix4f::Identity();
    Eigen::Affine3f rescale = Eigen::Scaling(0.5f * viewer->core().camera_base_zoom)
      * Eigen::Translation3f(viewer->core().camera_base_translation);
    Eigen::Affine3f view = Eigen::Affine3f( viewer->core().view * 1./viewer->core().camera_zoom ) * rescale.inverse();
    Eigen::Matrix4f proj = viewer->core().proj;
    if(viewer->core().orthographic)
    {
      view(2,3) -= 50; // Adjust depth for view transform
    }

    igl::opengl::glfw::imgui::EditTransform(view.matrix().data(), proj.data(), matrix.data(), viewer->core().orthographic);

    // Transform the slicing plane according to
    // ImGuizmo tool manipulations in the viewer
    Eigen::Affine3f model = Eigen::Affine3f(matrix) * rescale.inverse();
    Eigen::MatrixXd V = (SlicingPlugin::OV.rowwise().homogeneous()
      * model.matrix().cast<double>().transpose()).rowwise().hnormalized();

    if (data.V.rows() == V.rows())
      data.set_vertices(V);
  }

};

int main(int argc, char *argv[])
{
  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  std::string filename(TUTORIAL_SHARED_PATH "/cow.off");
  viewer.load_mesh_from_file(filename);

  // Custom menu
  SlicingPlugin menu;
  viewer.plugins.push_back(&menu);

  viewer.launch();
}
