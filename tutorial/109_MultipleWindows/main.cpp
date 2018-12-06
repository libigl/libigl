#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <iostream>
#include <list>
#include "tutorial_shared_path.h"

int main(int argc, char *argv[])
{
  Eigen::MatrixXd V1, V2;
  Eigen::MatrixXi F1, F2;

  // Load a mesh in OFF format
  igl::read_triangle_mesh(TUTORIAL_SHARED_PATH "/bunny.off", V1, F1);
  igl::read_triangle_mesh(TUTORIAL_SHARED_PATH "/cube.obj", V2, F2);

  // Init the viewer
  igl::opengl::glfw::Viewer viewer;
  std::list<igl::opengl::glfw::Viewer> nested_viewers;
  std::list<igl::opengl::glfw::imgui::ImGuiMenu> nested_menus;

  viewer.callback_key_down = [&](igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
  {
    std::cout<<"Key: "<<key<<" "<<(unsigned int)key<<std::endl;
    if (key == '1')
    {
      nested_viewers.emplace_back();

      // Attach a menu plugin
      nested_menus.emplace_back();
      // nested_menus.back().callback_draw_viewer_menu = [&]()
      // {
      //   ImGui::Button("Yay");
      // };
      nested_viewers.back().plugins.push_back(&nested_menus.back());

      nested_viewers.back().data().set_mesh(V2, F2);
      nested_viewers.back().data().set_face_based(true);
      nested_viewers.back().launch_with(&viewer);
      return true;
    }
    if (key == '2')
    {
      igl::opengl::glfw::Viewer subviewer;

      // Attach a menu plugin
      igl::opengl::glfw::imgui::ImGuiMenu submenu;
      subviewer.plugins.push_back(&submenu);

      subviewer.data().set_mesh(V2, F2);
      subviewer.data().set_face_based(true);
      subviewer.launch();
      return true;
    }
    return false;
  };

  // Attach a menu plugin
  igl::opengl::glfw::imgui::ImGuiMenu menu;
  viewer.plugins.push_back(&menu);

  // Plot the mesh
  viewer.data().set_mesh(V1, F1);
  // viewer.core.is_animating = true;
  viewer.launch();
}
