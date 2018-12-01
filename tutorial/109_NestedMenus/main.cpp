#include <igl/read_triangle_mesh.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <imgui/imgui.h>
#include <iostream>
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

  viewer.callback_key_down = [&](igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier) {
    std::cout<<"Key: "<<key<<" "<<(unsigned int)key<<std::endl;
    if (key == 'V')
    {
      igl::opengl::glfw::Viewer nested;

      // Attach a menu plugin
      // igl::opengl::glfw::imgui::ImGuiMenu menu;
      // viewer.plugins.push_back(&menu);

      nested.data().set_mesh(V2, F2);
      nested.data().set_face_based(true);
      nested.launch();
      return true;
    }
    return false;
  };

  // Attach a menu plugin
  // igl::opengl::glfw::imgui::ImGuiMenu menu;
  // viewer.plugins.push_back(&menu);

  // Plot the mesh
  viewer.data().set_mesh(V1, F1);
  viewer.launch();
}
