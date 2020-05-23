#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <imgui/imgui.h>
#include "slicing_plugin.h"
#include <igl/readSTL.h>

int main(int argc, char *argv[])
{

  Eigen::MatrixXd V;
  Eigen::MatrixXi F;
  Eigen::MatrixXd N;

  // Plot the mesh
  igl::opengl::glfw::Viewer viewer;
  std::string path = "/home/michelle/Documents/LIBIGL/LABELS_BUG/mishfork/libigl-example-project/3rdparty/libigl/tutorial/data/";
  path += "cow.off";
  viewer.load_mesh_from_file(path);
  std::cout << "The size of unsigned int: " << sizeof(unsigned int) << std::endl;

  // Custom menu
  SlicingPlugin menu;
  viewer.plugins.push_back(&menu);

  viewer.resize(1024, 1024);
  viewer.launch();
}