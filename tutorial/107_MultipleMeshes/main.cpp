#include "tutorial_shared_path.h"
#include <igl/opengl/glfw/Viewer.h>
#include <GLFW/glfw3.h>
#include <string>
#include <iostream>

int main(int argc, char * argv[])
{
  igl::opengl::glfw::Viewer viewer;
  const auto names = 
    {"cube.obj","sphere.obj","xcylinder.obj","ycylinder.obj","zcylinder.obj"};
  for(const auto & name : names)
  {
    viewer.load_mesh_from_file(std::string(TUTORIAL_SHARED_PATH) + "/" + name);
  }

  // Set colors of each mesh by selecting its index first
  viewer.selected_data_index = 0;
  viewer.data().set_colors(Eigen::RowVector3d(0.8,0.47,0.22));
  viewer.selected_data_index = 1;
  viewer.data().set_colors(Eigen::RowVector3d(0.6,0.01,0.11));
  viewer.selected_data_index = 2;
  viewer.data().set_colors(Eigen::RowVector3d(0.37,0.06,0.25));
  viewer.selected_data_index = 3;
  viewer.data().set_colors(Eigen::RowVector3d(1,1,1));

  viewer.callback_key_down = 
    [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
  {
    if(key == GLFW_KEY_BACKSPACE)
    {
      viewer.erase_mesh(viewer.selected_data_index);
      return true;
    }
    return false;
  };

  viewer.launch();
  return EXIT_SUCCESS;
}
