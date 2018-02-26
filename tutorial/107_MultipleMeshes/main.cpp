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
  std::vector<Eigen::RowVector3d> colors;
  int last_selected = -1;
  for(const auto & name : names)
  {
    viewer.load_mesh_from_file(std::string(TUTORIAL_SHARED_PATH) + "/" + name);
    colors.push_back(0.5*Eigen::RowVector3d::Random().array() + 0.5);
  }

  viewer.callback_key_down =
    [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
  {
    if(key == GLFW_KEY_BACKSPACE)
    {
      int old_index = viewer.selected_data_index;
      if (viewer.erase_mesh(viewer.selected_data_index))
      {
        colors.erase(colors.begin() + old_index);
        last_selected = -1;
      }
      return true;
    }
    return false;
  };

  // Refresh selected mesh colors
  viewer.callback_pre_draw =
    [&](igl::opengl::glfw::Viewer &)
  {
    if (last_selected != viewer.selected_data_index)
    {
      for (size_t i = 0; i < viewer.data_list.size(); ++i)
      {
        viewer.data_list[i].set_colors(colors[i]);
      }
      viewer.data_list[viewer.selected_data_index].set_colors(Eigen::RowVector3d(0.9,0.1,0.1));
      last_selected = viewer.selected_data_index;
    }
    return false;
  };

  viewer.launch();
  return EXIT_SUCCESS;
}
