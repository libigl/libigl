#include "tutorial_shared_path.h"
#include <igl/opengl/glfw/Viewer.h>
#include <GLFW/glfw3.h>
#include <string>
#include <iostream>
#include <map>

int main(int argc, char * argv[])
{
  igl::opengl::glfw::Viewer viewer;
  const auto names =
    {"cube.obj","sphere.obj","xcylinder.obj","ycylinder.obj","zcylinder.obj"};
  std::map<int, Eigen::RowVector3d> colors;
  int last_selected = -1;
  for(const auto & name : names)
  {
    viewer.load_mesh_from_file(std::string(TUTORIAL_SHARED_PATH) + "/" + name);
    colors.emplace(viewer.data().id, 0.5*Eigen::RowVector3d::Random().array() + 0.5);
  }

  viewer.callback_key_down =
    [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
  {
    if(key == GLFW_KEY_BACKSPACE)
    {
      int old_id = viewer.data().id;
      if (viewer.erase_mesh(viewer.selected_data_index))
      {
        colors.erase(old_id);
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
      for (auto &data : viewer.data_list)
      {
        data.set_colors(colors[data.id]);
      }
      viewer.data_list[viewer.selected_data_index].set_colors(Eigen::RowVector3d(0.9,0.1,0.1));
      last_selected = viewer.selected_data_index;
    }
    return false;
  };

  viewer.launch();
  return EXIT_SUCCESS;
}
