#include "tutorial_shared_path.h"
#include <igl/opengl/glfw/Viewer.h>
#include <igl/opengl/glfw/imgui/ImGuiMenu.h>
#include <igl/opengl/glfw/imgui/ImGuiHelpers.h>
#include <string>
#include <map>
#include <iostream>

std::vector<std::string> names = {"cube.obj","sphere.obj","xcylinder.obj","ycylinder.obj","zcylinder.obj"};
std::map<std::string, Eigen::RowVector3d> colors;
int last_colored_index = -1;

// Set colors of each mesh
void update_colors(igl::opengl::glfw::Viewer &viewer)
{
  for (int i = 0; i < (int) viewer.data_list.size(); ++i)
  {
    viewer.data_list[i].set_colors(colors[names[i]]);
  }
  viewer.data_list[viewer.selected_data_index].set_colors(Eigen::RowVector3d(0.9,0.1,0.1));
  last_colored_index = viewer.selected_data_index;
}

// Menu for selecting a mesh
class MyMenu : public igl::opengl::glfw::imgui::ImGuiMenu
{
  virtual void draw_viewer_menu() override
  {
    if (ImGui::Combo("Selected Mesh", (int *) &viewer->selected_data_index, names))
    {
      update_colors(*viewer);
    }
    if (last_colored_index != viewer->selected_data_index)
    {
      update_colors(*viewer);
    }
  }
};

int main(int argc, char * argv[])
{
  igl::opengl::glfw::Viewer viewer;
  for(const auto & name : names)
  {
    viewer.load_mesh_from_file(std::string(TUTORIAL_SHARED_PATH) + "/" + name);
    colors[name] = Eigen::RowVector3d::Random();
  }

  // Attach a custom menu
  MyMenu menu;
  viewer.plugins.push_back(&menu);

  // Color each mesh differently
  update_colors(viewer);

  viewer.callback_key_down =
    [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
  {
    // Delete
    if(key == '3')
    {
      int old_index = viewer.selected_data_index;
      if (viewer.erase_mesh(viewer.selected_data_index))
      {
        names.erase(names.begin() + old_index);
        update_colors(viewer);
      }
      return true;
    }
    return false;
  };

  viewer.launch();
  return EXIT_SUCCESS;
}
