#include "tutorial_shared_path.h"
#include <igl/opengl/glfw/Viewer.h>
#include <string>

int main(int argc, char * argv[])
{
  igl::opengl::glfw::Viewer v;
  for(const auto & name : 
    {"cube.obj","sphere.obj","xcylinder.obj","ycylinder.obj","zcylinder.obj"})
  {
    v.load_mesh_from_file(std::string(TUTORIAL_SHARED_PATH) + "/" + name);
  }
  v.selected_data_index = 0;
  v.selected_data().set_colors(Eigen::RowVector3d(0.8,0.47,0.22));
  v.selected_data_index = 1;
  v.selected_data().set_colors(Eigen::RowVector3d(0.6,0.01,0.11));
  v.selected_data_index = 2;
  v.selected_data().set_colors(Eigen::RowVector3d(0.37,0.06,0.25));
  v.selected_data_index = 3;
  v.selected_data().set_colors(Eigen::RowVector3d(1,1,1));
  v.launch();
  return EXIT_SUCCESS;
}
