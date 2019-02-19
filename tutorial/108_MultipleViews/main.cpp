#include "tutorial_shared_path.h"
#include <igl/opengl/glfw/Viewer.h>
#include <GLFW/glfw3.h>
#include <string>
#include <iostream>
#include <map>

int main(int argc, char * argv[])
{
  igl::opengl::glfw::Viewer viewer;

  viewer.load_mesh_from_file(std::string(TUTORIAL_SHARED_PATH) + "/cube.obj");
  viewer.load_mesh_from_file(std::string(TUTORIAL_SHARED_PATH) + "/sphere.obj");

  unsigned int left_view, right_view;
  int cube_id = viewer.data_list[0].id;
  int sphere_id = viewer.data_list[1].id;
  viewer.callback_init = [&](igl::opengl::glfw::Viewer &)
  {
    viewer.core().viewport = Eigen::Vector4f(0, 0, 640, 800);
    left_view = viewer.core_list[0].id;
    right_view = viewer.append_core(Eigen::Vector4f(640, 0, 640, 800));
    return false;
  };

  viewer.callback_key_down = [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
  {
    if(key == GLFW_KEY_SPACE)
    {
      // By default, when a core is appended, all loaded meshes will be displayed in that core.
      // Displaying can be controlled by calling viewer.data().set_visible().
      viewer.data(cube_id).set_visible(false, left_view);
      viewer.data(sphere_id).set_visible(false, right_view);
    }
    return false;
  };

  viewer.callback_post_resize = [&](igl::opengl::glfw::Viewer &v, int w, int h) {
    v.core( left_view).viewport = Eigen::Vector4f(0, 0, w / 2, h);
    v.core(right_view).viewport = Eigen::Vector4f(w / 2, 0, w - (w / 2), h);
    return true;
  };

  viewer.launch();
  return EXIT_SUCCESS;
}
