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
  
  int leftView, rightView;
  int cubeID = viewer.data_list[0].id, sphereID = viewer.data_list[1].id;
  viewer.callback_init = 
	  [&](igl::opengl::glfw::Viewer &)
  {
    viewer.core().viewport = Eigen::Vector4f(0, 0, 640, 800);
	leftView = viewer.core_list[0].id;
	
	rightView = viewer.append_core(Eigen::Vector4f(640, 0, 640, 800));

	return true;
  };

  viewer.callback_key_down =
    [&](igl::opengl::glfw::Viewer &, unsigned int key, int mod)
  {
    if(key == GLFW_KEY_SPACE)
    {
		// by default, when a core is appended, all loaded meshes will be displayed in that core
		// displaying can be controlled by changing viewer.coreDataPairs
		viewer.data_list[viewer.mesh_index(cubeID)].is_visible = 0;
		viewer.data_list[viewer.mesh_index(cubeID)].set_visible(true, viewer.core_index(leftView));
		viewer.data_list[viewer.mesh_index(sphereID)].is_visible = 0;
		viewer.data_list[viewer.mesh_index(sphereID)].set_visible(true, viewer.core_index(rightView));
    }
    return false;
  };


  viewer.launch();
  return EXIT_SUCCESS;
}
