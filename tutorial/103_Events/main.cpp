#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <iostream>
#include "tutorial_shared_path.h"

Eigen::MatrixXd V1,V2;
Eigen::MatrixXi F1,F2;

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  std::cout<<"Key: "<<key<<" "<<(unsigned int)key<<std::endl;
  if (key == '1')
  {
    // Clear should be called before drawing the mesh
    viewer.data().clear();
    // Draw_mesh creates or updates the vertices and faces of the displayed mesh.
    // If a mesh is already displayed, draw_mesh returns an error if the given V and
    // F have size different than the current ones
    viewer.data().set_mesh(V1, F1);
    viewer.core().align_camera_center(V1,F1);
  }
  else if (key == '2')
  {
    viewer.data().clear();
    viewer.data().set_mesh(V2, F2);
    viewer.core().align_camera_center(V2,F2);
  }

  return false;
}


int main(int argc, char *argv[])
{
  // Load two meshes
  igl::readOFF(TUTORIAL_SHARED_PATH "/bumpy.off", V1, F1);
  igl::readOFF(TUTORIAL_SHARED_PATH "/fertility.off", V2, F2);
  std::cout<<R"(
1 Switch to bump mesh
2 Switch to fertility mesh
    )";

  igl::opengl::glfw::Viewer viewer;
  // Register a keyboard callback that allows to switch between
  // the two loaded meshes
  viewer.callback_key_down = &key_down;
  viewer.data().set_mesh(V1, F1);
  viewer.launch();
}
