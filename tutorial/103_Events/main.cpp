#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/bounding_box_diagonal.h>
#include <iostream>
#include "tutorial_shared_path.h"

Eigen::MatrixXd V;
Eigen::MatrixXi F;

// This function is called every time a keyboard button is pressed
bool key_down(igl::opengl::glfw::Viewer& viewer, unsigned char key, int modifier)
{
  std::cout<<"Key: "<<key<<" "<<(unsigned int)key<<std::endl;
  if (key == '1')
  {
    // Clear should be called before drawing the mesh
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
    viewer.core.radius_in_screen_space = false;
    Eigen::VectorXd radius(V.rows());
    radius.setConstant(0.005 / igl::bounding_box_diagonal(V));
    viewer.data().set_points(V, Eigen::RowVector3d(1, 0, 0), radius);
    viewer.core.align_camera_center(V,F);
  }
  else if (key == '2')
  {
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
    viewer.core.radius_in_screen_space = true;
    Eigen::VectorXd radius(V.rows());
    radius.setConstant(0.5);
    viewer.data().set_points(V, Eigen::RowVector3d(1, 0, 0), radius);
    viewer.core.align_camera_center(V,F);
  }

  return false;
}


int main(int argc, char *argv[])
{
  // Load two meshes
  igl::readOFF(TUTORIAL_SHARED_PATH "/bunny.off", V, F);
  igl::opengl::glfw::Viewer viewer;
  // Register a keyboard callback that allows to switch between
  // the two loaded meshes
  viewer.callback_key_down = &key_down;
  viewer.data().set_mesh(V, F);
  Eigen::VectorXd radius(V.rows());
  radius.setConstant(0.005 / igl::bounding_box_diagonal(V));
  viewer.data().set_points(V, Eigen::RowVector3d(1, 0, 0), radius);
  viewer.launch();
}
