#include <igl/readOFF.h>
#include <igl/opengl/glfw/Viewer.h>
#include <igl/bounding_box_diagonal.h>
#include <igl/avg_edge_length.h>
#include <igl/read_triangle_mesh.h>
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
    viewer.core.align_camera_center(V,F);
    viewer.core.radius_in_screen_space = false;
    Eigen::VectorXd radius(V.rows());
    radius.setConstant(igl::avg_edge_length(V, F) / 4.0 * viewer.core.camera_base_zoom);
    viewer.data().set_points(V, Eigen::RowVector3d(1, 0, 0), radius);
  }
  else if (key == '2')
  {
    viewer.data().clear();
    viewer.data().set_mesh(V, F);
    viewer.core.align_camera_center(V,F);
    viewer.core.radius_in_screen_space = true;
    Eigen::VectorXd radius(V.rows());
    radius.setConstant(0.05 * viewer.core.camera_base_zoom);
    viewer.data().set_points(V, Eigen::RowVector3d(1, 0, 0), radius);
  }

  return false;
}


int main(int argc, char *argv[])
{
  // Load two meshes
  igl::read_triangle_mesh(TUTORIAL_SHARED_PATH "/bumpy.off", V, F);
  igl::opengl::glfw::Viewer viewer;
  // Register a keyboard callback that allows to switch
  viewer.callback_key_down = &key_down;
  viewer.data().set_mesh(V, F);
  viewer.launch();
}
